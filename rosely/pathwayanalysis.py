import pandas as pd
from scipy.special import betainc
from collections import OrderedDict
from .neutralstats import locfdr
from .libanalysis import freqcdf
import requests, os, re, time
import numpy as np
import xml.etree.cElementTree as et

try:
    from Bio.KEGG import REST
except: pass

__all__ = ['enrichment_p_value', 'download_all_kegg_pathways', 'enrichment_analysis', 'draw_kegg_pathways',
           'load_all_kegg_pathways', 'getGeneId', 'drawPathway', 'downloadMyGeneInfo']

def download_all_kegg_pathways(species_code='mmu'):
    """
    
    """
    pathways_str = REST.kegg_list("pathway", species_code).read()
    pathways = {p.split('\t')[0]:{'name':p.split('\t')[1]} for p in pathways_str.rstrip().split('\n')}
    def get_genes_for(pathways):
        for pathway in pathways:
            pathways[pathway]['geneid'] = set(); pathways[pathway]['gene_symbol'] = set()
            pathway_file = REST.kegg_get(pathway).read()  # query and read each pathway
            # iterate through each KEGG pathway file, keeping track of which section
            # of the file we're in, only read the gene in each pathway
            current_section = None
            for line in pathway_file.rstrip().split("\n"):
                section = line[:12].strip()  # section names are within 12 columns
                if not section == "":
                    current_section = section
                if current_section == "GENE":
                    try:
                        gene_identifiers, _ = line[12:].split("; ")[:2]
                        geneid, gene_symbol = gene_identifiers.split()
                        pathways[pathway]['geneid'].add(int(geneid))
                        pathways[pathway]['gene_symbol'].add(gene_symbol)
                    except: pass#print('Discarded:', line); 
    
    get_genes_for(pathways)
    return pathways

def load_all_kegg_pathways(species_code='mmu'):
    import pickle
    fn = 'pathways_' + species_code + '.pkl'
    if os.path.isfile(fn):
        with open(fn, 'rb') as fh: pathways = pickle.load(fh)
    else:
        pathways = download_all_kegg_pathways(species_code='mmu')
        with open(fn, 'wb') as fh: pickle.dump(pathways, fh)
    return pathways

def enrichment_p_value(Ltop, n, Lpop, N):
    """
    n out of Ltop genes in top list, while N out of Lpop genes in background (population).
    Assume that Lpop is large and N/Lpop is an accurate estimate of probability for occurence.
    
    ss.binom.cdf(X, N, p) == betainc(N-X, X+1, 1-p)
    """
    p = 1 - betainc(Ltop-(n+0.25-1), n+0.25, 1 - N/Lpop) # 0.25 is added to avoid p > 0.5
    p = min(1, 2 * p)
    return p

def enrichment_analysis(top_genes, population, pathways, tosymbols=None, values=None, valuename='', key='geneid', nbins=30):
    Ltop = len(top_genes); Lpop = len(population)
    ps = {}; ns = {}; Ns = {}; enpath = {}; genes = {}; fold = {}
    for pw in pathways:
        symbol = set(pathways[pw][key])
        gs = list(symbol.intersection(top_genes))
        n = len(gs)
        N = len(symbol.intersection(population))
        if N > 1 and n / (1+Ltop) >= N / Lpop:
            ps[pw] = freqcdf(Ltop, n, Lpop, N, 1)[0]
            if not 0<=ps[pw]<=1: print('p value error:', ps[pw], Ltop, n, Lpop, N)
            ns[pw] = str(n)+' / '+str(Ltop)
            Ns[pw] = str(N)+' / '+str(Lpop)
            enpath[pw] = pathways[pw]['name'].split(' - ')[0]
            if tosymbols is not None:
                gs0 = gs
                if values is not None: vs = [values[g] for g in gs]
                gs = [tosymbols[g] for g in gs]
            elif key == 'geneid': gs = [str(int(g)) for g in gs]
            if values is not None:
                try: genes[pw] = ', '.join([g + ' ' + str(round(v, 2)) for (g, v) in zip(gs, vs)])
                except: print(gs0); print(gs); print(vs); raise
            else:
                genes[pw] = ', '.join(gs) 
            fold[pw] = round((n/Ltop) / (N/Lpop), 2)
    ps = pd.Series(ps)
    ps = (ps / ps.max()).to_dict()
    fdr = locfdr(ps, nbins=nbins, p0s=[0,])
    pathname = 'Pathways'
    if 'GO:' == pw[:3]: pathname = 'GO terms'
    if valuename != '': valuename = ' & ' + valuename
    enriched = pd.DataFrame(OrderedDict({pathname:enpath, 'p value':ps, 'LFDR':fdr, 
                                   'study':ns, 'background':Ns, 'enrich fold':fold,
                                   key+valuename:genes})).sort_values(by='p value')

    return enriched

def downloadMyGeneInfo(genes, scopes, species):
    import sys, io, mygene
    mg = mygene.MyGeneInfo()
    save_stdout = sys.stdout
    res = []
    print('Querying', scopes, end=' ')
    for i in range(0, len(genes), 1000):
        print(i+1, '~', min(i+1000, len(genes)), end=' ')
        while True:
            try: 
                sys.stdout = io.StringIO()
                r = mg.querymany(genes[i:i+1000], scopes=scopes, species=species)
                res.extend(r)
                sys.stdout = save_stdout
                break
            except KeyboardInterrupt: 
                sys.stdout = save_stdout; raise
            except: sys.stdout = save_stdout; time.sleep(1)
    print('')
    return res

def getGeneId(genes, scopes='ensemblgene,symbol', species='mouse', taxid=None):
    """
    Find out the Entrez gene id for genes from MyGeneInfo, which can be a list of ensemblgene id or gene symbols
    scoopes: set to 'symbol' if genes are gene symbols 
    
    return a pandas DataFrame of Entrez gene ids and gene symbols index by genes.
    """
    taxids = {'mouse':10090, 'human':9606}
    if taxid is None: taxid = taxids[species]
    idmap = {}; gs = list(genes)
    corrected = False
    for i in range(len(gs)):
        g = gs[i]; newid = g
        if len(g) > 4 and g[-4:].lower() == '-mar': newid = 'March' + g[:-4]
        if len(g) > 4 and g[-4:].lower() == '-sep': 
            newid = 'Sept'  + g[:-4]
            if newid == 'Sept15': newid = 'Sep15'
        if g != newid:
            if not corrected:  print('Symbol corrections: ', end='')
            print(g, '-->', newid, end='; ')
            corrected = True
        idmap[newid] = gs[i]
        gs[i] = newid
    if corrected: print('')
    raw = downloadMyGeneInfo(gs, scopes=scopes, species=taxid)
    for r in raw: 
        try: r['query'] = idmap[r['query']]
        except:
            for m in idmap:
                if r['query'] in m.split(','): r['query'] = idmap[m]
                if r['query'] in m.split(';'): r['query'] = idmap[m]
    ids = pd.DataFrame(raw)
    ids = ids[~ids['entrezgene'].isnull()]
    ids = ids.loc[ids.groupby('query')['_score'].idxmax()]
    ids = ids.set_index('query')
    df = pd.DataFrame()
    df['geneid'] = ids['entrezgene']
    df['symbol'] = ids['symbol']
    
    gs = list(df[df['geneid'].isnull()]['symbol'])
    if len(gs) and 'ensembl' in scopes:
        ids = pd.DataFrame(downloadMyGeneInfo(gs, scopes=scopes)).set_index('query')
        ids = ids[~ids['entrezgene'].isnull()][['entrezgene', 'symbol']]
        df.merge(ids, on='symbol')
    df = df[~df['geneid'].isnull()]
    df['geneid'] = df['geneid'].astype('int')
    return df

def draw_kegg_pathways(pathways, DEG_results, colorcolumn='Controlled z-score', folder=None, overlap_cutoff=0):
    """
    pathways: a pandas DataFrame indexed by path:id and ordered by decreasing statistical significance, with column 'Pathways'.
    DEG_results: a pandas DataFrame indexed by entrezgene and has a column `symbol` and a column for color values.
    overlap_cutoff: a cutoff value for displaying the gene symbols of genes co-localized in the same box of KEGG pathways.
    """
    if folder is not None and not os.path.isdir(folder): os.mkdir(folder)
    for (i, path) in enumerate(pathways.index):
        drawPathway(DEG_results, path, colorcolumn, cutoff=overlap_cutoff,
                    filename = os.path.join(folder, ("%02d. " % i) + pathways['Pathways'][path] + '_' + path.split(':')[-1]))

BLUE = np.array([0, 0, 255]); RED = np.array([255, 0, 0]); WHITE = np.array([255, 255, 255])
def value2color(v, vmax, pwr):
#     print(value, maxvalue)
    v = min(1, max(-1, v/vmax))
    v = np.sign(v) * abs(v) ** pwr
    if v > 0:
        color =  RED *  v   + WHITE * (1 - v)
    else:
        color = BLUE * (-v) + WHITE * (1 + v)
    return int(color[0]+0.5), int(color[1]+0.5), int(color[2]+0.5)

def drawPathway(data, pathwayId, colorcolumn = 'color value', cutoff=0, power = 3, filename=None):
    from PIL import Image, ImageDraw, ImageFont
    keggpath = os.path.join('kegg', pathwayId)
    if not os.path.isdir('kegg'): os.mkdir('kegg')
    if not os.path.isfile(keggpath + ".xml"):
        kgml = requests.get('http://rest.kegg.jp/get/'+pathwayId+'/kgml').text
        f = open(keggpath+'.xml','w')
        f.write(kgml)
        f.close()
    else:
        f = open(keggpath+'.xml', 'r')
        kgml = f.read()
        f.close()
    kgtree = et.fromstring(kgml)
    
    if not os.path.isfile(keggpath + ".png"):
        os.system("wget " + kgtree.get("image") + ' -q -O ' + keggpath + ".png")
    im = Image.open(keggpath + ".png")
    im = im.convert('RGBA')
    imdata = np.array(im)   # "imdata" is a height x width x 4 numpy array
    red, green, blue = imdata.T[:3]
    areas = (red < 222) & (blue < 222) & (green > 250)
    imdata[:,:,:3][np.transpose(areas)] = WHITE # Set green region to white.
    areas = (abs(red - green) > 20) | (abs(blue - green) > 20)
    imdata[:,:,:3][np.transpose(areas)] = [240, 240, 240] # Set green region to white.
    pathwaygenes = []
    for entry in kgtree.findall('entry'):
        if entry.get('type')=='gene':
            for gene in re.findall('\d+', entry.get('name')):
                pathwaygenes.append(int(gene))
    maxvalue = data.loc[pathwaygenes][colorcolumn].abs().max()
    xlen = len(imdata[0])
    genesfound = set(pathwaygenes).intersection(data.index)
    changed = []
    for entry in kgtree.findall('entry'):
        if entry.get('type')=='gene':
            samepos = False
            value = 0; nv = 0
            for gene in re.findall('\d+', entry.get('name')):
                gene = int(gene)
#                 print gene
                if gene in genesfound:
                    symbol = data.loc[gene]['symbol']
                    newval = data[colorcolumn][gene]
                    try:
                        if abs(newval) > abs(value): value = newval
                    except ValueError:  # Possible when a gene id matches multiples rows in data.
                        mx = max(newval.abs())
                        newval = list(newval[newval.abs()==mx])[0]
                        if abs(newval) > abs(value): value = newval
                    nv += 1
                    if samepos==False or (abs(newval) < cutoff and abs(value) == abs(newval)):
                        try: 
                            x = float(entry[0].get('x'))
                            y = float(entry[0].get('y'))
                            width = float(entry[0].get('width'))
                            height = float(entry[0].get('height'))
                        except:
                            if entry[0].get('type') == 'line': 
                                x1, y1, x2, y2 = [float(c) for c in entry[0].get('coords').split(',')]
                                x = (x1 + x2) / 2; y = (y1 + y2) / 2; width = 1; height = 1
                            else: print(entry[0].get('type'), 'is not implemented yet.'); raise 
                        elem = {'name':[symbol], 'x':x, 'y':y, 'w':width, 'h':height}
                        if samepos: changed[-1] = elem
                        else: changed.append(elem)
                    else:
                        if abs(newval) < cutoff and abs(value) > abs(newval): continue
                        if nv % 2 == 0: changed[-1]['name'][-1] += " " + symbol
                        else: changed[-1]['name'].append(symbol)
                    samepos = True
            if samepos:
                imdata[int(y-height/2+1):int(y+height/2), 
                       int(x-width /2+1):int(x+width /2), :3] = value2color(value, maxvalue, power)
                   
    for i in range(1,101):
        imdata[25:50, xlen-300-(i-1),:3] = value2color(-i, 100, power)
        imdata[25:50, xlen-300+(i-1),:3] = value2color( i, 100, power)
    im2 = Image.fromarray(imdata)
    draw = ImageDraw.Draw(im2)
    font = ImageFont.truetype(os.path.join(os.path.dirname(__file__), 'HelveticaCY.dfont'), 16)
    draw.text((xlen-350, 5), colorcolumn, (0,0,0), font=font)
    draw.text((xlen-433, 30), str(int(-maxvalue*100-0.5)/100.), value2color(-maxvalue, maxvalue, power), font=font)
    draw.text((xlen-195, 30), str(int( maxvalue*100+0.5)/100.), value2color( maxvalue, maxvalue, power), font=font)
    font = ImageFont.truetype(os.path.join(os.path.dirname(__file__), 'HelveticaCY.dfont'), 13)
    for gene in changed:
        x = gene['x']; w = gene['w']; y = gene['y']; h = gene['h']; nr = len(gene['name'])
        for i in range(nr):
            draw.text((x-w/2+5+(7-len(gene['name'][i]))*3, y-h*(nr/2-i)*0.7), gene['name'][i], (0,60,0), font=font)
    if filename is None: im2.show()
    else: im2.save(filename + ".png")
