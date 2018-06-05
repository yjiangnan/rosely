import pandas as pd
from scipy.special import betainc
from collections import OrderedDict
from .neutralstats import locfdr
import requests, os

try:
    from py2cytoscape.data.cyrest_client import CyRestClient
    from py2cytoscape.data.style import StyleUtil
    from Bio.KEGG import REST
except: pass

__all__ = ['enrichment_p_value', 'download_all_kegg_pathways', 'enrichment_analysis', 'draw_kegg_pathways',
           'delete_all_pathways', 'save_pathway_images', 'load_all_kegg_pathways']

def download_all_kegg_pathways(species_code='mmu'):
    """
    
    """
    pathways_str = REST.kegg_list("pathway", species_code).read()
    pathways = {p.split('\t')[0]:{'name':p.split('\t')[1]} for p in pathways_str.rstrip().split('\n')}
    def get_genes_for(pathways):
        for pathway in pathways:
            pathways[pathway]['gene_id'] = set(); pathways[pathway]['gene_symbol'] = set()
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
                        gene_id, gene_symbol = gene_identifiers.split()
                        pathways[pathway]['gene_id'].add(gene_id)
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

def enrichment_p_value(n, Ltop, N, Lpop):
    """
    n out of Ltop genes in top list, while N out of Lpop genes in background (population).
    Assume that Lpop is large and N/Lpop is an accurate estimate of probability for occurence.
    
    ss.binom.cdf(X, N, p) == betainc(N-X, X+1, 1-p)
    """
    p = 1 - betainc(Ltop-(n+0.25-1), n+0.25, 1 - N/Lpop) # 0.25 is added to avoid p > 0.5
    p = min(1, 2 * p)
    return p

def enrichment_analysis(top_genes, population, pathways, values=None, key='gene_symbol', nbins=30):
    Ltop = len(top_genes); Lpop = len(population)
    ps = {}; ns = {}; Ns = {}; enpath = {}; genes = {}; fold = {}
    for pw in pathways:
        symbol = set(pathways[pw][key])
        gs = symbol.intersection(top_genes)
        n = len(gs)
        N = len(symbol.intersection(population))
        if N > 1 and n / (1+Ltop) >= N / Lpop:
            ps[pw] = enrichment_p_value(n, Ltop, N, Lpop)
            ns[pw] = str(n)+' / '+str(Ltop)
            Ns[pw] = str(N)+' / '+str(Lpop)
            enpath[pw] = pathways[pw]['name'].split(' - ')[0]
            if values is not None:
                genes[pw] = ', '.join([g + ' ' + str(round(values[g], 2)) for g in list(gs)])
            else: genes[pw] = ', '.join(gs)
            fold[pw] = round((n/Ltop) / (N/Lpop), 2)
    fdr = locfdr(ps, nbins=nbins)
    enriched = pd.DataFrame(OrderedDict({'Pathway':enpath, 'p value':ps, 'LFDR':fdr, 
                                   'study':ns, 'background':Ns, 'fold':fold,
                                   'genes & log2 fold changes':genes})).sort_values(by='p value')

    return enriched

def draw_kegg_pathways(enriched, DEG_results, LFDR_cutoff=0.4, by='Controlled z-score', scale_by=10, update_color_only=False):
    cy = CyRestClient()
    if not update_color_only:
        all_suid = cy.network.get_all()
        delete_all_pathways(all_suid)
        
        for pathid in (enriched[enriched['LFDR'] < LFDR_cutoff * 1.1]).index:
            requests.get('http://localhost:1234/keggscape/v1/' + pathid.split(':')[1])
    all_suid = cy.network.get_all()
    enriched['SUID'] = -1
    pd.options.mode.chained_assignment = None
    column = 'color value'
    DEG_results.results[column]  = DEG_results.results[by] ** 3
    DEG_results.results[column] /= DEG_results.results[column].std() * scale_by
    for i in all_suid:
        net = cy.network.create(i)
        table = net.get_node_table()
        pathid = net.get_network_table()['KEGG_PATHWAY_ID'][0]
        enriched['SUID'][pathid] = i
        if by in table or column in table:
            meantable4cy = table
            table[column] = DEG_results.results[column]
        else: 
            meantable4cy = table.merge(DEG_results.results, left_on='KEGG_NODE_LABEL_LIST_FIRST', 
                                   right_index=True).groupby('SUID').mean()
        net.update_node_table(df=meantable4cy, network_key_col='SUID')
        
    cy = CyRestClient()
    # Vizual mapping
    my_kegg_style = cy.style.create('KEGG Style')
    new_defaults = {
        'NODE_FILL_COLOR': '#ffffff',
        'NODE_BORDER_WIDTH' : '1.0',
    }
    my_kegg_style.update_defaults(new_defaults)
    color_gradient = StyleUtil.create_3_color_gradient(min =-1, 
                                                       mid = 0,
                                                       max = 1, 
                                                       colors=('blue', 'white', 'red'))
    my_kegg_style.create_continuous_mapping(column=column, vp='NODE_FILL_COLOR', 
                                            col_type='Double', points=color_gradient)       

def delete_all_pathways(all_suid):
    cy = CyRestClient()
    for i in all_suid:
        net = cy.network.create(i)
        cy.network.delete(net)

def save_pathway_images(enriched, folder='pathways', LFDR_cutoff=0.4):
    cy = CyRestClient()
    for i, suid in enumerate(enriched[enriched['LFDR']<LFDR_cutoff]['SUID']):
        net = cy.network.create(suid)
        name = net.get_network_value('name').decode("utf-8")
        file = open(os.path.join(folder, "%02d. " % (i+1,) + name + ".pdf"), 'wb')
        file.write(net.get_pdf())
        file.close()
