# rosely
### Robust and sensitive tools for high-throughput data analysis

## Installation
Simply copy the folder in any directory. When using, add the directory to the path. For example, I have rosely under /scripts folder:

    import sys
    if 'rudolphLab/scripts' not in sys.path: sys.path.insert(0, '/Users/jiangnan/Documents/rudolphLab/scripts')
    from rosely import *

Importantly, ascertained t-test in rosely depends on a specific version of pyloess https://github.com/jcrotinger/pyloess for python3 support.
