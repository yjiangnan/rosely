# rosely
### Robust and sensitive tools for high-throughput data analysis

## Installation
Simply copy the folder in any directory. When using, add the directory to the path. For example, I have rosely under /scripts folder:
    import sys
    if 'rudolphLab/scripts' not in sys.path: sys.path.insert(0, '/Users/jiangnan/Documents/rudolphLab/scripts')
    from rosely import *

Importantly, ascertained $t$-test in rosely depends on pyloess that has only python2 support. If you use python3, execnet has to be installed on both python2 and python3:
    pip2 install execnet
    pip3 install execnet
Additionally, python2 should have numpy installed.
