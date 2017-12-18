'''
去除hashtag大于3的微博
'''

import pickle
import sys

find = True if sys.argv[1] == 'find' else False

if find:
    ipf = open('ht', 'r')
    removelist = []
    for idx, line in enumerate(ipf):
        ll = line.split()
        if len(ll) > 3:
            removelist.append(idx)
    pickle.dump(removelist, open('removelist.pkl', 'wb'))
else:
    removelist = pickle.load(open('removelist.pkl', 'rb'))
    ipf = open(sys.argv[2], 'r')
    opf = open('new/'+sys.argv[2], 'w')
    for idx, line in enumerate(ipf):
        if idx not in removelist:
            opf.write(line)
    opf.close()

