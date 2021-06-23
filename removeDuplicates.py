from pathlib import Path
from sequenceAnalyzer import FastAreader
import os
import sys

outDic = {}
inDic = {}

with open('sys.args[1]','r') as f:
    for line in f:
        removeList = []
        outDic = {}
        inDic = {}

        species = line.strip()
        if os.stat('{0}/outsidegenes.fa'.format(species)).st_size > 0:
            myReaderGenome = FastAreader('{0}/outsidegenes.fa'.format(species))
            for header, sequence in myReaderGenome.readFasta():
                outDic[header] = sequence
       	if os.stat('{0}/insidegenes.fa'.format(species)).st_size > 0:
 
            myReaderGenome = FastAreader('{0}/insidegenes.fa'.format(species))
            for header, sequence in myReaderGenome.readFasta():
                inDic[header] = sequence
        for head1 in outDic:
            for head2 in outDic:
                if head1 != head2:
                    start1 = int(head1.split('_')[-2])
                    stop1 = int(head1.split('_')[-1])

                    start2 = int(head2.split('_')[-2])
                    stop2 = int(head2.split('_')[-1])
                    if (start1 <= start2 and start2 <= stop1) or (start1 <= stop2 and stop2 <= stop1):
                        removeList.append(head2)
        for key in removeList:

            try:
                outDic.pop(key)
            except KeyError:
       	       	pass

       	        #print(key)

        removeList = []
        for head1 in inDic:
            for head2 in inDic:
                if head1 != head2:
                    start1 = int(head1.split('_')[-2])
                    stop1 = int(head1.split('_')[-1])

                    start2 = int(head2.split('_')[-2])
                    stop2 = int(head2.split('_')[-1])
                    if (start1 <= start2 and start2 <= stop1) or (start1 <= stop2 and stop2 <= stop1):
       	       	       	removeList.append(head2)
        for key in removeList:
            try:
                inDic.pop(key)
            except KeyError:
                pass
                #print(key)
        print(species,'\t',len(outDic),'\t',len(inDic))
        with open('{0}/outsideDeduped.fa'.format(species),'w') as f:
            for head in outDic:
                f.write('>{0}\n{1}\n'.format(head,outDic[head]))
       	with open('{0}/insideDeduped.fa'.format(species),'w') as f:
       	    for	head in	inDic:
       	       	f.write('>{0}\n{1}\n'.format(head,inDic[head]))

