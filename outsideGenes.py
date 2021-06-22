from sequenceAnalyzer import FastAreader
from pathlib import Path
import os
import glob
import sys


def main():
    with open(sys.args[1],'r') as f:
        for line in f:
            sp = line.strip().split('\t')
            dir = sp[0]
            
            fasta = '{0}/{1}'.format(dir,sp[1])
            annotation = '{0}/{1}'.format(dir,sp[2])

            blastBackDic = {}
            ieDic = {}
            geneDic = {}
            with open('{0}/blast_back.tsv'.format(dir),'r') as f:
                for line in f:
                    sp = line.split('\t')
                    fam = sp[0]
                    scaff = sp[1]
                    start = sp[8]
                    stop = sp[9]
                    if scaff not in blastBackDic:
                        blastBackDic[scaff] = [{'fam':fam,'start':int(start),'stop':int(stop)}]

                    else:
       	       	        blastBackDic[scaff].append({'fam':fam,'start':int(start),'stop':int(stop)})

            for i in range(1,100):
                if i < 10:
                    fam = '0{0}'.format(str(i))
                else:
                    fam = str(i)
                ieFile = Path("{0}/{1}.Pass.withcoords".format(dir,fam))
                if ieFile.is_file():
                    myReaderIEs = FastAreader('{0}/{1}.Pass.withcoords'.format(dir, fam))
                    for header, sequence in myReaderIEs.readFasta():
                        coords = header.split(' ')[-1]
                        scaff = coords.split(':')[0]
                        if '+' in coords:
                            start = coords.split('+')[1].split('-')[0]
                            stop = coords.split('+')[1].split('-')[1]
                        else:
                            
                            start = coords.split('-')[1]
                            stop = coords.split('-')[2]              


                        if scaff not in ieDic:
                            ieDic[scaff]=[{'start':int(start),'stop':int(stop)}]
                        else:
                            ieDic[scaff].append({'start':int(start),'stop':int(stop)})
            with open('{0}'.format(annotation),'r') as f:
                for line in f:
                    if '#' in line:
                        pass
                    else:
                        sp = line.split('\t')
                        if sp[2] == 'mRNA':
                            if sp[0] not in geneDic:
                                if int(sp[4])>int(sp[3]):
                                    geneDic[sp[0]] = [sorted([int(sp[3]),int(sp[4])])]
                                else:
                                    geneDic[sp[0]] = [sorted([int(sp[4]),int(sp[3])])]

                            else:
       	       	       	       	if int(sp[4])>int(sp[3]):

       	       	       	            geneDic[sp[0]].append(sorted([int(sp[3]),int(sp[4])]))
                                else:
                                    geneDic[sp[0]].append(sorted([int(sp[4]),int(sp[3])]))


                            
            insideDic = {}
            outsideDic = {}
            #print(geneDic)
            for scaff in blastBackDic:
                outsideDic[scaff] = []
                insideDic[scaff] = []
                #found = False
                for coords in blastBackDic[scaff]:
                    try:
                        if scaff in ieDic:

                            found = False
                            for ie in ieDic[scaff]:
                           # print(ie,coords)
                                if coords['start'] >= ie['start'] and coords['start'] <= ie['stop']:
                                    found = True

                                    break
       	       	                elif coords['stop'] >= ie['start'] and coords['stop'] <= ie['stop']:
                                    found = True 
                                    break
                            if found == False:
                                g = False
                                for gene in geneDic[scaff]:
                                    if coords['stop'] > coords['start']:
                                        if coords['start'] >= gene[0] and coords['stop'] <= gene[1]:
                                        #print(coords['start'], gene, coords['stop'])
                                            insideDic[scaff].append(coords)
                                            g = True
                                            break
                                    else:
                                        if coords['stop'] >= gene[0] and coords['start'] <= gene[1]:
                                            insideDic[scaff].append(coords)
                                            g = True
                                            break
                                if g == False:
                                    outsideDic[scaff].append(coords)
                        else:
                            g = False
                            for gene in geneDic[scaff]:
                                if coords['stop'] > coords['start']:
                                    if coords['start'] >= gene[0] and coords['stop'] <= gene[1]:
                                        insideDic[scaff].append(coords)
                                        g = True
                                        break
                                else:
                                    if coords['stop'] >= gene[0] and coords['start'] <= gene[1]:
                                        insideDic[scaff].append(coords)
                                        g = True
                                        break
                            if g == False:
                                outsideDic[scaff].append(coords)
                    except KeyError:
                        outsideDic[scaff].append(coords)
            

            fastaOutside = {}
            fastaInside = {}           
            myReaderGenome = FastAreader('{0}'.format(fasta))

            for header, sequence in myReaderGenome.readFasta():
                scaff = header.split(' ')[0]
                if scaff in outsideDic:
                    #print(scaff)
                    for coords in outsideDic[scaff]:
                        if coords['stop'] > coords['start']:
                            fastaOutside['{0}_{1}_{2}_{3}'.format(coords['fam'],scaff,str(coords['start']),str(coords['stop']))] = sequence[coords['start']:coords['stop']]
                        else:
       	       	       	    print(sequence[coords['stop']:coords['start']])

                            fastaOutside['{0}_{1}_{2}_{3}'.format(coords['fam'],scaff,str(coords['stop']),str(coords['start']))] = sequence[coords['stop']:coords['start']]
                if scaff in insideDic:

       	       	    #print(scaff)
    
                    for coords in insideDic[scaff]:
                        if coords['stop'] > coords['start']:

                            fastaInside['{0}_{1}_{2}_{3}'.format(coords['fam'],scaff,str(coords['start']),str(coords['stop']))] = sequence[coords['start']:coords['stop']]
                        else:
                            fastaInside['{0}_{1}_{3}_{2}'.format(coords['fam'],scaff,str(coords['start']),str(coords['stop']))] = sequence[coords['stop']:coords['start']]

            with open('{0}/outsidegenes.fa'.format(dir),'w') as f:
                for ie in fastaOutside:
                    f.write('>{0}\n{1}\n'.format(ie,fastaOutside[ie]))
            with open('{0}/insidegenes.fa'.format(dir),'w') as f:
       	       	for ie in fastaInside:
       	       	    f.write('>{0}\n{1}\n'.format(ie,fastaInside[ie]))

                                
            
            
if __name__ == "__main__":
    main() 


