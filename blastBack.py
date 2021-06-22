import os
from sequenceAnalyzer import FastAreader
import glob
from pathlib import Path
import sys

class Consensus():
    def __init__(self, NAME):
        self.NAME = NAME
        self.consensusDic = {}       
    def msa(self):
        for i in range(0,100):
            if len(str(i)) == 1:
                x = '0' + str(i)

            else:
                x = str(i)
            file = Path("{0}/{1}.Pass.withcoords".format(self.NAME,x))
            if file.is_file():
                 
                myReaderIEs= FastAreader('{0}/{1}.Pass.withcoords'.format(self.NAME, x))
                with open('{0}/{1}_IEs.fa'.format(self.NAME,x),'w') as f:
                    for header, sequence in myReaderIEs.readFasta():
                        f.write('>{0}\n{1}\n'.format(header,sequence[20:-20]))
  
                os.system('mafft {0}/{1}_IEs.fa > {0}/{1}_IEs_msa.fa'.format(self.NAME, x))
                self.consensusDic[x] = self.consensus(x)
        print(self.consensusDic)
        return self.consensusDic


    def consensus(self,x):
        myReaderIEs= FastAreader('{0}/{1}_IEs_msa.fa'.format(self.NAME, x))
        indexDic = {}
        total = 0
        consensus = ''
        for header, sequence in myReaderIEs.readFasta():
            print(len(sequence))
            total += 1
            for i in range(0,len(sequence)):
                if i not in indexDic:
                    indexDic[i] = {'a':0,'t':0,'g':0,'c':0,'-':0,'n':0}
                if sequence[i] not in indexDic[i]:
                    indexDic[i][sequence[i]] = 0   
                indexDic[i][sequence[i]] += 1
        for index in indexDic:
            for nuc in indexDic[index]:
                if float(indexDic[index][nuc])/float(sum(indexDic[index].values())) > .5 and nuc != '-':
                    consensus += nuc.upper()
        print(len(consensus))
        return consensus




def main():
    with open(sys.args[1],'r') as f:
        for line in f:
            try: 
                sp = line.strip().split('\t')
                dir = sp[0]
                fasta = sp[1]
                
                print(fasta)
                os.system('cp rename.bash {0}'.format(dir))
                os.chdir(dir)
                os.system('pwd')
                os.system('bash rename.sh')
       	        os.chdir('..')

                myAlignments = Consensus(dir)
                consensusDic = myAlignments.msa()
                with open('{0}/IEconsensi.fa'.format(dir),'w') as f:
                    for header in consensusDic:
       	       	        f.write('>{0}\n{1}\n'.format(header, consensusDic[header]))
       	        os.system('makeblastdb -dbtype nucl -in {1} -out {0}/genomeDB'.format(dir,fasta))                    
                os.system('blastn  -db {0}/genomeDB -query {0}/IEconsensi.fa -outfmt 6 -perc_identity 80 -out {0}/blast_back.tsv'.format(dir))
            except IndexError:
                print(dir)                  
  
            
            
if __name__ == "__main__":
    main() 




