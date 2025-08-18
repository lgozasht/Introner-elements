import glob
import os


#for fasta in reversed(glob.glob('assemblies/ncbi_dataset/data/GC*/GC*.fna')):
for fasta in glob.glob('assemblies/ncbi_dataset/data/GC*/GC*.fna'):

    print(fasta)
    
    accession = fasta.split('/')[3]


    if  os.path.isfile('assemblies/ncbi_dataset/data/{0}/done.txt'.format(accession)):
        continue

    print('working on {0}'.format(accession))

    gffName = fasta.split('/')[-1].replace('.fna','.gff')
    os.system('cp FIND_CLUSTERS_ALLEUKS.pl  assemblies/ncbi_dataset/data/{0}/FIND_CLUSTERS_ALLEUKS.pl'.format(accession))     
    if os.path.isfile('assemblies/ncbi_dataset/data/{0}/{1}'.format(accession, gffName)):
        continue
    elif os.path.isfile('assemblies/ncbi_dataset/data/{0}/genomic.gff'.format(accession)):
        os.system('mv assemblies/ncbi_dataset/data/{0}/genomic.gff assemblies/ncbi_dataset/data/{0}/{1}'.format(accession, gffName))

    if os.path.isfile('assemblies/ncbi_dataset/data/{0}/{1}'.format(accession, gffName)):
        os.chdir('assemblies/ncbi_dataset/data/{0}/'.format(accession)) 
        os.system('perl FIND_CLUSTERS_ALLEUKS.pl {0}'.format(gffName))
            
        fasta = fasta.split('/')[-1]
        os.system('rm -r {0}.exons-introns'.format(fasta.replace('.fna','')))
        os.system('rm -r {0}.intronsflanks'.format(fasta.replace('.fna','')))
        os.system('rm -r {0}.pro*'.format(fasta.replace('.fna','')))
        os.system('rm -r {0}.ends*'.format(fasta.replace('.fna','')))
        os.system('rm -r {0}.paralogs'.format(fasta.replace('.fna','')))
        os.system('rm -r {0}.intronmatches'.format(fasta.replace('.fna','')))
        os.system('rm -r {0}.NoGroup'.format(fasta.replace('.fna','')))
        os.system('rm -r FIND_CLUSTERS_ALLEUKS.pl')
        if len(glob.glob('{0}*.Pass'.format(fasta.replace('.fna','')))) == 0:
            print('no introners')
            os.system('rm -r {0}'.format(fasta))
            os.system('rm -r {0}'.format(gffName))

        with open('done.txt','w') as f:
            f.write(' ')
       
        os.chdir('../../../../')
