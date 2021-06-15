import sys
from ete3 import Tree


"""

Removes parenthisis in node names that perturb the integrity of tree structure 

"""


def replacer(s, newstring, index, nofail=False):
    # raise an error if index is outside of the string
    if not nofail and index not in range(len(s)):
        raise ValueError("index outside given string")

    # if not erroring, but the index is still not in the correct range..
    if index < 0:  # add it to the beginning
        return newstring + s
    if index > len(s):  # add it to the end
        return s + newstring

    # insert the new string between "slices" of the original
    return s[:index] + newstring + s[index + 1:]

def clean(line,index,subtract):

    for i in range(len(line)):
        k = i + index
        try:

            if line[k] == '(' and line[k-1] != ',' and line[k-1] != '(' and line[k+1] != '(':
                with open('log.txt','a') as f:
                    subsection = line[k:k+30]
                    if subsection.find(')') < 0:
                        f.write('Warning: {0}\n'.format(line[k-50:k+50]))
                    f.write('replacing {0}\n'.format(line[k-5:k+30]))

                    line = replacer(line,'',k)
                    subsection = line[k:k+50]

                    line = replacer(line,'',k+subsection.find(')')) 
                    print(line[k-10:k+30])
                   
                    f.write('with {0}\n\n'.format( line[k-5:k+30])) 
                    subtract += 2

        except IndexError:
            print(k)
    return line


with open(sys.args[1],'r') as f:
    for line in f:
#        editLine = line.replace("\'","").replace(" ","_").replace(":","_").replace('_(T..Lally_760)_','_').replace('_(T.R._Lally_760)_','_').replace('_(Nippo-oleifera_Group)_','_').replace('_(subg._Metrosideros)_','_').replace('_(subg._Mearnsia)_','_').replace('_(28.11.2014)_','_').replace('_(b)_','_').replace('_(i)_','_').replace('_(g)_','_').replace('_(d)_','_').replace('_(e)_','_').replace('_(c)_','_').replace('_(h)_','_').replace('_(a)_','_').replace('_(M.E._Trudgen_MET_5369)_','_').replace('_(in__Viridiplantae)_','_').replace('_(type_1)_','_').replace('_(type_2)_','_').replace('_(Short_3969)_','_').replace('_(Watanabe_9)_','_').replace('_(E.Wimm.)_','_').replace('_(HUEFS)_','_').replace('_(BGM)_','_').replace('_(L.A.Craven_2357)_','_').replace('_(S.D._Hopper_1786)_','_').replace('_(Poeae_type)_','_').replace('_(Bory)_','_').replace('_(Laegaard_70382)_','_').replace('_(R.T.Wills_1423)_','_').replace('_(W.G.Trapnell)_','_').replace('_(W.G.Trapnell_269)_','_').replace('_(Miq.)_','_').replace('_(Engl.)_Engl.,','_').replace('_(Engl.)_','_').replace('_(species_A)_','_').replace('_(species_B)_','_').replace('B_(Long_33796)_','_').replace('_(Linis_1459-05)_','_').replace('_(Buck_23824)_','_').replace('s_(D.G.Long)_','_').replace('_(robust_form)_','_').replace('_(clade_I)_','_').replace('_(clade_H)_','_').replace('_(clade_A)_','_').replace('_(clade_B)_','_').replace('_(clade_D)_','_').replace('_(clade_C)_','_').replace('_(clade_E)_','_').replace('_(clade_F)_','_').replace('_(in__Rhodophyta)_','_').replace('_(in__Bacillariophyta)_','_').replace('_(in__Vampyrellida)_','_').replace('_(in__Phaeophyceae)_','_').replace('_(in__Xanthophyceae)_','_').replace('_(SA)_','_').replace('_(lineage_3)_','_').replace('_(lineage_2)_','_').replace('_(lineage_1)_','_').replace('_(lineage_4)_','_').replace('_(in__Ciliophora)_','_').replace('_(Haemamoeba)_','_').replace('_(Novyella)_','_').replace('_(Vinckeia)_','_').replace('_(FVO)_','_').replace('_(2000708)_','_').replace('_(Laverania)_','_').replace('_(Paraplasmodium)_','_').replace('_(Sauramoeba)_','_').replace('_(Huffia)_','_').replace('_(Giovannolaia)_','_').replace('_(Bennettinia)_','_').replace('_(Lacertamoeba)_','_').replace('_(Carinamoeba)_','_')
         
        cleanedLine = clean(line,0,0)

        with open('{0}.clean'.format(sys.args[1]),'w') as f3:

            f3.write(cleanedLine)
        break

t = Tree('{0}'.format(sys.args[1]),format=8)


