import math
import re
class Protein_Specification:
    def __init__(self,filename):
        self.filename=filename
    #main function to find everything asked in question and extra by taking dictionary as input
    def mainproteinspecs(self,coordinates):
        #cetregs function to calculate centroid,geometric size, k(use to find maximal distant residue),surface area, volume, compactness
        def centregs(xyzlist):
            sumx=0
            sumy=0
            sumz=0
            for i in range(0,len(xyzlist),3):
                sumx=sumx+xyzlist[i]
                sumy=sumy+xyzlist[i+1]
                sumz=sumz+xyzlist[i+2]
            centre={}
            n=int(len(xyzlist)/3)
            #centroid coordinates are mean of all x , mean of all y and mean of all z coordinates of residues
            centre["x"]=sumx/n
            centre["y"]=sumy/n
            centre["z"]=sumz/n
            print("\tCentre of protein is: ",centre)
            dfc={}
            j=1
            #for loop to calculte distance from every point to centroid
            for i in range(0,len(xyzlist),3):
                    dfc[j]=pow((pow(centre["x"]-xyzlist[i],2)+pow(centre["y"]-xyzlist[i+1],2)+pow(centre["z"]-xyzlist[i+2],2)),0.5)
                    j=j+1
            gs=list(dfc.values())
            k=gs.index(max(gs))
            #maximal residue will be find out by this k index
            geometricsize=max(gs) #geometric size can be find out by using maximum of euclidean distance between centroid and residues 
            volume=(4*(math.pi*pow(geometricsize,3)))/3  #geometric size can be considered as radius of sphere and volume and surface area can be find out
            surfacearea=(4*(math.pi*pow(geometricsize,2)))
            if(volume!=0):
                compactness=surfacearea/volume
            else:   #compactness will become infinity
                compactness=float('inf')
            specification=[geometricsize,volume,surfacearea,compactness,k,centre]
            return(specification)
        print("\nAll residues specification:\n")
        xyzlist=[]
        #xyzlist is list of coordinates only from coordinates dictionary which will be given as input to centregs function
        for i in coordinates.keys():
            try:
                xyzlist=xyzlist+[coordinates[i]["H"]]
            except:
                xyzlist=xyzlist+[coordinates[i]["P"]]
        allre=centregs(xyzlist)
        amdr=[]
        for i in coordinates.keys():
            if str(allre[4]+1) in i:
                amdr=amdr+[i]
                break
        #amdr is maximal distant residue
        print("\tMaximal distant residue is : ",amdr[0][0:3]," sequence number is: ",amdr[0][4:])
        print("\tGeometric size of protein(all residues): ",allre[0])
        print("\tCompactness of protein(all residues): ",allre[3],"\n")
        print("\nHydrophobic residues specification:\n")
        xyzlist=[]
        #list of hydrophobic residues coordinates from coordinates dictionary
        for i in coordinates.keys():
            try:
                xyzlist=xyzlist+[coordinates[i]["H"]]
            except:
                continue
        hydre=centregs(xyzlist)
        v=xyzlist[3*(hydre[4]+1)-1]
        hmdr=[]
        #hmdr is maximal distant residue
        for i in coordinates.keys():
            try:
                if coordinates[i]["H"]==v:
                        hmdr=hmdr+[i]
                        break
            except:
                continue
        print("\tMaximal distant residue is : ",hmdr[0][0:3]," sequence number is: ",hmdr[0][4:])
        print("\tGeometric size of protein(hydrophobic residues): ",hydre[0])
        print("\tCompactness of protein(hydrophobic residues): ",hydre[3],"\n")
        fractionhydvol=hydre[1]/allre[1]   #fraction of hydrophobic volume
        print("\nMain Specifications are: ")
        print("\ti).Centre of protein is: ",allre[5])
        print("\ti).Geometric size of protein(all residues): ",allre[0])
        print("\tii).Fraction of hydrophobic core volume is: ",fractionhydvol)
        print("\tiii).Compactness of protein(all residues): ",allre[3],"\n")
    def dictpdb(self):
        pdbcoordinates={} #main dictionary with sequential amino acids, hydrophobic or polar residues and coordinates as values
        i=1
        aahp = {"ALA":"P","GLY":"P","THR":"P","SER":"P","ASN":"P","GLN":"P",
                "ASP":"P","GLU":"P","HIS":"P","ARG":"P","LYS":"P","PRO":"P",           #dictionary defining amino acid to polar or hydrophobic
                "CYS":"H","MET":"H","PHE":"H","ILE":"H","LEU":"H","VAL":"H",
                "TRP":"H","TYR":"H"}
        pdb_2lwz = open("S:\\Projects\\Protein Specification\\"+self.filename+".pdb", 'r')  #assignment pdb file can also be taken as user input by giving path of file
        for data in pdb_2lwz:
            if(re.search(r'^ATOM\s+\d+\s+CA\s+', data)):  #using regex searched for the carbon alpha atoms lines
                list = re.split('[\s]+', data)            #each line with carbon alpha atom is taken as a list with elements are recognised with space seperation
                for j in aahp.keys():                     #loop to find corresponding amino acid nature wheter hydrophobic or polar from aahp dictionary
                    if j==list[3]:
                        y=aahp[j]
                        break
                pdbcoordinates[list[3]+"x"+str(i)]={}
                pdbcoordinates[list[3]+"x"+str(i)][y]=float(list[6])          #dictionary values with xyz coordinates is taken as input in each loop
                pdbcoordinates[list[3]+"y"+str(i)]={}
                pdbcoordinates[list[3]+"y"+str(i)][y]=float(list[7])
                pdbcoordinates[list[3]+"z"+str(i)]={}
                pdbcoordinates[list[3]+"z"+str(i)][y]=float(list[8])
                i=i+1
                if list[5]=="57":                                            #end of first model by putting this if statement and break in comment we can work on all models
                    break
        print(pdbcoordinates)
        return(pdbcoordinates) 
def main():
    filename="2lwz"
    obj=Protein_Specification(filename)   
    obj.mainproteinspecs(obj.dictpdb())             # main function is called
if __name__=="__main__":
    main()
                                         