#!/usr/bin/python3
carbon_sp3={"c_sp3":["CT1","CT2","CT3","CP1","CP2","CP3","CT","CT1X","CT2X","CT3X"],\
            "c_sp2":["C","CC","CA","CPH1","CPH2","CPT","CY","CD","CE1","CE2","CST"],\
            "c_sp":["CN"],\
            "n_sp3":["N","NH1","NH2","NH3","NY","NP","NR1"],\
            "n_sp2":["NR2","NR3"],\
            "n_sp":["NC","NC2"],\
            "o_sp3":["OH1","OS","OT"],\
            "o_sp2":["O","OB","OC"],\
            "o_sp":[],\
            "s":["S","SM"],\
            "h":["H","HC","HA","HT","HP","HB","HR1","HR2","HR3","HS","HE1","HE2","HA1","HA2","HA3","HF1","HF2"]}
from sys import argv

class alternateList(object):
    def __init__(self,typelist=None):
        self.typelist=typelist

    def ReturnAlternateAtomList(self,name=""):
        for i in self.typelist:
            if name in self.typelist[i]:
                return self.typelist[i]


class bond(object):
    def __init__(self,first="",last="",fc=None,eb=None):
        self.first=first
        self.last=last
        self.fc=fc
        self.eb=eb
    
    def CompareAndAssigne(self,target=None,ForceCheck=False):
        returnval=False
        if self.first==target.first and self.last==target.last:
            self.fc=target.fc
            self.eb=target.eb
            returnval=True
        elif self.first==target.last and self.last==target.first:
            self.fc=target.fc
            self.eb=target.eb
            returnval=True
        return returnval#self.first, self.last, self.fc, self.eb

    def Summary(self):
        print ("First: %s Last: %s fc: %12.5f eb: %12.5f"%(self.first, self.last, self.fc, self.eb))

    def ConvertToString(self):
        return "%s %s %12.5f %12.5f"%(self.first, self.last, self.fc, self.eb)
    
class angle(object):
    def __init__(self,first="",middle="",last="",fc=None,eb=None,u1="",u2=""):
        self.first=first
        self.middle=middle
        self.last=last
        self.fc=fc
        self.eb=eb
        self.u1=u1
        self.u2=u2
    
    def CompareAndAssigne(self,target=None,ForceCheck=False):
        returnval=False
        if self.first==target.first and self.last==target.last and self.middle==target.middle:
            self.fc=target.fc
            self.eb=target.eb
            self.u1=target.u1
            self.u2=target.u2
            returnval=True
        elif self.first==target.last and self.last==target.first and self.middle==target.middle:
            self.fc=target.fc
            self.eb=target.eb
            self.u1=target.u1
            self.u2=target.u2
            returnval=True
        return returnval#self.first, self.middle, self.last, self.fc, self.eb


    def Summary(self):
        print ("First: %s Middle: %s Last: %s fc: %12.5f eb: %12.5f u1: %s u2: %s"%(self.first, self.middle, self.last, self.fc, self.eb, self.u1, self.u2))

    def ConvertToString(self):
        return "%s %s %s %12.5f %12.5f %s %s"%(self.first, self.middle, self.last, self.fc, self.eb, self.u1, self.u2)


class dihedral(object):
    def __init__(self,first="",second="",third="",fourth="",fc=None,eb=None,ag=None):
        self.first=first
        self.second=second
        self.third=third
        self.fourth=fourth
        self.fc=fc
        self.eb=eb
        self.ag=ag
    
    def CompareAndAssigne(self,target=None,ForceCheck=False):
        returnval=False
        if self.first==target.first and self.second==target.second and self.third==target.third and self.fourth==target.fourth:
            self.fc=target.fc
            self.eb=target.eb
            self.ag=target.ag
            returnval=True
        elif self.first==target.fourth and self.second==target.third and self.third==target.second and self.fourth==target.first:
            self.fc=target.fc
            self.eb=target.eb
            self.ag=target.ag
            returnval=True
        elif target.first=='X' and self.second==target.second and self.third==target.third and target.fourth=='X':
            self.fc=target.fc
            self.eb=target.eb
            self.ag=target.ag
            returnval=True
        elif target.first=='X' and self.second==target.third and self.third==target.second and target.fourth=='X':
            self.fc=target.fc
            self.eb=target.eb
            self.ag=target.ag
            returnval=True
        elif target.second=='X' and self.first==target.first and self.fourth==target.fourth and target.third=='X':
            self.fc=target.fc
            self.eb=target.eb
            self.ag=target.ag
            returnval=True
        elif target.second=='X' and self.first==target.fourth and self.fourth==target.first and target.third=='X':
            self.fc=target.fc
            self.eb=target.eb
            self.ag=target.ag
            returnval=True
        return returnval#self.first, self.second, self.third, self.fourth, self.fc, self.eb, self.ag

    def Summary(self):
        print ("First: %s Second: %s Third: %s Fourth: %s fc: %12.5f eb: %12d ag: %12.5f"%(self.first, self.second, self.third, self.fourth, self.fc, self.eb, self.ag))
    
    def ConvertToString(self):
        return "%s %s %s %s %12.5f %12d %12.5f"%(self.first, self.second, self.third, self.fourth, self.fc, self.eb, self.ag)


##x='x'
##y='y'
##a=bond(x,y,0,0)
##b=bond(y,x,10,10)
##print a.CompareAndAssigne(b)
##print a

class parameters(object):
    def __init__(self,unknown=None,known=None,change=None,missing=None):
        self.unknown=unknown
##        self.targetparams=[]
##        self.referenceparams=[]
        self.known=known
        self.change=change
        self.missing=missing
        
        self.targetparams=self.GetAllParams(self.unknown)
        self.referenceparams=self.GetAllParams(self.known)
        self.assigneparam,self.nonassigneparam=self.AllAssigne()
        self.Draft(self.change,self.missing)
##        self.failparams=self.Draft()
##        self.PrintUnassignedDihedrals()
        
        
 
    def GetAllParams(self,targets):
        bonds=[]
        angles=[]
        dihedrals=[]

        for target in targets:
            of=open(target)
            lines=of.readlines()
            i=0
            
            while i<len(lines):
                if lines[i].strip().startswith('BONDS'):
                    i=i+1
                    while 1==1:                    
                        if lines[i].strip().startswith('!'):
                            pass                                           
                        elif lines[i].strip()=='':
                            break
                        else:
                            params=lines[i].strip().split()
    ##                        print params[0],params[1],params[2],params[3]
                            a=bond(first=params[0].strip(),last=params[1].strip(),fc=float(params[2].strip()),eb=float(params[3].strip()))
                            bonds.append(a)
                        i=i+1

                if lines[i].strip().startswith('ANGLES'):
                    i=i+1
                    while 1==1:                    
                        if lines[i].strip().startswith('!'):
                            pass                                           
                        elif lines[i].strip()=='':
                            break
                        else:
                            lines[i]=lines[i].split('!')[0].strip()
                            params=lines[i].strip().split()
                            first=params[0].strip()
                            middle=params[1].strip()
                            last=params[2].strip()
                            fc=float(params[3].strip())
                            eb=float(params[4].strip())
                            if len(params)==7:
                                u1=str(params[5].strip())
                                u2=str(params[6].strip())
                            else:
                                u1=""
                                u2=""
    ##                        print params[0],params[1],params[2],params[3],params[4]
                            a=angle(first=first,middle=middle,last=last,fc=fc,eb=eb,u1=u1,u2=u2)
                            angles.append(a)
                        i=i+1

                if lines[i].strip().startswith('DIHEDRALS'):
                    i=i+1
                    while 1==1:                    
                        if lines[i].strip().startswith('!'):
                            pass                                           
                        elif lines[i].strip()=='':
                            break
                        else:
                            params=lines[i].strip().split()
    ##                        print params[0],params[1],params[2],params[3],params[4],params[5]
                            a=dihedral(first=params[0].strip(),second=params[1].strip(),third=params[2].strip(),fourth=params[3].strip(),fc=float(params[4].strip()),eb=int(params[5].strip()),ag=float(params[6].strip()))
                            dihedrals.append(a)
                        i=i+1
                           
                i=i+1
                         
        return bonds,angles,dihedrals

    def AllAssigne(self):
        assigneparam=[]
        nonassigneparam=[]
        for i in range(len(self.targetparams)):
            tmpassigneparam=[]
            tmpnonassigneparam=[]
            for j in self.targetparams[i]:
                for k in self.referenceparams[i]:
                    x=j.CompareAndAssigne(k)
                    if x==True:
                        tmpassigneparam.append(j)                    
                        break
                if x==False:
                    if i==0:
                        print ("BONDS")
                        first=alternateList(carbon_sp3).ReturnAlternateAtomList(j.first)
                        last=alternateList(carbon_sp3).ReturnAlternateAtomList(j.last)
                        for f in first:
                            for l in last:
                                a=bond(f,l,0,0)
                                for k in self.referenceparams[i]:
                                    x=a.CompareAndAssigne(k)
                                    if x==True:
                                        print (a.first,a.last,a.fc,a.eb)
                                        a=bond(j.first,j.last,a.fc,a.eb)
                                        print (a.first,a.last,a.fc,a.eb)
                                        print ("-------------------------------------------")
                                        tmpassigneparam.append(a)
                                        break
                                if x==True:
                                    break
                            if x==True:
                                break
                        if x==False:
                            tmpnonassigneparam.append(j)

                    elif i==1:
                        print ("ANGLES")
                        first=alternateList(carbon_sp3).ReturnAlternateAtomList(j.first)
                        last=alternateList(carbon_sp3).ReturnAlternateAtomList(j.last)
                        middle=alternateList(carbon_sp3).ReturnAlternateAtomList(j.middle)
                        for f in first:
                            for l in last:
                                for m in middle:
                                    a=angle(f,m,l,0,0,"","")
                                    for k in self.referenceparams[i]:
                                        x=a.CompareAndAssigne(k)
                                        if x==True:
                                            print (a.first,a.middle,a.last,a.fc,a.eb,a.u1,a.u2)
                                            a=angle(j.first,j.middle,j.last,a.fc,a.eb,a.u1,a.u2)
                                            print (a.first,a.middle,a.last,a.fc,a.eb,a.u1,a.u2)
                                            print ("-------------------------------------------")
                                            tmpassigneparam.append(a)
                                            break
                                    if x==True:
                                        break
                                if x==True:
                                    break
                            if x==True:
                                break
                        if x==False:
                            tmpnonassigneparam.append(j)

                    elif i==2:
                        print ("DIHEDRALS")
                        first=alternateList(carbon_sp3).ReturnAlternateAtomList(j.first)
                        second=alternateList(carbon_sp3).ReturnAlternateAtomList(j.second)
                        third=alternateList(carbon_sp3).ReturnAlternateAtomList(j.third)
                        fourth=alternateList(carbon_sp3).ReturnAlternateAtomList(j.fourth)
                        for f in first:
                            for s in second:
                                for t in third:
                                    for fo in fourth:
                                        a=dihedral(f,s,t,fo,0,0,0)
                                        for k in self.referenceparams[i]:
                                            x=a.CompareAndAssigne(k)
                                            if x==True:
                                                print (a.first,a.second,a.third,a.fourth,a.fc,a.eb,a.ag)
                                                a=dihedral(j.first,j.second,j.third,j.fourth,a.fc,a.eb,a.ag)
                                                print (a.first,a.second,a.third,a.fourth,a.fc,a.eb,a.ag)
                                                print ("-------------------------------------------")
                                                tmpassigneparam.append(a)
                                                break
                                        if x==True:
                                            break
                                    if x==True:
                                        break
                                if x==True:
                                    break
                            if x==True:
                                break
                        if x==False:
                            tmpnonassigneparam.append(j)
                    
                    else:
                        tmpnonassigneparam.append(j)


                        
                            
            assigneparam.append(tmpassigneparam)
            nonassigneparam.append(tmpnonassigneparam)
        return assigneparam,nonassigneparam


    def CreateParameterBlock(self,header,parameter):
        lines=[]
        lines.append(header+'\n')
        for p in parameter:
            lines.append(p.ConvertToString()+'\n')
        return lines

    def Draft(self,output,fail):
        outfile=open(output,'w')
        failfile=open(fail,'w')
        ptypes=['BONDS','ANGLES','DIHEDRALS']
        for i in range(len(self.assigneparam)):
            x=self.CreateParameterBlock(ptypes[i],self.assigneparam[i])            
            outfile.writelines(x)
            outfile.write('\n')
        for i in range(len(self.nonassigneparam)):
            x=self.CreateParameterBlock(ptypes[i],self.nonassigneparam[i])            
            failfile.writelines(x)
            failfile.write('\n')
        outfile.close()
        failfile.close()
        return

    def PrintUnassignedDihedrals(self):
        backbones=[]
        for i in self.nonassigneparam[2]:
            backbones.append((i.second,i.third))
        unique=set(backbones)
        for u in unique:
            print ("u")
        return
            
        




  
                    
                    
if __name__=='__main__':
    ana=parameters(unknown=[argv[1]],\
    known=[argv[2]],\
    change=argv[3],\
    missing=argv[4])
    

    

##    for bond in ana.targetparams[0]:
##        bond.Summary()
##    for angle in ana.targetparams[1]:    
##        angle.Summary()
##    for dihedral in ana.targetparams[2]:    
##        dihedral.Summary()

##    for bond in ana.referenceparams[0]:
##        bond.Summary()
##    for angle in ana.referenceparams[1]:    
##        angle.Summary()
##    for dihedral in ana.referenceparams[2]:    
##        dihedral.Summary()

        
##    for bond1 in ana.targetparams[0]:
##        m=0
##        for bond2 in ana.referenceparams[0]:
##            if float(bond1.CompareAndAssigne(bond2)[2]) !=0:
##                print bond1.CompareAndAssigne(bond2)
##                m=m+1
##                break
##        if m==0:
##            print '*********'
##            bond1.Summary()
##            print '*********'
##            
        
