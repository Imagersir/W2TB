import numpy as np
import sympy as sp
from tqdm import tqdm

def Schmidt_orth(eigenvector):
    num = eigenvector.shape[1]
    for i in range(num):
        for j in range(i):
            eigenvector[:,i] = eigenvector[:,i] - eigenvector[:,j]*np.dot(eigenvector[:,j].conj().T, eigenvector[:,i])/(np.dot(eigenvector[:,j].conj().T,eigenvector[:,j]))
        eigenvector[:, i] = eigenvector[:, i]/np.linalg.norm(eigenvector[:, i])
    return eigenvector
def step_abs(x,m):
    return np.array(np.abs(x) >= m, dtype = np.int32)
def blockdiag(lm):
    n=len(lm)
    tm=[]
    ll=[]
    cl=[]
    for m in lm:
        (l,c)=m.shape
        ll.append(l)
        cl.append(c)
    for i in range(n):
        m=lm[i]
        t=[]
        for j in range(n):
            if i==j:
                t.append(m)
            else:
                t.append(np.zeros((ll[i],cl[j]),dtype=np.complex64))
        tm.append(t)
    return np.block(tm)
class tb_h:
    def print_inf(self):
        print("The projected orbitals of each atom:")
        print(self.orinf)
        print("The number of the atoms:")
        if self.spinors:
            for i in range(self.num_atoms):
                if i < round(self.num_atoms/2):
                    print([i,self.anlist[i],'↑'])
                else:
                    print([i,self.anlist[i],'↓'])
        else:
            for i in range(self.num_atoms):
                print([i,self.anlist[i]])
    def __init__(self,fermi=0,seedname="wannier90",dn=2):
        self.dn=dn
        self.seedname=seedname
        self.k1,self.k2,self.k3=sp.symbols("kx ky kz", real=True)
        self.orinf={}
        self.num_wann=0
        self.num_atoms=0
        self.anlist=[]
        self.factor=[]
        self.spinors=False
        self.R=[]
        wannier_hr=seedname+"_hr.dat"
        with open(wannier_hr, "r") as f1:
            lines = f1.readlines()
            f1.close()
        with open(seedname+".win", "r") as f2:
            contents = f2.readlines()
            f2.close()
        self.atom_quantities = {}
        projections_dict = {}
        bv = []
        in_projections = False
        in_atoms_cart = False
        in_unit_cell_cart = False

        for line in contents:
            stripped_line = line.strip()

            # Check for num_wann
            if stripped_line.startswith('num_wann'):
                self.num_wann = int(stripped_line.split('=')[1].strip())
            # Check for spinors
            if stripped_line.startswith('spinors'):
                spin = stripped_line.split('=')[1].strip()
                if spin=='.true.':
                    self.spinors = True
                if spin=='T':
                    self.spinors = True
                if spin=='True':
                    self.spinors = True
                if spin=='.T.':
                    self.spinors = True

            # Projections section
            elif stripped_line == 'begin projections':
                in_projections = True
            elif stripped_line == 'end projections':
                in_projections = False
            elif stripped_line.startswith('begin unit_cell'):
                in_unit_cell_cart = True
            elif stripped_line.startswith('end unit_cell'):
                in_unit_cell_cart = False
            elif in_projections:
                projection_parts = stripped_line.split(':')
                if len(projection_parts) == 2:
                    atom_type, projection_str = projection_parts
                    atom_type = atom_type.strip()
                    projections = [proj.strip() for proj in projection_str.split(';')]
                    if 'p' in projections:
                        projections.remove('p')
                        projections.extend(['pz','px','py'])
                    if 'd' in projections:
                        projections.remove('d')
                        projections.extend(['dz2','dxz','dyz','dxy','dx2y2'])
                    projections_dict[atom_type]=projections
            elif in_unit_cell_cart:
                    if len(stripped_line.split())==3:
                        bv.append([float(x) for x in stripped_line.split()])
            # Atoms_cart section for atom quantities
            elif stripped_line.startswith('begin atoms'):
                in_atoms_cart = True
            elif stripped_line.startswith('end atoms'):
                in_atoms_cart = False
            elif in_atoms_cart:
                atom_line_parts = stripped_line.split()
                if len(atom_line_parts) > 1:
                    self.num_atoms+=1
                    atom_type = atom_line_parts[0]
                    self.anlist.append(atom_type)
                    if atom_type not in self.atom_quantities:
                        self.atom_quantities[atom_type] = 0
                    self.atom_quantities[atom_type] += 1
        if self.spinors:
            self.anlist+=self.anlist
            self.num_atoms*=2
        self.basis_vector=bv
        self.orinf=projections_dict
        self.V=np.dot(bv[0], np.cross(bv[1], bv[2]))
        self.fermi=fermi
        if abs(self.V)<0.01:
            raise Exception("Basis vectors may be wrong because the volume of cell is too small!")
        self.rec=[np.cross(bv[1], bv[2]) * 2 * np.pi / self.V,np.cross(bv[2], bv[0]) * 2 * np.pi / self.V,np.cross(bv[0], bv[1]) * 2 * np.pi / self.V]
        for i in range(self.num_wann):
            self.factor.append([])
            self.R.append([])
            for j in range(self.num_wann):
                self.factor[len(self.factor) - 1].append([])
                self.R[len(self.R) - 1].append([])
        for l in lines:
            ls=l.split()
            if len(ls) == 7:
                self.factor[int(ls[3]) - 1][int(ls[4]) - 1].append(round(float(ls[5]),dn) + 1j * round(float(ls[6]),dn))
                self.R[int(ls[3]) - 1][int(ls[4]) - 1].append([int(ls[0]), int(ls[1]), int(ls[2])])
        self.a2o=[]
        osum=0
        for i in range(self.num_atoms):
            na = self.anlist[i]
            tl=[]
            no=len(self.orinf[na])
            for j in range(no):
                tl.append(j+osum)
            osum+=no
            self.a2o.append(tl)
        self.H = [[[0 for i in range(7)] for j in range(7)] for k in range(7)]
        self.print_inf()
        print('Constructing the Hmatrix')
        for r1 in range(-3,4):
            for r2 in range(-3,4):
                for r3 in range(-3,4):
                    rt=[r1,r2,r3]
                    Htemp=np.zeros((self.num_wann, self.num_wann),dtype='complex')
                    for i in range(self.num_wann):
                        for j in range(self.num_wann):
                            Rl=self.R[i][j]
                            if rt in Rl:
                                k=Rl.index(rt)
                                f=self.factor[i][j][k]
                                if (rt==[0,0,0])&(i==j):
                                    f-=self.fermi
                                Htemp[i,j] = f
                    self.H[r1+3][r2+3][r3+3]=np.round(Htemp,self.dn)
        print('The Hmatrix has been constructed')
    def output_hr(self):
        print('Start to construct the file hr.dat')
        ff = open(self.seedname+"_new_hr.dat", "w+")         # 返回一个文件对象 
        ff.write(" writen by W2TB\n") 
        ff.write("   %i"%self.num_wann)
        nr=7*7*7
        ff.write("   %i"%nr)
        for i in range(nr):
            if i%15==0:
                ff.write('\n')
            ff.write(" %i"%1)
        db=100.0/nr/self.num_wann
        with tqdm(total=100) as pbar:
            for r1 in range(-3,4):
                for r2 in range(-3,4):
                    for r3 in range(-3,4):
                        for i in range(self.num_wann):
                            for j in range(self.num_wann):
                                hr=np.real(self.H[r1+3][r2+3][r3+3][i,j])
                                hi=np.imag(self.H[r1+3][r2+3][r3+3][i,j])
                                ff.write("\n %i %i %i %i %i %f %f "%(r1,r2,r3,i+1,j+1,hr,hi))
                            pbar.update(db)
        
        ff.close()
        print('The file hr.dat has been constructed')
    def select_Hblock(self,na,nb,r1=[0,0],r2=[0,0],r3=[0,0],obl=[],dn=None,nm=0,k=[]):
        if dn==None:
            dn=self.dn
        if obl==[]:
            oa=self.a2o[na]
            ob=self.a2o[nb]
        else:
            oa=self.a2o[na][obl[0][0]:obl[0][-1]+1]
            ob=self.a2o[nb][obl[1][0]:obl[1][-1]+1]
        if k==[]:
            Hsum=sp.zeros(len(oa),len(ob))
        else:
            Hsum=np.zeros((len(oa),len(ob)),dtype=np.complex64)
        for ri in range(r1[0],r1[-1]+1):
            for rj in range(r2[0],r2[-1]+1):
                for rk in range(r3[0],r3[-1]+1):
                    temp=step_abs(self.H[ri+3][rj+3][rk+3][oa[0]:oa[-1]+1,ob[0]:ob[-1]+1],nm)
                    Htemp=np.round(self.H[ri+3][rj+3][rk+3][oa[0]:oa[-1]+1,ob[0]:ob[-1]+1]*temp,dn)
                    if k==[]:
                        eikr=sp.exp(sp.I*(ri*self.k1+rj*self.k2+rk*self.k3))
                        h=Htemp.tolist()
                        Hsym=sp.Matrix(h)*eikr
                    else:
                        eikr=np.round(np.exp(2j*np.pi*(ri*k[0]+rj*k[1]+rk*k[2])),dn)
                        Hsym=Htemp*eikr
                    Hsum+=Hsym
        return Hsum
    def show_Hmatrix(self,r1=[0,0],r2=[0,0],r3=[0,0],rr=[],cr=[],mm=[],m0=0,dn=None,nm=0,k=[]):
        if dn==None:
            dn=self.dn
        if rr==[]:
            rr=range(self.num_atoms)
        if cr==[]:
            cr=rr
        Hlist=[]
        if m0!=0:
            if m0!=2:
                for i in rr:
                    htemp=[]
                    ni=len(self.a2o[i])
                    for j in cr:
                        nj=len(self.a2o[j])
                        htemp.append(self.select_Hblock(i,j,r1=r1,r2=r2,r3=r3,obl=[[0,0],[0,0]],dn=dn,nm=nm,k=k))
                    if m0!=1:
                        for j in cr:
                            nj=len(self.a2o[j])
                            htemp.append(self.select_Hblock(i,j,r1=r1,r2=r2,r3=r3,obl=[[0,0],[1,nj]],dn=dn,nm=nm,k=k))
                    Hlist.append(htemp)
            if m0!=1:
                for i in rr:
                    htemp=[]
                    ni=len(self.a2o[i])
                    if m0!=2:
                        for j in cr:
                            nj=len(self.a2o[j])
                            htemp.append(self.select_Hblock(i,j,r1=r1,r2=r2,r3=r3,obl=[[1,ni],[0,0]],dn=dn,nm=nm,k=k))
                    for j in cr:
                        nj=len(self.a2o[j])
                        htemp.append(self.select_Hblock(i,j,r1=r1,r2=r2,r3=r3,obl=[[1,ni],[1,nj]],dn=dn,nm=nm,k=k))
                    Hlist.append(htemp)
        else:
            for i in rr:
                htemp=[]
                for j in cr:
                    htemp.append(self.select_Hblock(i,j,r1=r1,r2=r2,r3=r3,obl=[],dn=dn,nm=nm,k=k))
                Hlist.append(htemp)
        for m in mm:
            mi=m[0]
            mj=m[1]
            mat=m[2].tolist()
            if k==[]:
                Hlist[mi][mj]=sp.Matrix(mat)
                Hlist[mj][mi]=Hlist[mi][mj].H
            else:
                Hlist[mi][mj]=np.array(mat)
                Hlist[mj][mi]=Hlist[mi][mj].conj().T
        if k==[]:
            Hs=sp.Matrix(Hlist)
        else:
            Hs=np.block(Hlist)
        return Hs
    def get_eigs(self, kv, R1r, R2r, R3r, udiag=[], rr=[] ,mm=[], m0=0):
        if rr==[]:
            rr=range(self.num_atoms)
        Hs=self.show_Hmatrix(r1=R1r,r2=R2r,r3=R3r,rr=rr,cr=rr,m0=m0,mm=mm,dn=self.dn,nm=0,k=kv)
        if (udiag!=[]):
            Htt=Hs+np.diag(udiag)
        else:
            Htt=Hs
        eig,V=np.linalg.eig(Htt)
        idx = np.argsort(eig)
        eig = eig[idx]
        V=V[:,idx]
        return [eig,V]
    def band_cal(self,kpath,klabels,udiag=[],kn=40,R1r=[-2,2],R2r=[-2,2],R3r=[-2,2],rr=[],mm=[],m0=0):
        if rr==[]:
            rr=range(self.num_atoms)
        k_point = []
        kline=[]
        kticks=[0.0]
        ls=0
        for i in range(len(kpath) - 1):
            interval = np.array(kpath[i + 1]) - np.array(kpath[i])
            k1_vector = interval[0] * np.array(self.rec[0])
            k2_vector = interval[1] * np.array(self.rec[1])
            k3_vector = interval[2] * np.array(self.rec[2])
            k_vec = k1_vector + k2_vector + k3_vector
            tl=np.sqrt(np.dot(k_vec,k_vec))
            interval = interval / kn
            for j in range(kn + 1):
                k_point.append(np.array(kpath[i]) + j * interval)
                kline.append(j*1.0*tl/kn+ls)
            ls+=tl
            kticks.append(ls)
        solution = []
        wl=[]
        no=0
        for i in rr:
            no+=len(self.a2o[i])
        if (m0==0)|(m0==3):
            for i in range(no):
                solution.append([])
                wl.append([])
        elif m0==1:
            for i in range(len(rr)):
                solution.append([])
                wl.append([])
        elif m0==2:
            for i in range(no-len(rr)):
                solution.append([])
                wl.append([])
        #print('Constructing the matrix')
        kcount=0.0
        for kv in k_point:
            kcount+=100.0/len(k_point)
            slvs=self.get_eigs(kv.tolist(),R1r,R2r,R3r,udiag,rr,mm,m0)
            eig = slvs[0]
            vecs= slvs[1]    
            for i in range(len(eig)):
                solution[i].append(eig[i])
                v=vecs[:,i]
                wl[i].append(v)
        bandsave=[kline,solution,kticks,klabels,wl]
        return bandsave
