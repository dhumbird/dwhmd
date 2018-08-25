#!/usr/bin/python3

import sys
import subprocess
from lammps import lammps

def writeRecord(f, i, name, ix, x, y, z, vx, vy, vz, m):
    f.write("ATOM  ")                                   #record name (string 6)
    f.write("{0:5d}".format(i))                         #serial number (integer 5)
    f.write(" ")                                        #skip column 12
    f.write('{0:>4} '.format(name))                     #atom name (string 4)
    f.write("UNK")                                      #residue name (string 3)
    f.write(" ")                                        #skip column 21
    f.write(" ")                                        #chain identifier (char)
    f.write("{0:4d}".format(ix))                        #residue sequence number (integer 4)
    f.write(" ")                                        #code for insertion of residues (char)
    f.write("   ")                                      #skip columns 28-30
    f.write ('{0:8.3f}{1:8.3f}{2:8.3f}'.format(x,y,z))  #x,y,z position (angstroms) (8.3, 8.3, 8.3)
    f.write('{0:6.2f}'.format(1.0))                     #occupancy factor (6.2) -- previously used for color defs
    ek=(vx**2+vy**2+vz**2)*m/2.0/9648.531
    f.write('{0:6.2f}'.format(ek))                      #temperature factor (6.2)
    f.write('\n')                                       #end of record


name={}
name[28]="Si"
name[19]="F"
launch = True

args = sys.argv
datfile=(args[1:])
Nmax=0
if len(datfile)==1:
    pdbfile = datfile[0].replace('dat','pdb')
else:
    pdbfile = datfile[0][:datfile[0].find("_")]+'.pdb'
for f in datfile:
    with open(f) as fp:
        for i, line in enumerate(fp):
            if i == 2:
                Nmax=max(Nmax,int(line.split()[0]))
            elif i > 2:
                break

with open(pdbfile,'w') as f:
    f.write("HEADER  PROTEIN\n")
    f.write("COMPND  "+pdbfile+"\n")
    f.write("AUTHOR  ALAN SMITHEE\n")

    for dat in datfile:
        if len(datfile)>1:
            f.write("MODEL ")
            f.write('{0:4d}\n'.format(datfile.index(dat)+1))
        lmp = lammps(cmdargs=["-l", "dat2pdb.log"])
        lmp.commands_list(["units metal","atom_style atomic","atom_modify map array","read_data "+dat])

        pos={}
        vel={}
#        mass={}
        types={}

        x = lmp.extract_atom("x",3) 
        v = lmp.extract_atom("v",3) 
        t = lmp.extract_atom("type",0) 
        ids = lmp.extract_atom("id",0) 
        mass = lmp.extract_atom("mass",2) 
        N=lmp.get_natoms()
        for i in range(0,N):
            pos[ids[i]]=(x[i][0], x[i][1], x[i][2])
            vel[ids[i]]=(v[i][0], v[i][1], v[i][2])
            types[ids[i]]=t[i]
            #mass[ids[i]]=m[i]

        for i in range(1,Nmax+1):
            if i in pos:
                writeRecord(f, i, name[round(mass[types[i]])], ids[i], pos[i][0], pos[i][1], pos[i][2], 
                    vel[i][0], vel[i][1], vel[i][2], mass[types[i]])
            else:
                f.write("ATOM  {0:5d}\n".format(i))
        lmp.close()
        if len(datfile)>1:
            #for i in range(N,Nmax):
            #    f. write("ATOM  {0:5d}\n".format(i))
            f.write("ENDMDL\n")  
    f.write("END\n")    

if launch:
    subprocess.Popen(["vmd",pdbfile])
    sys.exit(0)