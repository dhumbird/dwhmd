#!/usr/bin/python3

import numpy as np
import sys
import math

name={}
mass={}
name[14]="Si"
mass[14]=28.066
name[9]="F"
mass[9]=18.9984
name[18]="Ar"
mass[18]=39.948
name[17]="Cl"
mass[17]=35.453

class Atom:
    def __init__(self):
        self.r=np.zeros(3) 
        self.v=np.zeros(3) 
        self.specie=''
        self.id=0
        self.fix=False

class Config(object):
    def __init__(self, file=None):
        self.N=0
        self.atoms={}
        if file != None:
            with open(file) as f:
                for c in f.read().splitlines():
                    if any(item in c for item in ['#', '%']):
                        if c.find("#N ")>-1:
                            self.N=int(c[3:])
                            for i in range(self.N):
                                self.atoms[i] = Atom() 
                    else:
                        line=(c.replace("<","")).replace(">","").split()
                        self.atoms[int(line[0])].id=int(line[0])+1
                        self.atoms[int(line[0])].specie=int(line[1])
                        self.atoms[int(line[0])].r=np.array([float(line[2]), float(line[3]), float(line[4])])
                        self.atoms[int(line[0])].v=np.array([float(line[5]), float(line[6]), float(line[7])])
                        self.atoms[int(line[0])].fix=(line[8]=='1')

args = sys.argv
cfgfile=str(args[1])
datfile = cfgfile.replace('cfg','dat')

with open(cfgfile) as info:
    for c in info.read().splitlines():
        if c.find("#dimensions ")>-1:
            line = c.split()
            dx=float(line[1])
            dy=float(line[2])
            dz=float(line[3])

cfg = Config(cfgfile)

species = [cfg.atoms[i].specie for i in cfg.atoms]
atype={key: None for key in set(species)}
t=1
atype[14]=t
for k in atype.keys():
    if atype[k] is None:
        atype[k]=t+1
        t+=1

with open(datfile,'w') as f:
    f.write ('hum2lmp: convert '+cfgfile+' to dat\n\n')
    f.write (str(cfg.N)+' atoms\n')
    f.write (str(len(set(species))) +" atom types\n\n")

    f.write ('{0:.5f} {1:.5f} xlo xhi\n'.format(dx*-0.5, dx*0.5))
    f.write ('{0:.5f} {1:.5f} ylo yhi\n'.format(dy*-0.5, dy*0.5))
    f.write ('{0:.5f} {1:.5f} zlo zhi\n'.format(dz*-0.5, dz*0.5+5))

    f.write ('\nMasses\n\n')
    for k in atype.keys():
        f.write (str(atype[k])+' '+str(mass[k])+'\n')
        
    f.write ('\nAtoms\n\n')
    for i in cfg.atoms:
        f.write ('{0:d} {1:d} {2:.5f} {3:.5f} {4:.5f}\n'
                 .format(cfg.atoms[i].id, atype[cfg.atoms[i].specie], 
                  cfg.atoms[i].r[0], cfg.atoms[i].r[1], cfg.atoms[i].r[2]))

    f.write ('\nVelocities\n\n')
    for i in cfg.atoms:
        if cfg.atoms[i].fix:
            f.write ('{0:d} {1:.5f} {2:.5f} {3:.5f}\n'.format(cfg.atoms[i].id, 0, 0, 0))
        elif math.sqrt(pow(cfg.atoms[i].v[0],2)+pow(cfg.atoms[i].v[1],2) + pow(cfg.atoms[i].v[2],2))< 0.001:
            f.write ('{0:d} {1:.5f} {2:.5f} {3:.5f}\n'.format(cfg.atoms[i].id, 0.001, 0.001, 0.001))
        else:
            f.write ('{0:d} {1:.5f} {2:.5f} {3:.5f}\n'
                     .format(cfg.atoms[i].id, cfg.atoms[i].v[0], cfg.atoms[i].v[1], cfg.atoms[i].v[2]))
