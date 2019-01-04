from lammps import lammps

atomtype={}
atomtype[1]="Si"
atomtype[2]="F"
atomtype[3]="Cl"

def init(lmp, file, log):
  timestep = 0.001
  lmp.commands_list(["units metal",
                     "atom_style atomic",
                     "atom_modify map array",
                     "read_data "+file,
                     "pair_style      tersoffHG",
                     "variable vm atom sqrt(vx*vx+vy*vy+vz*vz)",
                     "group unfixed variable vm",
                     "group fixed subtract all unfixed",
                     "fix             fxSurf unfixed nve",
                     "timestep        "+str(timestep),
                     "neighbor        1.0 bin",
                     "neigh_modify    every 1 delay 10 check yes",
                     "compute pea all pe/atom",
                    ])
  log.write("* Initializing config from file "+file+".\n")

######################################
def run(lmp, timesteps):
  if lmp.extract_global("ntypes",0)==1:
      lmp.command("pair_coeff * * SiF.tersoffHG Si")
  else:
      lmp.command("pair_coeff * * SiCl.tersoffHG Si Si Cl")
  lmp.commands_list(["run "+str(timesteps)])
######################################

def thermo(lmp, seed, Temp):
    lmp.commands_list(["velocity unfixed create "+str(Temp)+" "+str(seed)+" rot no dist gaussian"])

######################################

def thermalV(seed, T, m):
  import random
  import math
  random.seed(a=seed)
  vx=0
  vy=0
  vz=0
  kT=math.sqrt(8.6173857E-5 * T / m)
  for q in range (0,12): 
    vx+=random.random()
    vy+=random.random()
    vz+=random.random()
  vx=(vx-6)*kT
  vy=(vy-6)*kT
  vz=(vz-6)*kT
  return vx, vy, vz


def addion(lmp, seed, log, e=0, T=300):
  import math
  ionType = 3
#  ionType = 2
  log.write("* Adding species "+str(ionType)+". (# Nmax) Seed: "+str(seed)+"\n")

  lmp.commands_list(["variable zmin equal bound(all,zmin)",
                     "variable zmax equal bound(unfixed,zmax)"])
  zmax = lmp.extract_variable("zmax",None,0)
  zmin = lmp.extract_variable("zmin",None,0)

  lmp.commands_list([
    "change_box all boundary p p fm z final "+str(zmin)+" "+str(zmax+3),    # Add 3A to the top boundary
    "region vac block EDGE EDGE EDGE EDGE "+str(zmax+2.8)+" "+ str(zmax+3), # Region vac is a 0.2A layer under the boundary
    "create_atoms "+str(ionType)+" random 1 "+str(seed)+" vac",             # Create atom in region vac
    "group ion region vac",                                                 # Add the new atom to group ion
    "group unfixed region vac",                                             # Add the new atom to NVE group unfixed
    "compute ionpe ion reduce sum c_pea",                                   # Track the PE of the atom (zero if no neighbors)
    "compute vIon ion reduce ave vz",
  ])
  mass = lmp.extract_atom("mass",2) 
  vx, vy, vz = thermalV(seed+1, 300, mass[ionType])
  if e==0:
    vz=-abs(vz)
  lmp.command("velocity ion set "+str(vx)+" "+str(vy)+" "+str(vz))          # Set its velocity
  vmag = math.sqrt(vx*vx + vy*vy + vz*vz)

  step = 0.1                                                                # advance the ion until it has a neighbor
  dx = -vx*step/vz
  dy = -vy*step/vz
  dz = -step
  run(lmp, 0)
  eng = lmp.extract_compute('ionpe',0,0)
  while eng==0:
    lmp.command("displace_atoms ion move "+str(dx)+" "+str(dy)+" "+str(dz))  
    run(lmp, 0)
    eng = lmp.extract_compute('ionpe',0,0)

  step = -0.01                                                              # back off the ion until it has no neighbor
  dx = -vx*step/vz
  dy = -vy*step/vz
  dz = -step
  run(lmp, 0)
  while eng!=0:
    lmp.command("displace_atoms ion move "+str(dx)+" "+str(dy)+" "+str(dz))  
    run(lmp, 0)
    eng = lmp.extract_compute('ionpe',0,0)

  step = 0.001                                                              # advance the ion until it has a neighbor
  dx = -vx*step/vz
  dy = -vy*step/vz
  dz = -step
  run(lmp, 0)
  while eng==0:
    lmp.command("displace_atoms ion move "+str(dx)+" "+str(dy)+" "+str(dz))  
    run(lmp, 0)
    eng = lmp.extract_compute('ionpe',0,0)

  lmp.command("displace_atoms ion move "+str(-dx)+" "+str(-dy)+" "+str(-dz)) # back off one step

######################################
def dump(lmp, file):
    lmp.commands_list(["write_data "+file])
    return file

######################################
def timeSeed():
    from time import clock_gettime
    return int(clock_gettime(0) * 1e5 % 1e9)

def dateTime():
  from time import localtime, asctime
  return asctime(localtime())

######################################
def thermalImpact(lmp, log):
  run(lmp, 0)
  eng = lmp.extract_compute('ionpe',0,0)
  vz = lmp.extract_compute('vIon',0,0)

  t=0
  steps=50
  while vz<0 and t<20000:
    run(lmp, steps)
    run(lmp, 0)
    eng = lmp.extract_compute('ionpe',0,0)
    vz = lmp.extract_compute('vIon',0,0)
    t+=steps

  if eng!=0: # make sure it's not going to come off
    run(lmp, 4000)
    run(lmp, 0)
    eng = lmp.extract_compute('ionpe',0,0)
    vz = lmp.extract_compute('vIon',0,0)

  if eng==0 and vz>0:
    log.write("* Species scattered. Ion deleted.\n")
    lmp.command("delete_atoms group ion")
    return False

  if t>=20000:
    log.write("* Species too slow. Ion deleted.\n")
    lmp.command("delete_atoms group ion")
    return False

  else:
    log.write("* Species bound. Total atoms "+str(lmp.get_natoms())+".\n")
    return True

######################################
def productSweep(lmp, log):
  from collections import defaultdict

  delete = False
  lmp.command("compute clu unfixed cluster/atom 3")
  lmp.command("compute all_id all property/atom id type")
  run(lmp,0)
  N=lmp.get_natoms()
  idlist = lmp.extract_compute('all_id',1,2)[0:N]
  typmap={}
  for i in range(0,len(idlist)):
    typmap[int(idlist[i][0])]=int(idlist[i][1])
  clulist = lmp.extract_compute('clu',1,1)[0:N]
  clumap=defaultdict(list)
  for c in range(0,len(clulist)):
    clumap[int(clulist[c])].append(int(idlist[c][0]))
  for key in clumap.keys():
    if len(clumap[key])<32:
      delete=True
      cluform={}
      cluform["Si"]=0
      cluform["F"]=0
      cluform["Cl"]=0
      for v in clumap[key]:
        cluform[atomtype[typmap[v]]]+=1
        lmp.command("group etchprod id "+str(v))
      cluname =""
      if cluform["Si"]>0: cluname+="Si"
      if cluform["Si"]>1: cluname+=str(cluform["Si"])
      if cluform["F"]>0: cluname+="F"
      if cluform["F"]>1: cluname+=str(cluform["F"])
      if cluform["Cl"]>0: cluname+="Cl"
      if cluform["Cl"]>1: cluname+=str(cluform["Cl"])
      log.write("* Deleting etch cluster "+cluname+"\n")
  if delete:
      lmp.command("delete_atoms group etchprod")
  return delete

