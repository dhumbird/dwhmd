#!/usr/bin/python3

import sys
#from lammps import lammps

def ionBombard(args):
	from lammps import lammps
	import dwhlmp as md
	
	datfile="temp_000000.dat"
	logfile="md.log"
	e_ion=0
	Tion=300
	T=300
	runs=1
	seeds=[]
	cont=False
	
	i=0
	while i < len(args):
		if args[i].find(".dat")>-1: datfile=args[i]
		elif (args[i]=="-e"): 
			i+=1
			e_ion=float(args[i])
		elif (args[i]=="-T"):
			i+=1
			T=float(args[i])
		elif (args[i]=="-Tion"):
			i+=1
			Tion=float(args[i])
		elif (args[i]=="-runs"):
			i+=1
			runs=int(args[i])
		elif (args[i]=="-log"):
			i+=1
			logfile=args[i]
		elif (args[i]=="-cont"):
			cont = True
		elif (args[i]=="-seeds"):
			i+=1
			for j in range(0,runs):
				seeds.append(int(args[i]))
				i+=1
		else: print("ionBombard: bad argument ",args[i])
		i+=1
	filebase = datfile[0:datfile.find("_")+1]
	A=datfile
	if len(seeds)==0:
		seed=md.timeSeed()
		for run in range(0,runs):
			seeds.append(seed+run)
	if cont:
		maxrun=0
		with open("md.log") as log:
			for line in log:
				if line.find("---run") > -1:
					run= int((line[line.find("run")+4:]).replace('-','').strip())
					maxrun=max(run,maxrun)
		with open("md.log",'a',1) as log:
			for run in range(maxrun+1,maxrun+runs+1):
				lmp = lammps(cmdargs=["-echo", "screen"])
				log.write("------------"+md.dateTime()+"--------run "+str(run)+"-----------------------\n")
				md.init(lmp,A,log)
#				seed=md.timeSeed()
				seed=seeds[run-maxrun-1]
				md.thermo(lmp, seed, T)
				md.addion(lmp, seed, log)
				bound = md.thermalImpact(lmp,log)
				sput=md.productSweep(lmp,log)
				if bound or sput:
					file = filebase+str(run).rjust(6,"0")+".dat"
					A = md.dump(lmp, file)
				lmp.close()
				run+=1
			log.write("---------------------------------------------------------------------------\n")
			log.write(md.dateTime()+" Requested runs ("+str(runs)+") completed.\n")   	
	else:
		with open("md.log",'w',1) as log:
			for run in range(1,runs+1):
				lmp = lammps(cmdargs=["-echo", "screen"])
				log.write("------------"+md.dateTime()+"--------run "+str(run)+"-----------------------\n")
				md.init(lmp,A,log)
#				seed=md.timeSeed()
				seed=seeds[run-1]
				md.thermo(lmp, seed, T)
				md.addion(lmp, seed, log)
				bound = md.thermalImpact(lmp,log)
				sput=md.productSweep(lmp,log)
				if bound or sput:
					file = filebase+str(run).rjust(6,"0")+".dat"
					A = md.dump(lmp, file)
				lmp.close()
				run+=1
			log.write("---------------------------------------------------------------------------\n")
			log.write(md.dateTime()+" Requested runs ("+str(runs)+") completed.\n")   	

#########

def EasyRun(args):
	from lammps import lammps
	import dwhlmp as md
	
	datfile="temp_000000.dat"
	logfile="md.log"
	runtime=1
	To=0
	dt=0.001
	bhb_dt=0
	cfgevery=0
	track=-1
	tso=False
	printe=False
	dump_end=False
	quiet=False
	sput=1
	quench=False 


	i=0
	while i < len(args):
		if args[i].find(".dat")>-1: datfile=args[i]
		elif (args[i]=="-time"): 
			i+=1
			runtime=float(args[i])
		elif (args[i]=="-T"):
			i+=1
			To=float(args[i])
		elif (args[i]=="-tso"):
			tso=True
		elif (args[i]=="-cfgevery"):
			i+=1
			cfgevery=int(args[i])
		elif (args[i]=="-log"):
			i+=1
			logfile=args[i]
		elif (args[i]=="-dt"):
			i+=1
			dt=float(args[i])
		elif (args[i]=="-bhb_dt"):
			i+=1
			bhb_dt=float(args[i])
		elif (args[i]=="-printe"):
			printe=True
		elif (args[i]=="-track"):
			i+=1
			track=int(args[i])
		elif (args[i]=="-dump_end"):
			dump_end=True
		elif (args[i]=="-quiet"):
			quiet=True
		elif (args[i]=="-nosput"):
			sput=0
		elif (args[i]=="-sweep"):
			sput=2
		elif (args[i]=="-quench"):
			i+=1
			quench=True
			To=float(args[i])
		else: print("run: bad argument ",args[i])
		i+=1

		filebase = datfile[0:datfile.find("_")+1]

		steps = int(runtime/dt)
		with open(logfile,'w',1) as log:
			lmp = lammps(cmdargs=["-echo", "screen"])
			log.write("------------"+md.dateTime()+"--------run -----------------------\n")
			md.init(lmp,datfile,log)
			md.run(lmp,steps)
#			md.thermo(lmp, seed, T)
#			md.addion(lmp, seed, log)
#			bound = md.thermalImpact(lmp,log)
#			sput=md.productSweep(lmp,log)
#			if bound or sput:
			file = filebase+str(steps).rjust(6,"0")+".dat"
			A = md.dump(lmp, file)
#			lmp.close()
			log.write("---------------------------------------------------------------------------\n")
#			log.write(md.dateTime()+" Requested runs ("+str(runs)+") completed.\n")   	

#  config cfg(scfg);
#  if (quench)
#    cfg.RunQuenched(runtime, cfgevery, dt, printe, dump_end, track, quiet, tso, sput, To);
#  else{
#    if (To==0)
#      cfg.Run(runtime, cfgevery, dt, printe, dump_end, track, quiet, tso, sput);
#    else
#      cfg.Thermo(To, bhb_dt, cfgevery, dt, printe, dump_end, tso,runtime);


###################################################################
###################### COMMAND-LINE HANDLER #######################
###################################################################
if len(sys.argv) <= 2:
	print ("supply some more arguments.")

else: 
	arg=sys.argv[1]
	if (arg=="-run"): EasyRun(sys.argv[2:])
  # else if (arg=="-addfix") MainAddFixed(argc,argv);
  # else if (arg=="-reset") MainReset(argc, argv);
  # else if (arg=="-delete") MainDelete(argc, argv);
  # else if (arg=="-seto") MainShiftOrigin(argc,argv);
  # else if (arg=="-rand") MainRandomize(argc,argv);
  # else if (arg=="-init") MainInit(argc,argv);
  # else if (arg=="-edep") MainEnergyDep(argc, argv);
  # //else if (arg=="-sat") MainSat(argc,argv);
  # //else if (arg=="-pass") MainPass(argc,argv);
	elif (arg=="-addion"): print("MainAddion ", sys.argv[2:])
  # //else if (arg=="-path") MainPath(argc,argv);
	elif (arg=="-bomb"): ionBombard(sys.argv[2:])
  # else if (arg=="-mixion") MainMixIon(argc,argv);
  # else if (arg=="-cgmin") MainCGmin(argc,argv);
  # else if (arg=="-sweep") MainMisc(argc,argv);
	else: print("bad first argument ",arg)


