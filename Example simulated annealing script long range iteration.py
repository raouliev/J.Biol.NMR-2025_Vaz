xplor.requireVersion("2.34")

# slow cooling protocol in torsion angle space for ASR.

xplor.parseArguments()
import os
import protocol
command = xplor.command

# define directory and number of structures
str_dir = "../structures/anneal_ASR_longrange_roundX/"
numberOfStructures=100

# general setup
outFilename = str_dir+"SCRIPT_STRUCTURE.pdb"
os.system("mkdir "+str_dir)
#starting structure, extended/template/etc.
protocol.loadPDB("../pdb/PRE_template.pdb",deleteUnknownAtoms=True)

#ASR specific - fixes retinal in all-trans conformation (later in script)
ret_sel = "resid 210 and (name CB or name HB1 or name HB2 or name CG or name HG1 or name HG2 or name CD or name HD1 or name HD2 or name CE or name HE1 or name HE2 or name N16 or name H16 or name C15 or name H15 or name C14 or name H14 or name C13 or name C20 or name H20A or name H20B or name H20C or name C12 or name H12 or name C11 or name H11 or name C10 or name H10 or name C9 or name C19 or name H19A or name H19B or name H19C or name C8 or name H8 or name C7 or name H7 or name C6 or name C5 or name C18 or name H18A or name H18B or name H18C or name C4 or name H4A or name H4B or name C3 or name H3A or name H3B or name C2 or name H2A or name H2B or name C1 or name C16 or name H16A or name H16B or name H16C or name C17 or name H17A or name H17B or name H17C)"
#everything other than retinal
nonret_sel = "not ("+ret_sel+")"

from potList import PotList
potList = PotList()

from simulationTools import MultRamp, StaticRamp, InitialParams
rampedParams=[]
highTempParams=[]

#only needed if there is a reference structure, backbone difference potential
#from posDiffPotTools import create_PosDiffPot
#refRMSD = create_PosDiffPot("refRMSD","name CA or name C or name N",
#                            pdbFile="../pdb/reference_structure.pdb",
#                            cmpSel="not name H*")

# set up NOE potential
noe=PotList('noe')
potList.append(noe)
from noePotTools import create_NOEPot
#distance restraints included here: h-bonds, PRE's, disabmiguated restraints
for (name,scale,file) in [('longrange_unambig_roundX',10,"../restraints/longrate_unambig_roundX.tbl"),
                          ('PREs',10,"../restraints/PREs.tbl"),
                          ('hbond',10,"../restraints/hbond_regular.tbl")
                          ]:
    pot = create_NOEPot(name,file)
    pot.setScale(scale)
    pot.setPotType("soft")
    noe.append(pot)
rampedParams.append( MultRamp(2,30, "noe.setScale( VALUE )") )

# set up dihedral angles
from xplorPot import XplorPot
dihedralRestraintFilename="../restraints/dihedral_nosegid.tbl"
protocol.initDihedrals(dihedralRestraintFilename, scale=20)
potList.append( XplorPot('CDIH') )
highTempParams.append( StaticRamp("potList['CDIH'].setScale(20)") )
rampedParams.append( StaticRamp("potList['CDIH'].setScale(200)") )
# set custom values of threshold values for violation calculation
potList['CDIH'].setThreshold( 5 )

# hbda - distance/angle bb hbond term
protocol.initHBDA('../restraints/hbda_asr.tbl')
potList.append( XplorPot('HBDA') )
potList['HBDA'].setScale(200)

# setup parameters for atom-atom repulsive term. (van der Waals-like term)
potList.append( XplorPot('VDW') )
rampedParams.append( StaticRamp("protocol.initNBond()") )
rampedParams.append( MultRamp(.004,4,
                              "command('param nbonds rcon VALUE end end')") )
# nonbonded interaction only between CA atoms
highTempParams.append( StaticRamp("""protocol.initNBond(cutnb=100,
                                                        rcon=0.004,
                                                        tolerance=45,
                                                        repel=1.2,
                                                        onlyCA=1)""") )

potList.append( XplorPot("BOND") )
potList.append( XplorPot("ANGL") )
potList['ANGL'].setThreshold( 5 )
rampedParams.append( MultRamp(0.4,1,"potList['ANGL'].setScale(VALUE)") )
potList.append( XplorPot("IMPR") )
potList['IMPR'].setThreshold( 5 )
rampedParams.append( MultRamp(0.1,1,"potList['IMPR'].setScale(VALUE)") )

# Give atoms uniform weights, except for the anisotropy axis
protocol.massSetup()

# IVM setup
from ivm import IVM
dyn = IVM()
#retinal always remains fixed
dyn.group(ret_sel)
protocol.torsionTopology(dyn)

minc = IVM()
minc.group(ret_sel)
protocol.initMinimize(minc)
protocol.cartesianTopology(minc)

# object which performs simulated annealing
from simulationTools import AnnealIVM
init_t  = 20000.     # Need high temp and slow annealing to converge. Can be raised when first introducing PRE's if helices don't form typical rhodopsin bundle
mid_t = 1000.
hi_t = 10000.
cool = AnnealIVM(initTemp =init_t,
                 finalTemp=1000,
                 tempStep =250.0,
                 ivm=dyn,
                 rampedParams = rampedParams)
cart_cool = AnnealIVM(initTemp =mid_t,
	              finalTemp=25,
		      tempStep =25.0,
                      ivm=minc,
                      rampedParams = rampedParams)

def calcOneStructure(loopInfo):
    """ this function calculates a single structure, performs analysis on the
    structure, and then writes out a pdb file, with remarks.
    """

#randomization of torsion angles only needed when introducing PRE's, after secondary structure is formed (see fig. 4)
    # generate a new structure with randomized torsion angles
#    from monteCarlo import randomizeTorsions
#    randomizeTorsions(dyn,sel=nonret_sel)
#    protocol.fixupCovalentGeom(sel=nonret_sel,maxIters=100,useVDW=1,maxViols=3)

    # initialize parameters for high temp dynamics
    InitialParams( rampedParams )
    InitialParams( highTempParams )

    # high temp dynamics
    protocol.initDynamics(dyn,
                          potList=potList, # potential terms to use
                          bathTemp=hi_t,
                          initVelocities=1,
                          numSteps=144000,   # whichever comes first
                          printInterval=1000)
    dyn.setETolerance( init_t/1000 )  #used to det. stepsize. default: t/1000 
    dyn.run()

    # torsion angle annealing
    InitialParams( rampedParams )
    protocol.initDynamics(dyn,
                          potList=potList,
                          numSteps=144000,
                          printInterval=100)
    cool.run()
    # final torsion angle minimization
    protocol.initMinimize(dyn,
                          printInterval=50)
    dyn.run()

    # Cartesian annealing
    protocol.initDynamics(minc,
                          potList=potList,
                          numSteps=144000,
                          printInterval=100)
    cart_cool.run()
    # final all-atom minimization
    protocol.initMinimize(minc,
                          potList=potList,
                          dEPred=0.01)
    minc.run()

    #do analysis and write structure when function returns
    pass

from simulationTools import StructureLoop, FinalParams
StructureLoop(numStructures=numberOfStructures,
              doWriteStructures=True,
              pdbTemplate=outFilename,
              structLoopAction=calcOneStructure,
              genViolationStats=True,
              averageTopFraction=0.5, #report stats on best 50% of structs
              averageContext=FinalParams(rampedParams),
# disabled, no reference structure potential              
#              averageCrossTerms=refRMSD,
              averageSortPots=[potList['BOND'],potList['ANGL'],potList['IMPR'],noe,potList['CDIH']],
              averagePotList=potList).run()
