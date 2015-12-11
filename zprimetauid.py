# Example that shows how to make jets from PFA (reconstructed particles) using  LCIO format
# This example makes invariant mass of 2 reconstructed jets and makes Z peak from Z-> tau tau
# This is a simple example that uses Durham algorithm to force exactly 2 jets
# Use Z to bbar data 
#  Get 4 file in 2 thread (see HepSim manual):
#  > hs-get gev250ee_pythia6_zpole_tautau%rfull001 gev250ee_pythia6_zpole_tautau 2 50 
#
# author: S.Chekanov (ANL)


from org.lcsim.event import *
from org.lcsim.util import *
from hep.physics.vec import BasicHepLorentzVector,HepLorentzVector 
from java.util import *
from java.io import *
from org.lcsim.lcio import LCIOReader
from hep.io.sio import SIOReader
from hep.lcio.implementation.sio import SIOLCReader
from hep.lcio.implementation.io import LCFactory
from hep.lcio.event import *
from hep.lcio.io import *
from jhplot import *  # import graphics
from hephysics.particle import LParticle 
from hep.physics.jet import FixNumberOfJetsFinder
from math import *

# make list of files..
import glob
files=glob.glob("/users/kotwal/fullSim/fullsim/output/tev1mm_pythia6_zprime1tev_tautau_0001_pandora.slcio")
factory = LCFactory.getInstance()
reader = factory.createLCReader()
reader.open(files)

#====================================Modules=========================================#
#------------------------------------------------------------------------------
def delR(v1,v2): #angular distance where v1 and v2 are in LParticle format
   d = sqrt( (v1.phi()-v2.phi())**2 + (v1.pseudoRapidity()-v2.pseudoRapidity())**2 )
   return d
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def COS(v1,v2): #rationale behind the recursion is v1 and v2 atleast would differ by name always, so this will return cos and not dot
   c = (v1.x()*v2.x())  + (v1.y()*v2.y())  + (v1.z()*v2.z())
   if (v1 != v2):
      c = c*1./sqrt(COS(v1,v1)*COS(v2,v2))
   return c
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def Et(v): #transverse energy with v in LParticle format #Although the built-in function works fine, checked
   Et_sq = (v.e())**2 - (v.z())**2
   et = sqrt(abs(Et_sq))
   return et
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def track_check(v): #check if v (LParticle) a track or not.
   if ((v.getCharge() != 0)):# or (v.getName() == '2112')):
      tc = True
   else:
      tc = False
   return tc
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def prong(V):
   n = 0
   for i in range(len(V)):
      if (track_check(V[i])): n = n+1
   return n
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# jet V being an list of its members in LParticle format 
#and jet being the combined jet particle in LParticle format
def fcent(V, jet): #Central energy fraction
   E_deno = 0
   Et_num = 0        
   global core
   for i in range(len(V)):
      if (delR(V[i], jet) < core):
         E_deno = E_deno + V[i].e()
         if ((delR(V[i], jet) < (core*0.5))):
            Et_num = Et_num + Et(V[i])
            #1#else: print V[i].getName()
   if (E_deno == 0):
      c = "NAN"
   else:
      c = (Et_num*1./E_deno)
   return c
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# jet V being an list of its members in LParticle format 
#and jet being the combined jet particle in LParticle format
def ftrack(V, jet): #Leading track momentum fraction
   Et_deno = 0
   pt_num = 0
   global core
   for i in range(len(V)):
      if (delR(V[i], jet) < core):
         Et_deno = Et_deno + Et(V[i])
         if (track_check(V[i])):
            if (V[i].perp() > pt_num):
               pt_num = V[i].perp()
   if (Et_deno == 0):
      t = "NAN"
   else:
      t = (pt_num*1./Et_deno)
   return t
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# jet V being an list of its members in LParticle format 
#and jet being the combined jet particle in LParticle format
def Rtrack(V, jet): #Track Radius
   num = 0
   den = 0
   for i in range(len(V)):
      if (track_check(V[i])):
         r = delR(V[i], jet)
         if (r < 0.4): #core + isolation region
            num = num + r*V[i].perp()
            den = den + V[i].perp()
   if (den == 0):
      rt = "NAN"
   else:
      rt = num*1./den
   return rt
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# jet V being an list of its members in LParticle format 
#and jet being the combined jet particle in LParticle format
def Niso_track(V, jet): #Number of tracks in the isolated region
   n = 0
   global core
   for i in range(len(V)):
      if (track_check(V[i])):
         r = delR(V[i], jet)
         if ((r < 0.4) and (r > core)):
            n = n + 1
   return n
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# jet V being an list of its members in LParticle format 
#and jet being the combined jet particle in LParticle format
def max_delR(V, jet): #the max delR of a track in  core region from the jet direction
   del_r = 0
   global core
   for i in range(len(V)):
      if (track_check(V[i])):
         r = delR(V[i], jet)
         if ((r < core) and (r > del_r)):
            del_r = r
   return del_r
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# jet V being an list of its members in LParticle format 
#and jet being the combined jet particle in LParticle format
def mtrack(V, jet): #Inv. mass using all tracks in iso and core region, assuming each with pi0 mass
   mpi0=0.139
   jetsum = LParticle()
   for i in range(len(V)):
      if (track_check(V[i])):
         r = delR(V[i], jet)
         if (r < 0.4):
            E = sqrt(mpi0**2 + V[i].x()**2 + V[i].y()**2 + V[i].z()**2)
            W = LParticle(V[i].x(), V[i].y(), V[i].z(), E)
            jetsum.add(W)
###            print jetsum.calcMass()
   return jetsum.calcMass()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def PTCUT(jet):
   Flag=False
   global TRACK_PT_CUT
   pT_max = 0.0
   hard = 0
   for ss in range(len(jet)):
      if (pT_max < jet[ss].perp()):
         pT_max = jet[ss].perp()
   if (pT_max > TRACK_PT_CUT):
      Flag=True
   return Flag
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#def pi0rec(,jet)
#------------------------------------------------------------------------------
#====================================================================================#
##h1=H1D("M(jj)",30,30,120)
hfcent=H1D("fcent",30,0,1.0)
hfcentT=H1D("fcent_truth",30,0,1.0)
hftrack=H1D("ftrack",50,0,1.0)
hftrackT=H1D("ftrack_truth",50,0,1.0)
hRtrack=H1D("Rtrack",50,0,0.4)
hRtrackT=H1D("Rtrack_truth",50,0,0.4)
hNiso_track=H1D("Niso_track", 20, 0, 10)
hNiso_trackT=H1D("Niso_track_truth", 20, 0, 10)
hMaxDeltaR=H1D("MaxDeltaR", 50, 0, 0.2)
hMaxDeltaRT=H1D("MaxDeltaR_truth", 50, 0, 0.2)
hMtrack=H1D("Mtrack", 50, 0, 2.5)#10)
hMtrackT=H1D("Mtrack_truth", 50, 0, 10)
hpT=H1D("pT", 100, 0, 200)
hpT_uncut=H1D("pT_uncut", 100, 0, 200)
hnjettracks=H1D("no. of tracks/jet after cut", 20, 0, 10)
hnjettracks_uncut=H1D("no. of tracks/jet without cut", 20, 0, 10)
hneventtracks=H1D("no. of tracks/event after cut", 100, 0, 50)
hneventtracks_uncut=H1D("no. of tracks/event without cut", 100, 0, 50)
hjetpT=H1D("pT of jets in the event after cut", 100, 0, 200)
hjetpT_uncut=H1D("pT of jets in the event without cut", 100, 0, 200)
hleptaud=H1D("cosine angle of recon. jet from had dking tau direction", 50, -1, 1.01)
# h1.doc()  # look at API

# jet description
# http://java.freehep.org/freehep-physics/apidocs/hep/physics/jet/package-summary.html
# this is Jade
from hep.physics.jet import DurhamJetFinder,FixNumberOfJetsFinder
fjet=FixNumberOfJetsFinder(2) # request 2 jets. The Jade algorithm at work 
#fjet=DurhamJetFinder(0.05) 
core=0.1 #(size of the core of a jet)
pl = 0
nEvent=0
neventtracks=0
neventtracks_uncut=0
TRACK_PT_CUT = 15.0 #hardest track in a jet candidate should be more than 15 GeV to be considered
lepdk_taueventcount = 0
taulepinivis1 = 0
iter1 = 0
finalchargedpion=0
taunum=0
lepdk = 0
hadk=0
event_rej = 0
one_lep=0
taucount=0
onegamma=0
N = input("enter the number of prong(s): ") #no of prong(s) decay
tau3p = 0
jetcount=0
realjetcount=0
tau1p = 0
while(1):
     evt=reader.readNextEvent()
     if (evt == None): break
#     print nEvent ,"ends"
     nEvent=nEvent+1
#     if (nEvent==2): break #print 
     print "# Event: ",nEvent
     strVec = evt.getCollectionNames()
     if nEvent == 1:
            for col in  strVec:
                           print col        
#.........rejecting events with tau_lep-invis ....................#
     tau_lepdk = 0
     col1 = evt.getCollection("MCParticle")
     nMC = col1.getNumberOfElements()
     turthparticles=ArrayList()
     print "nMC: ", nMC, evt
     for ss in range(nMC):
          par=col1.getElementAt(ss)
          pdg=par.getPDG()
          flaglep=False
          photorecoil=False
        #  if(abs(pdg) == 2112):
        #     parents=par.getParents()
        #     print "neutron with %i parents:" %(len(parents))
        #     for kk in range(len(parents)):
        #        print parents[kk].getPDG(), "parent"
          if(abs(pdg) == 15):
             tau_d = par.getDaughters()
             l = len(tau_d)
             if (l > 1): #when tau starts decaying
                for jj in range(len(tau_d)):
                    daughter = tau_d[jj].getPDG()
                    if ((abs(daughter)==15)): 
                       photorecoil=True
        #               print "photon recoil"
                       break # tau just emits photon
                    if ((abs(daughter)>=11) and (abs(daughter)<=14)):
                       flaglep=True #looking for lepton decay
                       one_lep = one_lep + 1
                       tau_lepdk = tau_lepdk + 1
                       break      
                if ((photorecoil==False)):
                   taucount=taucount+1
                   print taucount, nEvent
                   for xx in range(len(tau_d)):
                      print tau_d[xx].getPDG()
       #            else: print "D: ", daughter
                if ((flaglep==False)and(photorecoil==False)):
                  # taucount=taucount+1
                 #  print taucount, nEvent
#                   print par.getMomentum()[0], par.getMomentum()[1], par.getMomentum()[2]
                   leptau=LParticle('tauhad-dir')
                   leptau.setPxPyPz(par.getMomentum()[0], par.getMomentum()[1], par.getMomentum()[2])
                   truejet=[]
                   for aa in range(len(tau_d)):
                      daughtau=tau_d[aa]
                      if ((abs(daughtau.getPDG()))!=16):
                         truejetmem=LParticle(str(daughtau.getPDG()))
                         truejetmem.setPxPyPzE(daughtau.getMomentum()[0], daughtau.getMomentum()[1], daughtau.getMomentum()[2], daughtau.getEnergy())
                         truejetmem.setCharge(daughtau.getCharge())
                         truejet.append(truejetmem)
#                   print leptau
#!                   if ((PTCUT(truejet)) and (prong(truejet)==N)):
#!                      FCENT_t = fcent(truejet, leptau)
#!                      FTRACK_t = ftrack(truejet, leptau)
#!                      RTRACK_t = Rtrack(truejet, leptau)
#!                      NISO_t = Niso_track(truejet, leptau)
#!                      MAXDELR_t = max_delR(truejet, leptau)
#!                      MTRACK_t = mtrack(truejet, leptau)
#!                      if (FCENT_t != "NAN"):
#!                         hfcentT.fill(FCENT_t)
#!                      if (FTRACK_t != "NAN"):
#!                         hftrackT.fill(FTRACK_t)
#!                      if (RTRACK_t != "NAN"):
#!                         hRtrackT.fill(RTRACK_t)
#!                      hNiso_trackT.fill(NISO_t)
#!                      hMaxDeltaRT.fill(MAXDELR_t)
#!                      hMtrackT.fill(MTRACK_t)
#                if (flaglep):
#                   tau_lepdk = tau_lepdk + 1
#                   break
#     print tau_lepdk
     if (tau_lepdk==2):
        event_rej = event_rej + 1
        continue #skipping the event as both the taus dk leptonically
     else:
        realjetcount=realjetcount+1
#-----------------------------------------------------------------#










#     print nEvent
     col = evt.getCollection("PandoraPFOCollection")
     nPFA = col.getNumberOfElements()

     particles=ArrayList() # list of particles
     particleCharge=[]  
     #print "num: ", nPFA
     PDGlist = []
     for i in range(nPFA):
          njettracks=0
          njettracks_uncut=0
          pa=col.getElementAt(i)
          charge=pa.getCharge()
          p4=pa.getMomentum()  
          ee=pa.getEnergy() 
          typep=pa.getType() # type of PF 
        #  print typep
          PDGlist.append(typep)
          p=BasicHepLorentzVector(ee,p4[0],p4[1],p4[2])
          particles.add(p) # add particle to the list 
          particleCharge.append(charge)
#          if (abs(typep) == 211):
#             print p
#          print "PID: ",typep, "\t", p, "\t", type(particles)
     #print particles.size()
     if (particles.size()>1): # ask for >1 particles
       fjet.setEvent(particles)
##       print "i:", particles
##       print particleCharge
      ## print "Jet - ", fjet
       # print "Nr of jets=", fjet.njets()  
       alljets=[] # make a new list with jets 
#       print "# of jets ", fjet.njets()
       for i in range(fjet.njets()):
		 # print "Jet=",i," Nr if particles in jet =",fjet.nParticlesPerJet(i) 
         pjet=fjet.jet(i)   # make  HepLorentzVector
###         print fjet.particlesInJet(i)
###         print "jet", pjet
         ee=pjet.t()  # energy
         p3=pjet.v3() # 3-vector  
         #print " ",ee,p3
         pl_uncut=LParticle(p3.x(),p3.y(),p3.z(),ee)
         hjetpT_uncut.fill(pl_uncut.perp())
#================MAW===========================================================
         JetIds = [] #storing the PID of particles in the jet
         JETset=[] #converting fjet into LParticle array
         for j in range(fjet.nParticlesPerJet(i)):   
            for s in range(len(particleCharge)):
               if (fjet.particlesInJet(i)[j] == particles[s]):
                  particles.remove(particles[s])
                  CHARGE=particleCharge[s]
                  particleCharge.pop(s)
                  pid_j =  PDGlist[s]
                  PDGlist.pop(s)
                  if (CHARGE != 0.0):
                     njettracks_uncut=njettracks_uncut+1
                  break
            JetIds.append(pid_j)
            mom = fjet.particlesInJet(i)[j].v3()
            JETelement=LParticle(str(pid_j))
            JETelement.setPxPyPzE(mom.x(), mom.y(), mom.z(), fjet.particlesInJet(i)[j].t())
            JETelement.setCharge(CHARGE)
            JETset.append(JETelement)
         #1#print i, ": ", JetIds
         if ((len(JetIds) == 1) and (JetIds[0] == 22)):
            onegamma=onegamma+1
#~~~~~~~~~~~~~~~~~~~~~~ PT CUT ON HARDEST TRACK IN A TAU JET CANDIDATE ~~~~~~~~
#         pT_max = 0.0
#         hard = 0
#!         for ss in range(len(JETset)):
#            if (pT_max < JETset[ss].perp()):
#               pT_max = JETset[ss].perp()
#               hard = ss
#!            hpT_uncut.fill(JETset[ss].perp())
         pl=LParticle(p3.x(),p3.y(),p3.z(),ee) # convert to a class with convinient kinematic methods
         if ((tau_lepdk!=0) and (COS(pl,leptau) < 0.98)): continue
         else:jetcount=jetcount+1
         if (prong(JETset)==3):
            tau3p = tau3p + 1
         if (prong(JETset)==1):
            tau1p = tau1p + 1
         if ((PTCUT(JETset)) and (prong(JETset)==N)):
            for b in range(len(JETset)):
               hpT.fill(JETset[b].perp())
               if (JETset[b].getCharge() != 0.0):
                  njettracks=njettracks+1
#            print "prong: ", njettracks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##         print fjet.particlesInJet(i)
#         print fjet.nParticlesPerJet(i)
##         print JETset 
##         print JETset[0].getCharge(), "\t", JETset[1].getCharge()
###            pl=LParticle(p3.x(),p3.y(),p3.z(),ee) # convert to a class with convinient kinematic methods  
###            if ((tau_lepdk!=0) and (COS(pl,leptau) < 0)):
###                continue
###            else:
###                jetcount=jetcount+1
#!               COSleptau = COS(pl,leptau)

#               print COSleptau
#!               hleptaud.fill(COSleptau)
#!            hjetpT.fill(pl.perp())
            #1#print JETset[hard].getName() , "center"
         #print type(fcent(JETset, pl)), type(pl.calcMass())
#!            hnjettracks.fill(njettracks)
#!            hnjettracks_uncut.fill(njettracks_uncut)
#!            FCENT = fcent(JETset, pl)
#!            FTRACK = ftrack(JETset, pl)
#!            RTRACK = Rtrack(JETset, pl)
#!            NISO = Niso_track(JETset, pl)
#!            MAXDELR = max_delR(JETset, pl)
#!            MTRACK = mtrack(JETset, pl)
##         print FCENT, "\t", FTRACK
#         if (FCENT == -100):
#            print "pjet ", pjet

#!            if (FCENT != "NAN"):
#!               hfcent.fill(FCENT)
#!            if (FTRACK != "NAN"):
#!               hftrack.fill(FTRACK)
#!            if (RTRACK != "NAN"):
#!               hRtrack.fill(RTRACK)
#!            hNiso_track.fill(NISO)
#!            hMaxDeltaR.fill(MAXDELR)
#!            hMtrack.fill(MTRACK)
#            if ((RTRACK != 'NAN') and (FTRACK == 'NAN')):
#               print "*****************"
#               for i in JETset:
#                  if (track_check(i)):
#                     print "defaulter: ",(delR(i, pl)), i.perp(), RTRACK 
#!     neventtracks=neventtracks+njettracks
#!     neventtracks_uncut=neventtracks_uncut+njettracks_uncut
#!     hneventtracks.fill(neventtracks)
#!     hneventtracks_uncut.fill(neventtracks_uncut)
#     iter1 = iter1+1
#     print "event-", iter1, " ends"
#     nEvent=nEvent+1
#==============================================================================
reader.close();
print tau1p, "1p<                  >3p", tau3p, "----------------", jetcount, "-----------", realjetcount, "-----------", taucount
print "events:%i event_rej:%i one_lep:%i iter:%i"%(nEvent, event_rej, one_lep, iter1)
exit()
from java.awt import Color
c1=HPlot("Canvas", 1000, 1200, 2, 3)
# c1.doc() # look at API
c1.visible(1)
c1.cd(1,1)
c1.setMarginLeft(100)
c1.setNameX("fcent")
c1.setNameY("Events")
c1.setRangeX(0,1.01)
hfcent.setColor(Color.red)
hfcent.setTitle("Det_level %i jets"%(hfcent.entries()))
hfcentT.setColor(Color.green)
hfcentT.setTitle("Truth level %i jets"%(hfcentT.entries()))
c1.draw(hfcent)
c1.draw(hfcentT)

c1.cd(1,2)
c1.setMarginLeft(100)
c1.setNameX("ftrack")
c1.setNameY("Events")
c1.setRangeX(0,1.01)
hftrack.setColor(Color.red)
hftrack.setTitle("Det_level %i jets"%(hftrack.entries()))
hftrackT.setColor(Color.green)
hftrackT.setTitle("Truth level %i jets"%(hftrackT.entries()))
c1.draw(hftrack)
c1.draw(hftrackT)

c1.cd(1,3)
c1.setMarginLeft(100)
c1.setNameX("rtrack")
c1.setNameY("Events")
c1.setRangeX(0,0.41)
hRtrack.setColor(Color.red)
hRtrack.setTitle("Det_level %i jets"%(hRtrack.entries()))
hRtrackT.setColor(Color.green)
hRtrackT.setTitle("Truth level %i jets"%(hRtrackT.entries()))
c1.draw(hRtrack)
c1.draw(hRtrackT)

c1.cd(2,1)
c1.setMarginLeft(100)
c1.setNameX("Niso_track")
c1.setNameY("Events")
c1.setRangeX(0,10.1)
hNiso_track.setColor(Color.red)
hNiso_track.setTitle("Det_level %i jets"%(hNiso_track.entries()))
hNiso_trackT.setColor(Color.green)
hNiso_trackT.setTitle("Truth level %i jets"%(hNiso_trackT.entries()))
c1.draw(hNiso_track)
c1.draw(hNiso_trackT)

c1.cd(2,2)
c1.setMarginLeft(100)
c1.setNameX("MaxDeltaR")
c1.setNameY("Events")
c1.setRangeX(0,0.21)
hMaxDeltaR.setColor(Color.red)
hMaxDeltaR.setTitle("Det_level %i jets"%(hMaxDeltaR.entries()))
hMaxDeltaRT.setColor(Color.green)
hMaxDeltaRT.setTitle("Truth level %i jets"%(hMaxDeltaRT.entries()))
c1.draw(hMaxDeltaR)
c1.draw(hMaxDeltaRT)

c1.cd(2,3)
c1.setMarginLeft(100)
c1.setNameX("Mtrack")
c1.setNameY("Events")
c1.setRangeX(0, 2.51)#10.1)
hMtrack.setColor(Color.red)
hMtrack.setTitle("Det_level %i jets"%(hMtrack.entries()))
hMtrackT.setColor(Color.green)
hMtrackT.setTitle("Truth level %i jets"%(hMtrackT.entries()))
c1.draw(hMtrack)
c1.draw(hMtrackT)


c1.export("tauId_%i_prong_core_%f_Zprimedata.pdf"%(N,core))




