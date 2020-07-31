#!/bin/env python
#
#  LUND FORMAT
# https://gemc.jlab.org/gemc/html/documentation/generator/lund.html#lund
"""
::: Header :::
   1   2 3 4     5     6 7 8 9 10
 Npart - - -  Bea Pol  - - - -  -

::: Particle :::

 1  2    3           4      5 6   7  8  9  10 11  12 13 14
ID  -  type       Part ID   - -  px py pz   -  -  vx vy vz
      1=active

"""
# all the not useful fields will be filled with 99

from ROOT import TGenPhaseSpace, TLorentzVector, TVector3, TDatabasePDG, TMath

import math
from random import random as rndm
from array import array

dbpdg = TDatabasePDG()

#  motherParticle
PDGmother=223
MotherMass = dbpdg.GetParticle(PDGmother).Mass()

#  daughters
PDGs = [ 211, -211, 111 ]
DaughterMass = [ dbpdg.GetParticle(i).Mass() for i in PDGs ]

# decay pi0
gPDGs = [ 22, 22 ]
gMasses = [ dbpdg.GetParticle(i).Mass() for i in gPDGs ]

# Exclusive collision Settings
# ============================

BeamPolarization = -1 # not useful

BeamEnergy = 10.6
TargetMass = dbpdg.GetParticle( 2212 ).Mass()

beam   = TLorentzVector( 0, 0, BeamEnergy, BeamEnergy )
target = TLorentzVector( 0, 0, 0, TargetMass )
CME    = beam + target

FSPDGs   = [ 11 , 2212, PDGmother ]
FSMasses = [  dbpdg.GetParticle(i).Mass() for i in FSPDGs ]

# electron polar angle cuts
ThetaMin=math.radians(5.0)
ThetaMax=math.radians(35.0)

# vertex 
VxMin = 0.
VxMax = 0.

VyMin = 0.
VyMax = 0.

VzMin = -3.
VzMax =  3.


NEVENTS=100000

## FileName
fname = "phasespace_"
fname += dbpdg.GetParticle(PDGmother).GetName()
for p in PDGs:
  fname += "_" +  dbpdg.GetParticle(p).GetName()
fname += ".lund"

def storeParticle(  pdg, vm, vertex, state ):
  particle = [99]*14
  particle[0] = counter;
  particle[2] = state
  particle[3] = pdg

  particle[6] = vm.Px()
  particle[7] = vm.Py()
  particle[8] = vm.Pz()
  particle[11] = vertex[0] 
  particle[12] = vertex[1] 
  particle[13] = vertex[2] 
  return particle


def writeParticles( outf, particles ):
  # write mother
  for i,p in enumerate( particles ):
    p[0] = i+1
    for h in p:
      outf.write(str(h))
      outf.write(" ")
    outf.write("\n")

def applyCuts( particles ):
  test = True
  pdic = {}
  for p in particles:
    if p[3] in pdic:
      pdic[ p[3]+1000 ] = TVector3( p[6], p[7], p[8] ) 
    else:
      pdic[ p[3] ] = TVector3( p[6], p[7], p[8] ) 


  # check electron angle
  electron = pdic[ 11 ]
  if electron.Theta()  < ThetaMin:
    return False
  if electron.Theta() > ThetaMax:
    return False 

  #check pions 
  pi = pdic[ 211 ]
  if pi.Theta()  < ThetaMin:
    return False
  if pi.Theta() > ThetaMax:
    return False 

  pi = pdic[ -211 ]
  if pi.Theta()  < ThetaMin:
    return False
  if pi.Theta() > ThetaMax:
    return False 

  # check photon
  pi = pdic[ 22 ]
  if pi.Theta() > ThetaMax:
    return False 
  pi = pdic[ 1022 ]
  if pi.Theta() > ThetaMax:
    return False 

  return test

# generate
with open(fname,"w") as fout:

  # generate event
  Event = TGenPhaseSpace()
  Event.SetDecay( CME, len(FSMasses),array('d',FSMasses))

  # loop on events
  for i in range(1,NEVENTS+1):

    # generate event
    counter = 1
    w = Event.Generate()


    # generate vertex
    vertex = (
              VxMin + rndm()*(VxMax - VxMin),
              VyMin + rndm()*(VyMax - VyMin),
              VzMin + rndm()*(VzMax - VzMin)    
             )
   
    # write header 
    header = [99]*10
    header[0] = len(FSPDGs) + len(PDGs) + len(gPDGs)  # daughters + mother + electron + proton
    header[4] = BeamPolarization
    header[6] = w  
    
    # array to store all the particles in the event    
    parts = []

    # the reaction
    for l,j in enumerate(FSPDGs):
      state = 1
      if j == PDGmother:
        state = 0
      parts.append( storeParticle( j, Event.GetDecay(l) , vertex, state ) )

    # prepare decay
    mother   = Event.GetDecay(2)
    Decay = TGenPhaseSpace()
    Decay.SetDecay(mother,len(DaughterMass),array('d',DaughterMass))

    Decay.Generate()
    for j,p in enumerate(PDGs):
      vp = Decay.GetDecay(j)
      state = 1
      if p == 111:
        state = 0
      parts.append( storeParticle( p, vp, vertex, state ) )
      
    # decay the pi0
    Decay.SetDecay( Decay.GetDecay(2),len(gMasses),array('d',gMasses))
    Decay.Generate()
    for j,p in enumerate(gPDGs):
      vp = Decay.GetDecay(j)
      counter += j
      parts.append( storeParticle(  p, vp, vertex, 1  ) )
    

    if applyCuts( parts ) == False:
      continue

    # write the event
    for h in header:
      fout.write(str(h))
      fout.write(" ")
    fout.write("\n")

    writeParticles( fout, parts )

  fout.close()
