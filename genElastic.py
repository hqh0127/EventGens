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

from ROOT import TGenPhaseSpace, TLorentzVector, TVector3, TDatabasePDG, TMath, TH2F

import math
from random import random as rndm
from array import array

dbpdg = TDatabasePDG()

# Exclusive collision Settings
# ============================

BeamPolarization = -1 # not useful

BeamEnergy = 6.5
TargetMass = dbpdg.GetParticle( 2212 ).Mass()

beam   = TLorentzVector( 0, 0, BeamEnergy, BeamEnergy )
target = TLorentzVector( 0, 0, 0, TargetMass )
CME    = beam + target

FSPDGs   = [ 11 , 2212  ]
FSMasses = [  dbpdg.GetParticle(i).Mass() for i in FSPDGs ]

# electron polar angle cuts
ThetaMin=math.radians(4.0)
ThetaMax=math.radians(40.0)

# proton polar angle cuts
PThetaMin=math.radians(35.0)
PThetaMax=math.radians(120.0)

# vertex 
VxMin = -.171
VxMax = -.171

VyMin =  0.1
VyMax =  0.1

VzMin = -2.5
VzMax =  2.5


NEVENTS= 35000

## FileName
fname = "elastic_E6.5_vx-0171vy01"
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

  # check proton angle
  proton = pdic[ 2212 ]
  if proton.Theta()  < PThetaMin:
    return False
  if proton.Theta() > PThetaMax:
    return False 

  return test


h2p = TH2F("hp","", 180,-180,180, 180,-180,180 )
h2t = TH2F("ht","", 180, 0,180, 90,0,45 )
# generate
with open(fname,"w") as fout:

  # generate event
  Event = TGenPhaseSpace()
  Event.SetDecay( CME, len(FSMasses),array('d',FSMasses))

  # loop on events
  i = 0
  while i < NEVENTS: 

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
    header[0] = len(FSPDGs)   # daughters + mother + electron + proton
    header[4] = BeamPolarization
    header[9] = w  
    
    # array to store all the particles in the event    
    parts = []

    # the reaction
    for l,j in enumerate(FSPDGs):
      state = 1
      #if j == PDGmother:
        #state = 0
      parts.append( storeParticle( j, Event.GetDecay(l) , vertex, state ) )

    if applyCuts( parts ) == False:
      continue


    #ve = Event.GetDecay(0)
    #vp = Event.GetDecay(1)
    ve = TVector3( parts[-2][6],parts[-2][7],parts[-2][8])
    vp = TVector3( parts[-1][6],parts[-1][7],parts[-1][8])

    h2p.Fill( math.degrees(vp.Phi()), math.degrees(ve.Phi()) ) 
    h2t.Fill( math.degrees(vp.Theta()), math.degrees(ve.Theta()) ) 

    # write the event
    for h in header:
      fout.write(str(h))
      fout.write(" ")
    fout.write("\n")

    writeParticles( fout, parts )
    i += 1

  fout.close()

h2t.Draw("colz")
