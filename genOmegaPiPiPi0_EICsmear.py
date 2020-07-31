#!/bin/env python
#
#  SIMPLE GENERATOR FORMAT
# https://eic.github.io/software/pythia6.html
"""
SIMPLE Event FILE
============================================
I, ievent, nParticles
============================================
I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)
============================================
0         0         1
============================================
I: line index, runs from 1 to nrTracks
K(I,1):	status code KS (1: stable particles 11: particles which decay 55; radiative photon)
K(I,2):	particle KF code (211: pion, 2112:n, ….)
K(I,3):	line number of parent particle
K(I,4):	normally the line number of the first daughter; it is 0 for an undecayed particle or unfragmented parton
K(I,5):	normally the line number of the last daughter; it is 0 for an undecayed particle or unfragmented parton.
P(I,1):	px of particle
P(I,2):	py of particle
P(I,3):	pz of particle
P(I,4):	Energy of particle
P(I,5):	mass of particle
V(I,1):	x vertex information
V(I,2):	y vertex information
V(I,3):	z vertex information
=============== Event finished ===============
"""
# all the not useful fields will be filled with 99

from ROOT import TGenPhaseSpace, TLorentzVector, TVector3, TDatabasePDG, TMath, TH2F

import math
import numpy as np
from random import random as rndm
from array import array

dbpdg = TDatabasePDG()

# Exclusive collision Settings
# ============================

BeamPolarization = -1 # not useful

BeamEnergy = 6.5
TargetMass = dbpdg.GetParticle( 2212 ).Mass()
LeptonMass = dbpdg.GetParticle( 11 ).Mass()


beam   = TLorentzVector( 0, 0, -math.sqrt(BeamEnergy**2-LeptonMass**2), BeamEnergy )
target = TLorentzVector( 0, 0, 0, TargetMass )
CME    = beam + target

PDGmother = 223
FSPDGs   = [ 11 , 2212, PDGmother  ] # collision produces omega
FSMasses = [  dbpdg.GetParticle(i).Mass() for i in FSPDGs ]

#  daughters
PDGs = [ 211, -211, 111 ]
DaughterMass = [ dbpdg.GetParticle(i).Mass() for i in PDGs ]

# decay pi0
gPDGs = [ 22, 22 ]
gMasses = [ dbpdg.GetParticle(i).Mass() for i in gPDGs ]

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


NEVENTS= 1000

## FileName
fname = "phasespace_"
fname += dbpdg.GetParticle(FSPDGs[-1]).GetName()
fname += "_E{:g}".format(BeamEnergy);
for p in PDGs:
  fname += "_" +  dbpdg.GetParticle(p).GetName()
fname += ".txt"

def storeParticle(  pdg, vm, vertex, state, par=0, daug1=0, daug2=0 ):
  particle = [0]*14
  particle[0] = 0;
  particle[1] = state
  particle[2] = pdg

  particle[3] = par
  particle[4] = daug1
  particle[5] = daug2

  particle[6] = vm.Px()
  particle[7] = vm.Py()
  particle[8] = vm.Pz()
  particle[9] = vm.E()
  particle[10] = dbpdg.GetParticle( pdg ).Mass()
  particle[11] = vertex[0]
  particle[12] = vertex[1]
  particle[13] = vertex[2]
  return particle


def writeParticles( outf, particles ):
  # write mother
  for i,p in enumerate( particles ):
    p[0] = i+1
    outf.write("{:<4d}".format(p[0]))
    for h in p[1:6]:
        outf.write("{:>7d}".format(h))
    for h in p[6:]:
        outf.write("{:> 12.6g}".format(h))
    outf.write("\n")

def applyCuts( particles ):
  #return True
  test = True
  pdic = {}
  for p in particles:
    if p[1] > 1:
      continue
    if p[2] in pdic:
      pdic[ p[2]+1000 ] = TVector3( p[6], p[7], p[8] )
    else:
      pdic[ p[2] ] = TVector3( p[6], p[7], p[8] )


  # check electron angle
  electron = pdic[ 11 ]
  if math.pi - electron.Theta()  < ThetaMin:
    #print("cut electron theta min: {:.2f}".format(math.degrees(math.pi - electron.Theta())))
    return False
  if math.pi - electron.Theta() > ThetaMax:
    #print("cut electron theta max: {:g}".format(math.degrees(math.pi - electron.Theta())))
    return False

  # check proton angle
  proton = pdic[ 2212 ]
  if math.pi - proton.Theta()  < PThetaMin:
    #print("cut proton theta min: {:g}".format(math.degrees(math.pi - proton.Theta())))
    return False
  if math.pi - proton.Theta() > PThetaMax:
    #print("cut proton theta max: {:g}".format(math.degrees(math.pi - proton.Theta())))
    return False

  #check pions
  pi = pdic[ 211 ]
  if math.pi - pi.Theta()  < ThetaMin:
    #print("cut pi+ theta min: {:g}".format(math.degrees(math.pi - pi.Theta())))
    return False
  if math.pi - pi.Theta() > ThetaMax:
    #print("cut pi+ theta max: {:g}".format(math.degrees(math.pi - pi.Theta())))
    return False

  pi = pdic[ -211 ]
  if math.pi - pi.Theta()  < ThetaMin:
    #print("cut pi- theta min: {:g}".format(math.degrees(math.pi - pi.Theta())))
    return False
  if math.pi - pi.Theta() > ThetaMax:
    #print("cut pi- theta max: {:g}".format(math.degrees(math.pi - pi.Theta())))
    return False

  # check photon
  pi = pdic[ 22 ]
  if math.pi - pi.Theta() > ThetaMax:
    #print("cut gamma1 theta max: {:g}".format(math.degrees(math.pi - pi.Theta())))
    return False
  pi = pdic[ 1022 ] # should detect the second gamma from pi0 decay
  if math.pi - pi.Theta() > ThetaMax:
    #print("cut gamma2 theta max: {:g}".format(math.degrees(math.pi - pi.Theta())))
    return False

  return test

# generate
with open(fname,"w") as fout:

  # generate event
  Event = TGenPhaseSpace()
  Event.SetDecay( CME, len(FSMasses),array('d',FSMasses))

  # loop on events
  i = 0
  fout.write("""SIMPLE Event FILE
============================================
I, ievent, nParticles
============================================
I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)
============================================\n""")
  parts = []
  while i < NEVENTS:

    # generate event
    w = Event.Generate()

    # generate vertex
    vertex = (
              VxMin + rndm()*(VxMax - VxMin),
              VyMin + rndm()*(VyMax - VyMin),
              VzMin + rndm()*(VzMax - VzMin)
             )

    # write header
    header = [99]*3
    header[0] = 0
    header[1] = i
    header[2] = 6   # 6 final state particles: e p pi+ pi- 2*gamma

    # array to store all the particles in the event
    parts.clear()
    # the reaction

    parts.append( storeParticle(11, beam, vertex, 21, 0, 3, 4) )
    parts.append( storeParticle(2212, target, vertex, 21, 0, 5, 5+len(FSPDGs)-2) )
    scatlepton = storeParticle( 11, Event.GetDecay(0) , vertex, 1, 1 )
    scathardon = storeParticle( 2212, Event.GetDecay(1) , vertex, 1, 2 )
    boson = storeParticle(22, beam-Event.GetDecay(0), vertex, 21, 1)
    boson[10] = (beam-Event.GetDecay(0)).Mag()
    parts.extend([scatlepton, boson, scathardon])
    
    index_hadron = 5+len(FSPDGs)-2
    parts.append( storeParticle(FSPDGs[-1], Event.GetDecay(2), vertex, 11, 2, index_hadron+1, index_hadron+len(PDGs)) )
    
    # prepare decay
    mother   = Event.GetDecay(2)
    Decay = TGenPhaseSpace()
    Decay.SetDecay(mother, len(DaughterMass),array('d',DaughterMass))

    Decay.Generate()
    for j,p in enumerate(PDGs):
      vp = Decay.GetDecay(j)
      state = 1
      if p == 111:
        state = 11
      if state == 1:
        parts.append( storeParticle( p, vp, vertex, state, index_hadron, 0, 0 ) )
      else:
        parts.append( storeParticle( p, vp, vertex, state, index_hadron, index_hadron+ j+2, index_hadron+j+1+len(gPDGs) ) )
      
    # decay the pi0
    Decay.SetDecay( Decay.GetDecay(2),len(gMasses),array('d',gMasses))
    Decay.Generate()
    for j,p in enumerate(gPDGs):
      vp = Decay.GetDecay(j)
      parts.append( storeParticle(  p, vp, vertex, 1, 9 ) )

    if applyCuts( parts ) == False:
      continue
    #print("success event {:>3d}".format(i))

    # write the event
    fout.write("{:1d}".format(header[0]))
    fout.write("{:>10d}".format(header[1]))
    fout.write("{:>10d}".format(header[2]))
    fout.write("\n============================================\n")

    writeParticles( fout, parts )
    i += 1
    fout.write("=============== Event finished ===============\n")

  fout.close()
