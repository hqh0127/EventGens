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


NEVENTS= 1000

## FileName
fname = "elastic_E6.5_vx-0171vy01"
fname += ".txt"

def storeParticle(  pdg, vm, vertex, state, par=0, daug1=0, daug2=0 ):
  particle = [0]*14
  particle[0] = counter;
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
    outf.write("{:1d}".format(p[0]))
    for h in p[1:6]:
        outf.write("{:>7d}".format(h))
    for h in p[6:]:
        outf.write("{:> 12.6g}".format(h))
    outf.write("\n")

def applyCuts( particles ):
  test = True
  pdic = {}
  for p in particles:
    if p[2] in pdic:
      pdic[ p[2]+1000 ] = TVector3( p[6], p[7], p[8] )
    else:
      pdic[ p[2] ] = TVector3( p[6], p[7], p[8] )


  # check electron angle
  electron = pdic[ 11 ]
  if math.pi - electron.Theta()  < ThetaMin:
    return False
  if math.pi - electron.Theta() > ThetaMax:
    return False

  # check proton angle
  proton = pdic[ 2212 ]
  if math.pi - proton.Theta()  < PThetaMin:
    return False
  if math.pi - proton.Theta() > PThetaMax:
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
  fout.write("""SIMPLE Event FILE
============================================
I, ievent, nParticles
============================================
I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)
============================================\n""")
  parts = []
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
    header = [99]*3
    header[0] = 0
    header[1] = i
    header[2] = len(FSPDGs)   # daughters + mother + electron + proton

    # array to store all the particles in the event
    parts.clear()
    # the reaction

    parts.append( storeParticle(11, beam, vertex, 21, 0, 3, 4) )
    parts.append( storeParticle(2212, target, vertex, 21, 0, 5, 0) )
    scatlepton = storeParticle( 11, Event.GetDecay(0) , vertex, 1, 1 )
    scathardon = storeParticle( 2212, Event.GetDecay(1) , vertex, 1, 2 )
    boson = storeParticle(22, beam-Event.GetDecay(0), vertex, 21, 1)
    boson[10] = (beam-Event.GetDecay(0)).Mag()
    parts.extend([scatlepton, boson, scathardon])
    
    #for l,j in enumerate(FSPDGs):
    #  state = 1
    #  #if j == PDGmother:
    #    #state = 0
    #  parts.append( storeParticle( j, Event.GetDecay(l) , vertex, state, l+1 ) )

    if applyCuts( parts ) == False:
      continue

    #ve = Event.GetDecay(0)
    #vp = Event.GetDecay(1)
    ve = TVector3( parts[-2][6],parts[-2][7],parts[-2][8])
    vp = TVector3( parts[-1][6],parts[-1][7],parts[-1][8])

    h2p.Fill( math.degrees(vp.Phi()), math.degrees(ve.Phi()) )
    h2t.Fill( math.degrees(vp.Theta()), math.degrees(ve.Theta()) )

    # write the event
    fout.write("{:1d}".format(header[0]))
    fout.write("{:>10d}".format(header[1]))
    fout.write("{:>10d}".format(header[2]))
    fout.write("\n============================================\n")

    writeParticles( fout, parts )
    i += 1
    fout.write("=============== Event finished ===============\n")

  fout.close()

h2t.Draw("colz")
