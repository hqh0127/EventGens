// read.cxx
//
// Created by TB on 6/13/11.
// Copyright 2011 BNL. All rights reserved.
//
// Example of how to read a file produced by BuildTree for a simple analysis.
// To run, in ROOT do:
// root [0] .L /path/to/read.cxx
// root [1] read("myInputFile.root", 10000 )
#include "eicsmear/functions.h"

void read(TString inFileNames) {

   // If the analysis solely uses TTree::Draw statements,
   // you don't need to load
   // the shared library. You will receive warnings such as
   // Warning in <TClass::TClass>: no dictionary for class Particle
   // is available
   // but these can be ignored. However if you wish to work with the event
   // objects themselves, the shared library must be loaded:
   // Load the shared library, if not done automaticlly:
   //   gSystem->Load("/path/to/libeicsmear.so" );
   gSystem->Load("libeicsmear");
   std::cout << "Using eic-smear version: " << erhic::EicSmearVersionString << std::endl;
   gSystem->Load("libeicsmeardetectors");
   auto libeicsmearpath=gSystem->DynamicPathName("libeicsmear");
   auto libdetectorpath=gSystem->DynamicPathName("libeicsmeardetectors");
   std::cout << "Using these eic-smear libraries : "<< std::endl
	     << libeicsmearpath << std::endl
	     << libdetectorpath << std::endl;


   // The TTrees are named EICTree.
   // Create a TChain for trees with this name.
   TFile* f = new TFile(inFileNames.Data());
   TTree* tree = (TTree*)f->Get("EICTree");

   // Add the file(s) we want to analyse to the chain.
   // We could add multiple files if we wanted.
   //tree->Add(inFileNames); // Wild cards are allowed e.g. tree.Add("*.root" );
// tree.Add(/path/to/otherFileNames ); // etc...

   // Create an object to store the current event from the tree.
   // This is how we access the values in the tree.
   // If you want to use generator-specific values, then
   // the event type should match the type in the TTree. Valid types are
   // EventPythia, EventPepsi, EventRapgap, EventDjangoh, EventMilou.
   // If you only need shared quantities like x, Q2 and the particle list
   // you can use EventBase and the macro will be general for any Monte Carlo.
   erhic::EventSimple* event(NULL);// = new EventPythia;
// EventBase* event(NULL);

   // Now associate the contents of the branch with the buffer.
   // The events are stored in a branch named event:
   tree->SetBranchAddress("event", &event ); // Note &event, not event.

   // Now we can do some analysis...

   // We record the largest particle pT we find here:
   double highestPt(-1. );
   double highestQ2(-1. );

   // Histograms for our analysis.
   TH2D Q2VsX("Q2VsX",
              "Q^{2} vs. Bjorken x;log_{10}(x);log_{10}(Q^{2})",
               100, -5., 0., 100, -.8, 5. );
   TH2D* theta_ep = new TH2D("theta_ep", ";cos(#theta_{e});cos(#theta_{p})", 100, -1, 1, 100, -1, 1);
   TH1D* hpthp = new TH1D("hpthp",";cos(#theta_{p});",100,-1,1);
   TH1D* hpthe = new TH1D("hpthe",";cos(#theta_{e});",100,-1,1);
   Q2VsX.SetStats(kFALSE);
   theta_ep->SetStats(kFALSE);

   int nEvents = tree->GetEntries();

   std::cout << nEvents << std::endl;
   // Loop over events:

   for(int i(0); i < nEvents; ++i ) {

      // Read the next entry from the tree.
      tree->GetEntry(i);

      // Fill the Q2 vs. x histogram:
      Q2VsX.Fill(TMath::Log10(event->GetX()),
                 TMath::Log10(event->GetQ2()));

      if (event->GetQ2()<0.1) continue;
      // The event contains a vector (array) of particles.
      int nParticles = event->GetNTracks();

      // We now know the number of particles in the event, so loop over
      // the particles:
      const Particle* proton, *electron;
      proton = nullptr;
      electron = nullptr;
      for(int j(0); j < nParticles; ++j ) {
         const Particle* particle = event->GetTrack(j);
         // Let's just select charged pions for this example:
         int pdg = particle->GetPdgCode();
         int KS = particle->GetStatus();
         if (pdg==11 && KS==1) electron = particle; 
         if (pdg==2212 && KS==1) proton = particle;

         // Update the highest pT:
         if(particle->GetPt() > highestPt ) {
            highestPt = particle->GetPt();
         } // if
      } // for
      if (proton && electron){
         theta_ep->Fill(TMath::Cos(electron->GetTheta()), TMath::Cos(proton->GetTheta()));
         hpthe->Fill(TMath::Cos(electron->GetTheta()));
         hpthp->Fill(TMath::Cos(proton->GetTheta()));
         if (event->GetQ2() > highestQ2 ) {
            highestQ2 = event->GetQ2();
         }
      }
   } // for

   std::cout << "The highest pT was " << highestPt << " GeV/c" << std::endl;
   std::cout << "The highest Q2 was " << highestQ2 << " GeV^2/c" << std::endl;

   TCanvas canvas;
   Q2VsX.Draw("colz" );
   canvas.Print("Q2VsX.png" );
   theta_ep->Draw("colz");
   canvas.Print("Theta_ep.png");
   hpthe->Draw();
   canvas.Print("hpthe.png");
   hpthp->Draw();
   canvas.Print("hpthp.png");
}
