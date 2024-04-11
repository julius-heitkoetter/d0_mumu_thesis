#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TROOT.h>

// Utility function to load necessary library
void loadPhysicsLibrary() {
    if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");
}

// Function to generate and fill histograms for different decays
void generateAndFillHistogram(
        TGenPhaseSpace& event, 
        TRandom3& rand, 
        const double* masses, 
        const double* reconstructed_masses, 
        int numParticles, 
        double mass_mu, 
        double mass_nu, 
        double mass_pi, 
        double mass_dstar, 
        double mass_D0,
        double pt_resolution, 
        TH1F* dimuon_mass_histogram,
        TH1F* delta_m_histogram,
        int n_events
    ) {


        for (int n = 0; n < n_events; ++n) {

            TLorentzVector Dstar;
            Dstar.SetPtEtaPhiM(rand.Exp(1.0 / 0.15), 0, 0, mass_dstar);

            /////////////////////
            /// D* -> D0 pi /////
            /////////////////////

            Double_t masses_dstar_to_piD0[2] = {mass_D0, mass_pi} ;
            event.SetDecay(Dstar, 2, masses_dstar_to_piD0);
            double weight1 = event.Generate();

            TLorentzVector* D0 = event.GetDecay(0);
            TLorentzVector* pi = event.GetDecay(1);
            TLorentzVector piMeasured;
            piMeasured.SetPtEtaPhiM(
                rand.Gaus(pi->Pt(), pi->Pt() * pt_resolution), 
                pi->Eta(),
                pi->Phi(),
                mass_pi
            );


            /////////////////////
            ///// D0 -> XX //////
            /////////////////////

            event.SetDecay(*D0, numParticles, masses);
            double weight2 = event.Generate();
            
            // Reconstruct dstar from reconstructed pion and reconstructed d0 products
            TLorentzVector dstarMeasured(0., 0., 0., 0.);
            dstarMeasured += piMeasured;

            // Reconstruct D0 from its reconstructed products
            TLorentzVector pMeasuredTotal(0., 0., 0., 0.);

            double weight3 = 1;

            for (int i = 0; i < numParticles; ++i) {
                TLorentzVector* p = event.GetDecay(i);

                //if (p->M() == mass_mu ) {
                if (true) {

                    TLorentzVector pMeasured;
                    pMeasured.SetPtEtaPhiM(
                        rand.Gaus(p->Pt(), p->Pt() * pt_resolution),
                        p->Eta(),
                        p->Phi(),
                        mass_mu
                    );
                    pMeasuredTotal += pMeasured;
                    dstarMeasured += pMeasured;

                } else {

                    // particle -> mu nu
                    TGenPhaseSpace product_decay_event;
                    Double_t masses_p_to_MuNu[2] = {reconstructed_masses[i], mass_nu} ;
                    product_decay_event.SetDecay(*p, 2, masses_p_to_MuNu);
                    double weight4 = product_decay_event.Generate();
                    weight3 *= weight4;
                    TLorentzVector* p_mu = product_decay_event.GetDecay(0);


                    TLorentzVector pMeasured;
                    pMeasured.SetPtEtaPhiM(
                        rand.Gaus(p_mu->Pt(), p_mu->Pt() * pt_resolution),
                        p_mu->Eta(),
                        p_mu->Phi(),
                        reconstructed_masses[i]
                    );
                    pMeasuredTotal += pMeasured;
                    dstarMeasured += pMeasured;
                }
            }
 

            dimuon_mass_histogram->Fill(pMeasuredTotal.M(), weight1 * weight2 * weight3);
            delta_m_histogram->Fill(dstarMeasured.M() - pMeasuredTotal.M(), weight1 * weight2 * weight3);
        }
        dimuon_mass_histogram->Scale(1/dimuon_mass_histogram->Integral());
        delta_m_histogram->Scale(1/delta_m_histogram->Integral());
    }

void plotHistograms(TH1F* histograms[], int n_hists, TString outfilename, TString xaxis_name) {


    TString title = "#font[61]{CMS}";
    TString title2 = "#font[52]{Preliminary Simulation}";
    TLatex tex;
    tex.SetTextFont(42);
    tex.SetTextSize(0.05);
    tex.SetTextAlign(11);
    tex.SetNDC();
    TString title3= "MC     13.6 TeV";

    TCanvas *c1 = new TCanvas("c1","c1");
    c1->cd();
    c1->SetMargin(0.17,0.06,0.13,0.07);
    histograms[0]->SetTitle("");
    histograms[0]->SetMaximum(0.04);
    histograms[0]->GetXaxis()->SetTitle(xaxis_name);
    histograms[0]->GetYaxis()->SetTitle("Fraction");
    histograms[0]->SetStats(kFALSE);
    histograms[0]->GetXaxis()->SetTitleSize(0.05);

    TLegend *legend = new TLegend(0.6,0.70,0.85,0.8);

    histograms[0]->Draw("HIST");
    histograms[0]->SetLineWidth(2);
    histograms[0]->SetLineColor(1);
    legend->AddEntry(histograms[0], "TODO", "l");
    for(int i=1; i<n_hists; i++){
        histograms[i]->SetLineWidth(2);
        histograms[i]->SetLineColor(i+1);
        histograms[i]->Draw("HIST same");
        legend->AddEntry(histograms[i], "TODO", "l");
    }

    legend->Draw("same");

    tex.SetTextSize(0.05);
    tex.DrawLatex(0.21,0.88,title);
    tex.SetTextSize(0.04);
    tex.DrawLatex(0.29,0.88,title2);
    tex.DrawLatex(0.70,0.94,title3);

    c1->SaveAs(outfilename);

}

TH1F* get_histogram(TString file_name, TString tree_name, TString variable_name) {
    TFile *fin = new TFile(file_name.Data());
    TTree *info = (TTree *)fin->Get("info");
    TTree *tree = (TTree *)fin->Get(tree_name.Data());

    TString cuts = "!(abs(mc_match)==413 && d0_d1_pt>4. && d0_d2_pt>4. && chan==2 && HLT_DoubleMu4_3_LowMass && dm<0.15 && dm>0.14 && d0_m<1.94 && d0_m>1.81 && dstar_vtx_prob>0.1 && d0_sl3d>3 && d0_alpha<0.1)";
    TH1F* h = new TH1F("h_" + variable_name, "Mass of D0; Mass (GeV/c^2); Entries", 50, 1.81, 1.94);
    tree->Draw(variable_name + ">>h_" + variable_name, cuts, "goff");
    //tree->Draw("d0_m>>h_d0_m");

    h->Scale(1/h->Integral());

    return h;
}

void background_analysis() {
    loadPhysicsLibrary();

    double mass_mu = 0.1057;
    double mass_nu = 0.0;
    double mass_pi = 0.1396;
    double mass_K = 0.4937;
    double mass_D0 = 1.865;
    double mass_dstar = 2.010;

    // Initial setup
    TGenPhaseSpace event;
    TRandom3 rand;
    double pt_resolution = 0.013;
    int n_events = 1000000;

    // Decay modes and their respective histograms
    int n_decays = 4;
    const char* decayNames1[] = {"h0", "h1", "h2", "h3"};
    const char* decayNames2[] = {"h4", "h5", "h6", "h7"};
    double decayMasses[][3] = {{mass_mu, mass_mu}, {mass_pi, mass_pi}, {mass_K, mass_K},  {mass_K, mass_pi}};
    double reconstructed_masses[][3] = {{mass_mu, mass_mu},{mass_mu, mass_mu},{mass_mu, mass_mu}, {mass_mu, mass_mu}};
    int decayParticles[] = {2, 2, 2, 2};

    TH1F* dimuon_mass_histograms[n_decays];
    TH1F* delta_m_histograms[n_decays];

    for (int i = 0; i < n_decays; ++i) {
        dimuon_mass_histograms[i] = new TH1F(decayNames1[i], decayNames1[i], 1000, 1.45, 2);
        delta_m_histograms[i] = new TH1F(decayNames2[i], decayNames2[i], 1000, 0.14, 0.155);

        generateAndFillHistogram(
            event, 
            rand, 
            decayMasses[i], 
            reconstructed_masses[i],
            decayParticles[i], 
            mass_mu, 
            mass_nu,
            mass_pi,
            mass_dstar,
            mass_D0,
            pt_resolution, 
            dimuon_mass_histograms[i],
            delta_m_histograms[i],
            n_events
        );
    }

    plotHistograms(dimuon_mass_histograms, n_decays, "reconstructed_D0_mass.pdf", "Dimuon Mass (GeV)");
    plotHistograms(delta_m_histograms, n_decays, "reconstructed_delta_m.pdf", "Delta M (GeV)");

    TH1F* inclusive_background_histogram_d0m = get_histogram("small-inclusive.root", "dzmmMcBkg", "d0_m");
    TH1F* inclusive_background_histograms_d0m[1] = {inclusive_background_histogram_d0m};
    plotHistograms(inclusive_background_histograms_d0m, 1, "inclusive_background_dimuon_mass.pdf", "Dimuon Mass (GeV)");

    TH1F* inclusive_background_histogram_dm = get_histogram("small-inclusive.root", "dzmmMcBkg", "dm");
    TH1F* inclusive_background_histograms_dm[1] = {inclusive_background_histogram_dm};
    plotHistograms(inclusive_background_histograms_dm, 1, "inclusive_background_delta_m.pdf", "Delta M (GeV)");
}

