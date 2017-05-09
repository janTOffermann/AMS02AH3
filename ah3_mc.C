//
//  ah3_mc.C
//
//
//  Created by Jan Offermann on 04/30/17.
//
//

#include <string>
#include <cctype>
#include <vector>
#include <algorithm>
#include <random>
#include <tuple>

#include <math.h>
#include <stdio.h>

#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TH1D.h"
#include "TF1.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TRandom3.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

/*
 * Simulating events where W (at rest) decays to
 * 3 anti-nucleons (2 anti-neutrons, anti-proton to
 * make anti-tritium, which quickly decays to
 * anti-helium-3). There may be other products
 * (other quarks + gluons etc.) but I will currently
 * assume the anti-nucleons to dominate.
 */

const Double_t mW = 80.385; // W boson mass, GeV
const Double_t mpbar = 0.9382720813; // antiproton mass, GeV
const Double_t mnbar = 0.939565560; // antineutron mass, GeV
const Double_t pc = 0.357; // coalescence momentum, GeV
const Double_t Np = 1.046 / 2.; // expectation value of # of pbar or nbar (factor of 1/2 introduced after rereading paper)

/*
 * The following information for the pbar/nbar momentum distribution is taken from
 * Abe et al, Physical Review D69 072003. This can be found at: https://journals.aps.org/prd/pdf/10.1103/PhysRevD.69.072003
 */
const Double_t prod_weights [36] = {14.51, 17.32, 13.75, 11.12, 10.75, 9.048, 7.669, 7.410, 6.587, 5.788, 5.344, 4.987, 4.278, 4.117, 3.633, 3.036, 2.568, 2.165, 1.931, 1.603, 0.871, 0.912, 0.775, 0.639, 0.511, 0.419, 0.358, 0.254, 0.173, 0.141, 0.095, 0.0688, 0.047, 0.0241, 0.0093, 0.0015};
const Double_t xp_low [36] = {0.014, 0.016, 0.022, 0.027, 0.033, 0.038, 0.044, 0.049, 0.055, 0.06, 0.066, 0.071, 0.077, 0.082, 0.088, 0.099, 0.11, 0.121, 0.143, 0.164, 0.186, 0.208, 0.23, 0.252, 0.274, 0.296, 0.318, 0.351, 0.384, 0.417, 0.45, 0.482, 0.526, 0.57, 0.658, 0.768};
const Double_t xp_range [36] = {0.002, 0.006, 0.005, 0.006, 0.005, 0.006, 0.005, 0.006, 0.005, 0.006, 0.005, 0.006, 0.005, 0.006, 0.011, 0.011, 0.011, 0.021, 0.021, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.033, 0.033, 0.033, 0.033, 0.032, 0.044, 0.044, 0.088, 0.11, 0.232};

std::mt19937 rng(1);    // random-number engine used (Mersenne-Twister in this case)
std::uniform_real_distribution<Double_t> uni(0.,1.); // guaranteed unbiased

Double_t weight_sum = 0.;
Double_t weight_distrib [36] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,};

/*
 * A helper function, for printing numbers
 * in scientific notation. It's a little bit
 * "hacky" but will get the job done.
 * (previously tested in code for PHY 312)
 */
TString ToScientific(Double_t val, Int_t decimals) {
    Double_t mag = val;
    Int_t tens = 0;
    
    TString ans = TString("");
    
    if (val == 0.) return TString("0");
    
    if(val < 1.){
        
        while (mag  < 1.) {
            mag = mag * 10.;
            tens --;
        }
    }
    
    else{
        
        while (mag  >= 10.) {
            mag = mag / 10.;
            tens ++;
        }
    }
    
    ans.Append(std::to_string(mag));
    if (decimals != 0) ans = TString(ans(0,ans.Index(".")+decimals+1));
    else ans = TString(ans(0,ans.Index(".")));
    ans.Append("\\times10^{");
    ans.Append(std::to_string(tens));
    ans.Append("}");
    return ans;
}

UInt_t Factorial(UInt_t input){
    
    if(input == 0) return 1;
    
    UInt_t ans = input;
    
    UInt_t temp = input - 1;
    
    while(temp > 0) {
        
        ans = ans * temp;
        temp--;
    }
    
    return ans;
}

/*
 * This helper function calculates some constant
 * values for distributing momenta. This could
 * technically be done outside of the code,
 * but it is convenient to call it once as a function.
 */
void Setup() {
    
    weight_sum = 0.;
    for (UInt_t i = 0; i < 36; i++){
        
        weight_sum += prod_weights[i];
        
        weight_distrib[i] = 0.;
        for (UInt_t j = 0; j < i+1; j++) weight_distrib[i] += prod_weights[j];
    }
}

std::pair<Double_t, UInt_t> XPLowerBound(Double_t input) {
    
    std::pair<Double_t, UInt_t> xp_index_pair;
    
    Double_t xp_lower_bound = 0.;
    UInt_t index = 0.;
    
    for (UInt_t i = 0; i < 36; i++) {
        
        if(input < weight_distrib[i]){
            
            xp_lower_bound = xp_low[i];
            index = i;
            break;
        }
    }
    
    xp_index_pair.first = xp_lower_bound;
    xp_index_pair.second = index;
    return xp_index_pair;
}

/*
 * We have a spectrum of x_p = 2p/(E_cm), will need to
 * convert to p from an x_p value. This simple function
 * will do that, to be used later on. Note that E_cm is
 * simply the mass of the W, as it is produced at rest
 * in the experiment from which we're drawing our data.
 */
Double_t XPToP(Double_t xp) {
    
    // E_cm = mW
    return mW * xp / 2.;
}

/*
 * Simulate a single event that produces N anti-nucleons,
 * check for coalescence. We will treat the decay of the W
 * as a series of 2-body interactions (i.e. nucleons are
 * created one-by-one).
 */
Double_t runOne(UInt_t N){
    
    Double_t ans = 0.;
    
    Double_t E; // energy, GeV
    Double_t p; // momentum magnitude, GeV
    Double_t px;
    Double_t py;
    Double_t pz;
    Double_t m; // mass of nucleon, GeV
    Double_t theta;
    Double_t phi;
    std::pair<Double_t, UInt_t> xp_index_pair; // used for finding momentum from distribution

    // speed of lab frame wrt W (starts at rest), used for Lorentz boosts between steps
    Double_t vx = 0.;
    Double_t vy = 0.;
    Double_t vz = 0.;

    TLorentzVector W; // lorentzvector for W boson
    TLorentzVector* nucleon; // lorentzvector for one nucleon
    std::vector<TLorentzVector> nucleon_vecs;
    std::vector<std::tuple<UInt_t, UInt_t, Double_t>> rel_p_vec;
    
    W.SetPxPyPzE(0.,0.,0.,mW);
    
    for (UInt_t i = 0; i < N; i++) {
        
        nucleon = new TLorentzVector();
        m = mpbar; // differences between proton and neutron are complicated, and treating both will introduce unnecessary complexity for N > 3
        
        // pick momentum magnitude
        xp_index_pair = XPLowerBound(uni(rng) * weight_sum);
        p = XPToP(xp_index_pair.first + uni(rng) * xp_range[xp_index_pair.second]);
        
        // pick angles theta and phi
        theta = acos(2. * uni(rng) - 1.);
        phi = 2. * M_PI * uni(rng);
        E = sqrt(m * m + p * p);
        
        if (E > W.E()){
            
            N = i;

            E = W.E();
            p = sqrt(E * E - m * m);
        }
        
        px = p * cos(phi) * sin(theta);
        py = p * sin(phi) * sin(theta);
        pz = p * cos(theta);
        
        nucleon->SetPxPyPzE(px, py, pz, E);
        
        vx = - W.Px() / W.M();
        vy = - W.Py() / W.M();
        vz = - W.Pz() / W.M();
        
        TLorentzRotation boost(vx, vy, vz);
        
        nucleon->Transform(boost);
        
        // Conserve 4-momentum.
        W = W - *nucleon;
        
        nucleon_vecs.push_back(*nucleon);
        
        
        /*
         * For some reason the boost occasionally results in nan, despite v < c.
         * For now, toss these results and redo. It's important to determine
         * if this the symptom of some serious issue, which could be throwing
         * off our results.
         */
        if (std::isnan(nucleon->Px() + nucleon->Py() + nucleon->Pz())){
            
            delete nucleon;
            return -1.;
        }
        
        delete nucleon;
    }
    
    // With the 3 momenta set, now check the relative 3-momenta. First, create a vector of all relative momenta
    for (UInt_t i = 0; i < nucleon_vecs.size(); i++) {
        
        for (UInt_t j = i + 1; j < nucleon_vecs.size(); j++) {
            
            Double_t rel_p = sqrt((nucleon_vecs[i].Px() - nucleon_vecs[j].Px()) * (nucleon_vecs[i].Px() - nucleon_vecs[j].Px()) + (nucleon_vecs[i].Py() - nucleon_vecs[j].Py()) * (nucleon_vecs[i].Py() - nucleon_vecs[j].Py()) + (nucleon_vecs[i].Pz() - nucleon_vecs[j].Pz()) * (nucleon_vecs[i].Pz() - nucleon_vecs[j].Pz()));
            
            rel_p_vec.push_back(std::make_tuple(i, j, rel_p));
//            std::cout << i << "\t" << j << "\t" << rel_p << std::endl;
        }
    }
//    nucleon_vecs.clear();
    
    // Now check for coalescence of any 3.
    for (UInt_t i = 0; i < rel_p_vec.size(); i++) {
        
        if(std::get<2>(rel_p_vec.at(i)) > pc) continue;
        
        // We now have rel_ab < pc. Seek a rel_ac.
        
        UInt_t a = std::get<0>(rel_p_vec.at(i));
        UInt_t b = std::get<1>(rel_p_vec.at(i));
        
        for (UInt_t j = 0; j < rel_p_vec.size(); j++) {
            
            if(std::get<0>(rel_p_vec.at(i)) > a) break;
            if(std::get<0>(rel_p_vec.at(i)) != a) continue;
            
            if(std::get<2>(rel_p_vec.at(j)) > pc) continue;
            
            // We now have rel_ac < pc. Now check rel_bc.
            UInt_t c = std::get<1>(rel_p_vec.at(j));
            
            for (UInt_t k = 0; k < rel_p_vec.size(); k++) {
                
                if(std::get<0>(rel_p_vec.at(k)) != b) continue;
                if(std::get<1>(rel_p_vec.at(k)) != c) continue;
                
                if(std::get<2>(rel_p_vec.at(k)) > pc) break;
                
                return 1;
            }
        }
    }

    return 0;
}

Double_t sim(UInt_t N_nucleons = 3, UInt_t Nevents = 1000){
    
    Double_t ans = 0.;
    
    Double_t events_run = 0.;
    Double_t events_success = 0.;
    Double_t events_total = (Double_t)Nevents;
    Double_t single_result = 0.;
    
    for (UInt_t i = 0; i < Nevents; i++) {
        
//        if (Nevents > 10000000 && i % 10000000 == 0) std::cout << Nevents << "\t" << i << std::endl;
        single_result = runOne(N_nucleons);
        if(single_result < 0.) {
            i--;
            continue;
        }
        events_success += single_result;
        events_run += 1.;
    }
    
    ans = events_success / events_run;
    return ans;
}

void Converge(const UInt_t n = 7, const UInt_t N_nucleons = 3){
    
    TString title = TString("\\text{Monte Carlo }\\epsilon \\text{ Computation, N}_{n} = ");
    TString canvas_title = TString("canvas_N");
    canvas_title.Append(std::to_string(N_nucleons));
    title.Append(std::to_string(N_nucleons));
    title.Append(";\\text{N}_{\\text{events}};\\epsilon \\equiv \\epsilon_{P}\\,\\,\\epsilon_{MC}");
    TString info = TString(std::to_string(n));
    info.Prepend("\\epsilon\\;(\\text{N}_{\\text{events}} \\text{ = }10^{");
    info.Append("},\\text{ N}_{n}\\text{ = }");
    info.Append(std::to_string(N_nucleons));
    info.Append(") = ");

    TCanvas* canvas;
    TGraph* plot;
    TPaveText* pave;
    
    Setup();
    
    Double_t Nevents_array [n];
    Double_t results_array [n];
    Double_t epsilon_poiss = exp(- Np)*pow(Np, N_nucleons) / Factorial(N_nucleons);
    
    for (UInt_t i = 0; i < n; i++) {
        Nevents_array[i] = 1.;
        for (UInt_t j = 0; j < i + 1; j++) {
            Nevents_array[i] = Nevents_array[i] * 10.;
        }
    }
    
    std::cout << "-- N = " << N_nucleons << " --" << std::endl;
    for (UInt_t i = 0; i < n; i++){
        
        results_array[i] = epsilon_poiss * sim(N_nucleons, (UInt_t)Nevents_array[i]);
        std::cout << i+1 <<"/" << n << std::endl;
    }
    
    for (UInt_t i = 0; i < n; i++) {
        std::cout << Nevents_array[i] << " | " << results_array[i] << std::endl;
    }
    
    canvas = new TCanvas(canvas_title, canvas_title, 1600, 800);
    plot = new TGraph(n, Nevents_array, results_array);
    plot->SetTitle(title);
    
    info.Append(ToScientific(results_array[n-1], 3));
    pave = new TPaveText(.2,.8,.5,.9, "NDC");
    pave->AddText(info);
    pave->SetTextSize(0.028);
    pave->SetFillColor(kWhite);
    pave->SetTextAlign(11);
//    pave->SetBorderSize(1);
    
    canvas->SetLogx();
    
    plot->Draw();
    pave->Draw();
}

void ConvergeMulti(const UInt_t n = 7, const UInt_t range = 3, const UInt_t start = 3){
    
    for (UInt_t i = 0; i < range; i++) {
        Converge(n, start + i);
    }
}








