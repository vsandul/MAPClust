#include <string>
#include <fstream>

#include "../../hclust-cpp-master/fastcluster.cpp" 
#include "../src/particle_cluster.h"
#include "../src/particle_cluster.cpp"
#include "../src/dbscan.cpp"
#include "../src/dsu.cpp"
#include "../src/ccl.cpp"
#include "../src/cluster_quantities.h"

int Sign(Double_t x) { return x/fabs(x);}
double Pt(Double_t px, Double_t py) { return sqrt(px*px+py*py);} 
double Pabs(Double_t px, Double_t py, Double_t pz) { return sqrt(px*px+py*py+pz*pz);}
double P0(Double_t px, Double_t py, Double_t pz, Double_t m) { return sqrt(px*px+py*py+pz*pz + m*m);}
double Y_Eta(Double_t p, Double_t pl) { return 0.5*log((p+pl)/(p-pl));}
double Phi(Double_t px, Double_t py) { return atan2(py, px);} //azimutal angle phi

bool IsLightNucleus(int pdg);
bool IsHeavyNucleus(int pdg);

vector<SingleParticle> ReadOscar(ifstream& ifstr);

Int_t setcolor(Float_t x, Float_t y, Float_t z);
void plot();

 
void smash_clusterize(){
        
    ifstream f1("particle_lists.oscar");
    vector<SingleParticle> not_nucleons;
    vector<SingleParticle> nucleons;
    string s;
    while(getline(f1, s)){
        if (s[0] == '#') continue;
        stringstream ss(s);
        double t; ss >> t;
        double x; ss >> x;
        double y; ss >> y;
        double z; ss >> z;
        double mass; ss >> mass;
        double p0; ss >> p0;
        double px; ss >> px;
        double py; ss >> py;
        double pz; ss >> pz;
        int pdg; ss >> pdg;
        int ID; ss >> ID; 
        int charge; ss >> charge;
        SingleParticle particle = {ID, x,y,z,  px/p0, py/p0, pz/p0};
        if (pdg == 2212 || pdg == 2112)
            nucleons.push_back(particle);
        else
            not_nucleons.push_back(particle);
    }
    f1.close();
    
    BunchOfParticles bunch(nucleons);
    bunch.Clusterize(5.5, "HClust_Single", 0.35, "HClust_Complete");
    map<int,vector<SingleParticle>> map_of_labels = bunch.GetLabelsMap();
    vector<int> labels = bunch.GetLabelsVector();
    cout << "Num. of clusters: " << map_of_labels.size() << endl;    
    
    TCanvas *cvs = new TCanvas("cvs", "", 800, 600);
    gStyle->SetCanvasColor(0);
    //gStyle->SetMarkerStyle(20);
    //gStyle->SetMarkerSize(7);
    int colors[] = {2,3,4,5,6,7, kPink,kMagenta+3, kAzure-3,kOrange+7,kCyan-3, 
        kPink-6,kSpring-1, kViolet+3,kBlue+4, kGreen+3,kPink-7 }; 
    
    TNtuple *not_nucl = new TNtuple("n", "n", "x:y:z");
    for (Int_t i = 0; i < not_nucleons.size(); i++){
        not_nucl->Fill(not_nucleons.at(i).X, not_nucleons.at(i).Y, not_nucleons.at(i).Z);
    }
    not_nucl->SetMarkerColor(kGray);
    not_nucl->SetMarkerStyle(4);
    not_nucl->Draw("x:y:z");    
       
    for (const auto& [labl, parts]:map_of_labels){
        if (parts.size() == 1){
             TNtuple *prot_and_neut = new TNtuple("n", "n", "x:y:z");
             prot_and_neut->Fill(parts.at(0).X, parts.at(0).Y, parts.at(0).Z);
             prot_and_neut->SetMarkerColor(43);
             prot_and_neut->SetMarkerStyle(4);
             prot_and_neut->Draw("x:y:z","","same");
        } else {
            int color = colors[parts.size()];
            TNtuple *nucl = new TNtuple("n", "n", "x:y:z");
            for (const auto& p:parts){
                nucl->Fill(p.X, p.Y, p.Z);
            }
            nucl->SetMarkerStyle(8);
            nucl->SetMarkerColor(color);
            nucl->Draw("x:y:z","","same");            
        }
    }
}




bool IsLightNucleus(int pdg) {
    if (pdg < 1000000020)
        return false;
    string s_pdg = to_string(pdg);
    if (s_pdg[6] == '0' && s_pdg[7] == '0' && (s_pdg[8] == '2' || s_pdg[8] == '3' || s_pdg[8] == '4'))
        return true;
    else
        return false;
}

bool IsHeavyNucleus(int pdg) {
    if (pdg >= 1000000040 && !IsLightNucleus(pdg)) 
        return true;
    return false;
}
