#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <tuple>
#include <sstream>
#include <cmath>

using namespace std;
//t x y z mass p0 px py pz pdg ID charge

//p0 px py pz pdg ch
//pt phi y pdg ch

typedef vector<tuple<double, double, double, int, int>> particles; //pt phi y pdg ch

void particle_list()
{
    cout<<" 211=pi+; -211=pi- \n"
          " 321=K+; -321=K-\n"
          " 2212=p; -2212=-p; \n\n";

    cout<<" Calculating: \n"
          "             pT spectra of pions, kaons, and protons \n"
          "             rapidity distribution of charged particles \n"
          "             v_2 of pions, kaons, and protons \n\n";
}

void fill(int n_bins, const double min, const double max, double spectra[], double value)
{
    for(int i=0 ; i<n_bins; i++)
        if(value>((i*((max-min)/n_bins))+min) && value<=(((i+1)*((max-min)/n_bins))+min))
        {
            spectra[i]++;
            return;
        }
}

int main() {
    particle_list();
    particles event;
    ifstream is(R"(/home/josef/School/smash/build/data/1/particle_lists.oscar)");
    ofstream os(R"(/home/josef/School/smash/build/data/1/output.dat)");

    const int n_bins = 20;
    const double pt_min = 0.;//GeV
    const double pt_max = 2.;//GeV
    const double y_min = -4.;
    const double y_max = 4.;

    double pt_spectra_pi[n_bins]; //pT spectra of pions
    double pt_spectra_K[n_bins]; //pT spectra of kaons
    double pt_spectra_p[n_bins]; //pT spectra of protons

    double y_spectra[n_bins]; //rapidity distribution of charged hadrons

    double v2_spectra_pi[n_bins]; //v_2(pT) of pions
    double v2_spectra_K[n_bins]; //v_2(pT) of kaons
    double v2_spectra_p[n_bins]; //v_2(pT) of protons

    int n_pi[n_bins]; //number of pions in pT bin
    int n_K[n_bins]; //number of kaons in pT bin
    int n_p[n_bins]; //number of protons in pT bin

    double p0, px, py, pz, pt, phi, y, skip; //particle properties
    int pdg, ch; //particle properties
    double Q2x_pi, Q2y_pi, Q2x_K, Q2y_K, Q2x_p, Q2y_p, psi2_pi, psi2_K, psi2_p;  //flow vector Q, event plane angle psi

    int N_eve = 0; //number of events

    for(int i=0; i<n_bins; i++) //initialization of  arrays
    {
        pt_spectra_pi[i]=0;
        pt_spectra_K[i]=0;
        pt_spectra_p[i]=0;
        y_spectra[i]=0;
        v2_spectra_pi[i]=0;
        v2_spectra_K[i]=0;
        v2_spectra_p[i]=0;
        n_pi[i]=0;
        n_K[i]=0;
        n_p[i]=0;
    }
    string line; //reading intro of the file
    getline(is, line); //#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge
    getline(is, line); //# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e
    getline(is, line); //# SMASH-2.2.1
    getline(is, line); //# event ...
    getline(is, line); //DATA
    do {
        N_eve++;
        if(!(N_eve%100))
            cout<<"Event number " << N_eve << "...";
        while (line[0] != '#' && !is.eof()) //while not "# event i end 0 impact   6.000 scattering_projectile_target yes"
        {//t x y z mass p0 px py pz pdg ID charge
            stringstream ss(line);
            ss >> skip;/*t*/ ss >> skip;/*x*/ ss >> skip;/*y*/ ss >> skip;/*z*/ ss >> skip;/*mass*/
            ss >> p0;
            ss >> px;
            ss >> py;
            ss >> pz;
            ss >> pdg;
            ss >> skip;/*ID*/
            ss >> ch;

            pt = sqrt((px * px) + (py * py)); //transverse momentum pT
            phi = atan2(py, px); //azimuthal angle phi
            y = log((p0 + pz) / (p0 - pz)) / 2; //rapidity

            event.push_back(make_tuple(pt, phi, y, pdg, ch)); //add particle to the event
            //cout<<p0<<"   "<<px<<"   "<<py<<"   "<<pz<<"   "<<pdg<<endl;
            //cout<<pt<<"   "<<phi<<"   "<<y<<"   "<<pdg<<"   "<<ch<<endl;
            getline(is, line); //next particle
        }

        //good
        Q2x_pi=0; Q2y_pi=0; Q2x_K=0; Q2y_K=0; Q2x_p=0; Q2y_p=0;

        for (auto i: event)
        {
            if (get<4>(i) != 0) //charged
                fill(n_bins, y_min, y_max, y_spectra, get<2>(i));
            if (abs(get<3>(i)) == 211 && fabs(get<2>(i))<=0.5) //pions, |y|<=0.5
            {
                fill(n_bins, pt_min, pt_max, pt_spectra_pi, get<0>(i));
                Q2x_pi+=get<0>(i) * cos(2*get<1>(i));
                Q2y_pi+=get<0>(i) * sin(2*get<1>(i));
            }
            if (abs(get<3>(i)) == 321 && fabs(get<2>(i))<=0.5) //kaons, |y|<=0.5
            {
                fill(n_bins, pt_min, pt_max, pt_spectra_K, get<0>(i));
                Q2x_K+=get<0>(i) * cos(2*get<1>(i));
                Q2y_K+=get<0>(i) * sin(2*get<1>(i));
            }
            if (abs(get<3>(i)) == 2212 && fabs(get<2>(i))<=0.5) //protons, |y|<=0.5
            {
                fill(n_bins, pt_min, pt_max, pt_spectra_p, get<0>(i));
                Q2x_p+=get<0>(i) * cos(2*get<1>(i));
                Q2y_p+=get<0>(i) * sin(2*get<1>(i));
            }
        }

        psi2_pi = atan2(Q2y_pi, Q2x_pi)/2;
        psi2_K = atan2(Q2y_K, Q2x_K)/2;
        psi2_p = atan2(Q2y_p, Q2x_p)/2;

        for(int b=0; b<n_bins; b++)
        {
            for(auto i : event)
            {
                if (abs(get<3>(i)) == 211 && fabs(get<2>(i))<=0.5
                    && get<0>(i)>(b*((pt_max-pt_min)/n_bins)) && get<0>(i)<((b+1)*((pt_max-pt_min)/n_bins))) //pions, |y|<=0.2, pt=bin
                {
                    v2_spectra_pi[b]+= cos(2*(get<1>(i)-psi2_pi));
                    n_pi[b]++;
                }

                if (abs(get<3>(i)) == 321 && fabs(get<2>(i))<=0.5
                    && get<0>(i)>(b*((pt_max-pt_min)/n_bins)) && get<0>(i)<((b+1)*((pt_max-pt_min)/n_bins))) //kaons, |y|<=0.5, pt=bin
                {
                    v2_spectra_K[b]+= cos(2*(get<1>(i)-psi2_K));
                    n_K[b]++;
                }

                if (abs(get<3>(i)) == 2212 && fabs(get<2>(i))<=0.5
                    && get<0>(i)>(b*((pt_max-pt_min)/n_bins)) && get<0>(i)<((b+1)*((pt_max-pt_min)/n_bins))) //protons, |y|<=0.5, pt=bin
                {
                    v2_spectra_p[b]+= cos(2*(get<1>(i)-psi2_p));
                    n_p[b]++;
                }
            }
        }

        event.clear();
        getline(is, line); //# event ...
        getline(is, line); //DATA
        if(!(N_eve%100))
            cout<<"done"<<endl;
    } while (!is.eof());

    for (int i = 0; i < n_bins; i++) //normalization and scaling
    {
        y_spectra[i]/= N_eve * ((y_max-y_min)/n_bins) * 2 * M_PI;
        pt_spectra_pi[i] /= N_eve * 2 * M_PI * ((i*((pt_max-pt_min)/n_bins))+(((pt_max-pt_min)/n_bins)/2)) * ((pt_max-pt_min)/n_bins) * (0.4); //1/N_eve*2pi*pt*y
        pt_spectra_K[i] /= N_eve * 2 * M_PI * ((i*((pt_max-pt_min)/n_bins))+(((pt_max-pt_min)/n_bins)/2)) * ((pt_max-pt_min)/n_bins) * (0.4); //1/N_eve*2pi*pt*y
        pt_spectra_p[i] /= N_eve * 2 * M_PI * ((i*((pt_max-pt_min)/n_bins))+(((pt_max-pt_min)/n_bins)/2)) * ((pt_max-pt_min)/n_bins) * (0.4); //1/N_eve*2pi*pt*y
        if(n_pi[i]!=0)
            v2_spectra_pi[i] /= n_pi[i];
        if(n_K[i]!=0)
            v2_spectra_K[i] /= n_K[i];
        if(n_p[i]!=0)
            v2_spectra_p[i] /= n_p[i];
    }


    for (double i : y_spectra)
        os  << i << "  ";
    os  << endl;
    for (double i : pt_spectra_pi)
        os  << i << "  ";
    os  << endl;
    for (double i : pt_spectra_K)
        os  << i << "  ";
    os << endl;
    for (double i : pt_spectra_p)
        os << i << "  ";
    os << endl;
    for (double i : v2_spectra_pi)
        os << i << "  ";
    os << endl;
    for (double i : v2_spectra_K)
        os << i << "  ";
    os << endl;
    for (double i : v2_spectra_p)
        os << i << "  ";
    os << endl;

    return 0;
}