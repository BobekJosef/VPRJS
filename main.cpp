#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <tuple>
#include <sstream>
#include <cmath>
#include <iomanip>


using namespace std;
//t x y z mass p0 px py pz pdg ID charge

//p0 px py pz pdg ch
//pt phi y pdg ch

typedef vector<tuple<double, double, double, int, int>> particles; //pt phi y pdg ch

void particle_list()
{
    cout<<" 211=pi+; -211=pi- \n"
          " 321=K+; -321=K-\n"
          " 2212=p; -2212=-p; \n";
}

void fill(int n_bins, const double min, const double max, double spectra[], double value)
{
    double dx=(max-min)/n_bins;
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
    double pt_spectra_pi[n_bins];// 1/(N_eve*2pi*pt) (dN/dpT)
    double pt_spectra_K[n_bins];
    double pt_spectra_p[n_bins];
    double y_spectra[n_bins];// 1/(N_eve*2pi) (dN/dy)
    double v2_spectra_pi[n_bins];
    double v2_spectra_K[n_bins];
    double v2_spectra_p[n_bins];
    double v2_eve_pi[n_bins];
    double v2_eve_K[n_bins];
    double v2_eve_p[n_bins];

    for(int i=0; i<n_bins; i++)
    {
        pt_spectra_pi[i]=0;
        pt_spectra_K[i]=0;
        pt_spectra_p[i]=0;
        y_spectra[i]=0;
        v2_spectra_pi[i]=0;
        v2_spectra_K[i]=0;
        v2_spectra_p[i]=0;
        v2_eve_pi[i]=0;
        v2_eve_K[i]=0;
        v2_eve_p[i]=0;
    }

    int N_eve = 0;

    string line;
    double p0, px, py, pz, pt, phi, y, skip, Q2x_pi, Q2y_pi, Q2x_K, Q2y_K, Q2x_p, Q2y_p, psi2_pi, psi2_K, psi2_p, v2_tmp_pi, v2_tmp_K, v2_tmp_p;
    int pdg, ch, n_pi, n_K, n_p;

    getline(is, line);//#!OSCAR2013 particle_lists t x y z mass p0 px py pz pdg ID charge
    getline(is, line);//# Units: fm fm fm fm GeV GeV GeV GeV GeV none none e
    getline(is, line);//# SMASH-2.2.1
    getline(is, line);//# event ...
    getline(is, line);//DATA
    do {
        N_eve++;
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

            pt = sqrt((px * px) + (py * py));
            phi = atan2(py, px);
            y = log((p0 + pz) / (p0 - pz)) / 2;

            event.push_back(make_tuple(pt, phi, y, pdg, ch));
            //cout<<p0<<"   "<<px<<"   "<<py<<"   "<<pz<<"   "<<pdg<<endl;
            //cout<<pt<<"   "<<phi<<"   "<<y<<"   "<<pdg<<"   "<<ch<<endl;
            getline(is, line);
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
            v2_tmp_pi=0; v2_tmp_K=0; v2_tmp_p=0;
            n_pi=0; n_K=0; n_p=0;

            for(auto i : event)
            {
                if (abs(get<3>(i)) == 211 && fabs(get<2>(i))<=0.5
                    && get<0>(i)>(b*((pt_max-pt_min)/n_bins)) && get<0>(i)<((b+1)*((pt_max-pt_min)/n_bins))) //pions, |y|<=0.2, pt=bin
                {
                    v2_tmp_pi+= cos(2*(get<1>(i)-psi2_pi));
                    n_pi++;
                }

                if (abs(get<3>(i)) == 321 && fabs(get<2>(i))<=0.5
                    && get<0>(i)>(b*((pt_max-pt_min)/n_bins)) && get<0>(i)<((b+1)*((pt_max-pt_min)/n_bins))) //kaons, |y|<=0.5, pt=bin
                {
                    v2_tmp_K+= cos(2*(get<1>(i)-psi2_K));
                    n_K++;
                }

                if (abs(get<3>(i)) == 2212 && fabs(get<2>(i))<=0.5
                    && get<0>(i)>(b*((pt_max-pt_min)/n_bins)) && get<0>(i)<((b+1)*((pt_max-pt_min)/n_bins))) //protons, |y|<=0.5, pt=bin
                {
                    v2_tmp_p+= cos(2*(get<1>(i)-psi2_p));
                    n_p++;
                }
            }

            //cout<<n_pi<<"  "<<n_K<<"  "<< n_p << endl;

            if (n_pi != 0) {
                v2_tmp_pi /= n_pi;
                v2_eve_pi[b]++;
            }
            if (n_K != 0) {
                v2_tmp_K /= n_K;
                v2_eve_K[b]++;
            }
            if (n_p != 0) {
                v2_tmp_p /= n_p;
                v2_eve_p[b]++;
            }

            v2_spectra_pi[b]+=v2_tmp_pi;
            v2_spectra_K[b]+=v2_tmp_K;
            v2_spectra_p[b]+=v2_tmp_p;
        }

        event.clear();
        getline(is, line);//# event ...
        getline(is, line);//DATA
        cout<<"done"<<endl;
    } while (!is.eof());

    for (int i = 0; i < n_bins; i++) //normalization and scaling
    {
        y_spectra[i]/= N_eve * ((y_max-y_min)/n_bins) * 2 * M_PI;
        pt_spectra_pi[i] /= N_eve * 2 * M_PI * ((i*((pt_max-pt_min)/n_bins))+(((pt_max-pt_min)/n_bins)/2)) * ((pt_max-pt_min)/n_bins) * (0.4); //1/N_eve*2pi*pt*y
        pt_spectra_K[i] /= N_eve * 2 * M_PI * ((i*((pt_max-pt_min)/n_bins))+(((pt_max-pt_min)/n_bins)/2)) * ((pt_max-pt_min)/n_bins) * (0.4); //1/N_eve*2pi*pt*y
        pt_spectra_p[i] /= N_eve * 2 * M_PI * ((i*((pt_max-pt_min)/n_bins))+(((pt_max-pt_min)/n_bins)/2)) * ((pt_max-pt_min)/n_bins) * (0.4); //1/N_eve*2pi*pt*y
        //v2_spectra_pi[i] /= v2_eve_pi[i];
        //v2_spectra_K[i] /= v2_eve_K[i];
        //v2_spectra_p[i] /= v2_eve_p[i];
        v2_spectra_pi[i] /= N_eve;
        v2_spectra_K[i] /= N_eve;
        v2_spectra_p[i] /= N_eve;
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