#include <bits/stdc++.h>
using namespace std;

#define IOS ios_base::sync_with_stdio(0); cin.tie(0); cout.tie(0);


double th_in,th_out,tc_in,tc_out;
double m_h,m_c;
double api_h,api_c;
double ft; // from graph in book
double k_h,k_c;
double den_h,den_c;
double cp_h,cp_c;
double q; //total heat exchange
double vis_h,vis_c;
double u_assumed = 227.13, u_calc; 
double A;
double lmtd;
double nc,nt; // no of tubes
double d_i, d_o, d_s;
double L;
double Nt;//no of tube passes
double re_tube,re_shell; // reynolds numbers
double pt; // pitch
double d_eq; //equivalent diameter
double hi,ho; // heat transfer coefficients 
double rdi,rdo; // fouling factors
double k_wall; 
double rel_error; // in percentage
int iteration = 1;
double p_f,p_r,total_p;

int std_tubes[]={16,32,45,56,76,112,132,166,208,252,288,326,398,460,518,574,644};
int ds_arr[]={8,10,12,13.25,15.25,17.25,19.25,21.25,23.25,25,27,29,31,33,35,37,39};

void calcQ(){
    q = m_h * cp_h * (th_in - th_out);
}

void calcM_c(){
    m_c = q/(cp_c * (tc_out - tc_in));
}

void calcA(){
    A = q/(ft * lmtd * u_assumed);
}

void calcnt(){
    nc = A/(3.1416 * d_o * L);
    int i = upper_bound(std_tubes,std_tubes+17,nc) - std_tubes;
    nt = std_tubes[i];
    d_s = ds_arr[i] * 0.0254;
}
void calcP(){
    double vt = (4 * m_c * Nt) / (den_c * nt * 3.1416 * pow(d_i,2));
    double f = 0.05321;
    p_f = (8*f*pow(Nt,2)*pow(m_c,2)*L) / (pow(3.1416,2)*pow(d_i,5)*pow(nt,2)*den_c);
    p_r = (2*Nt*pow(vt,2)*1000) / (den_c * 9.8);
    total_p = p_f + p_r;
}

void solve(){
    cout<<endl;
    cout<<endl;
    cout<<"NOTE : all the values are in SI units\n"<<endl;
    for(int i =0;i<25;i++){
        cout<<"-";
    }
    cout<<"Iteration - "<<iteration;
    for(int i =0;i<25;i++){
        cout<<"-";
    }
    cout<<endl;
    iteration++;
    cout<<endl;

    
    calcA();
    cout<<"\nArea of heat transfer : "<<A;

    d_o = 0.0254; // in m
    d_i = 0.017018;
    L = 6;
    Nt=2;

    calcnt();
    cout<<"\nStandard no of tubes : "<<nt;
    cout<<"\nd_s : "<<d_s;

    // heat transfer coefficient (tube side)
    double re_tube = (Nt * m_c)/(3.1416 * d_i * vis_c/4000 * nt);
    cout<<"\n\nReynolds number of tube side flow : "<<re_tube;
    double nuselt_tube = 0.023 * pow(re_tube,0.8) * pow(cp_c * vis_c * 0.001/k_c , 0.4);

    hi = nuselt_tube * k_c/d_i;
    cout<<"\nHeat transfer coefficient of tube side flow : "<<hi<<endl;


    // heat transfer coefficient (shell side)
    pt = 1.25 * 0.0254;
    d_eq = 4 * (pow(pt,2) - 3.1416 * pow(d_o,2)/4) / (3.1416 * d_o);
    double c = pt - d_o;
    double b = d_s/2;
    double as = c * b * d_s /pt;
    double vs = m_h / (den_h * as);
    re_shell = den_h * vs * d_eq / (vis_h * 0.001);
    cout<<"\nReynolds number of the shell side flow : "<<re_shell;
    double jh = 0.806455 + 0.36921 * pow(re_shell,0.541044);

    ho = jh * k_h * pow(cp_h * vis_h * 0.001 / k_h , 0.333)/d_eq;
    cout<<"\nHeat transfer coefficient of shell side flow : "<<ho<<endl;


    rdi = 0.00075; //tube side
    rdo = 0.00018; //shell side
    k_wall = 54;
    u_calc = pow(1/ho + rdo + pow(d_o/d_i,2) * ((d_o-d_i)/(2*k_wall) + rdi + 1/hi),-1);
    cout<<"\nAssumed U : "<<u_assumed;
    cout<<"\nCalculated U : "<<u_calc<<endl;

    rel_error = (abs(u_assumed - u_calc)/ u_assumed )* 100;
    if((int)rel_error < 5){
        cout<<"\nRelative error in overall heat transfer coefficient : "<<rel_error;
        cout<<"\nRelative error is less than 5 percent, so this U is acceptable."<<endl;
        double overdesign = (ceil(nc)==27)?((nt-30)/30)*100:((nt-nc)/nc)*100;
        // double overdesign = ((nt-nc)/nc) *100;
        cout<<"\nOverdesign : "<<overdesign<<endl;
        if(overdesign<10){
            cout<<"Overdesign is also less than 10 percent so this overdesign is acceptable."<<endl;
            calcP();
            cout<<"\nPressure drop : "<<total_p<<endl;
            if((int)total_p<68947){ //10 psi = 68947 Pa
                cout<<"Pressure drop is also less than 10 psi so this design is completely acceptable and proceed to mechanical design."<<endl;
            }
        }else{
            cout<<"\nOverdesign is also more than 10 percent so please change desgin parameters.";
        }
        cout<<"\n\n";
        for(int i = 0 ;i<62;i++){
            cout<<'-';
        }
        cout<<endl;
    }else {
        cout<<"\nRelative error in overall heat transfer coefficient : "<<rel_error;
        cout<<"\nRelative error is more than 5 percent, so make calculated 'U' as assumed 'U' for next iteration."<<endl;
        u_assumed = u_calc;
        solve();
    }
}

int main(){
    IOS
    // #ifndef ONLINE_JUDGE
    // freopen("input.txt","r",stdin);
    // // freopen("output.txt","w",stdout);
    // #endif
    cout<<"Enter mass flow rate of heavynaptha : ";
    cin>>m_h;
    cout<<"Enter inflow temperature of heavynaptha : ";
    cin>>th_in;    
    cout<<"Enter outflow temperature of heavynaptha : ";
    cin>>th_out;    
    cout<<"Enter inflow temperature of water : ";
    cin>>tc_in;
    cout<<"Enter outflow temperature of water : ";
    cin>>tc_out;

    api_h = 47.3; // in degree celsius
    api_c = 10; 

    k_h = 0.15;
    k_c = 0.619;

    cp_h = 2036.9947;
    cp_c = 4180;

    den_h = 791.387;
    den_c = 994.10;

    vis_h = 0.733; //in cp
    vis_c = 0.8; //in cp    

    ft = 0.9575;

    lmtd = ((th_in - tc_out) - (th_out - tc_in))/log((th_in - tc_out)/(th_out - tc_in));
    cout<<"\n\nLMTD of this system : "<<lmtd;
    calcQ();
    cout<<"\nTotal heat exchanged : " <<q;
    calcM_c();
    cout<<"\nMass flow rate of water : "<<m_c;

    solve();
}
