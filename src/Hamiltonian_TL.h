#include "tensors.h"
#include "Parameters_TL.h"
#include "Neighbors_TL.h"
#include "MFParameters_TL.h"

extern "C" {
    void zheev_(char *, char *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, double *, int *);
    void dgesdd_ (char *, int *, int *, double *, int *, double *, double *, int *, double *, int*, double *, int *, int *, int *);
    void zgesdd_ (char *, int *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, std::complex<double> *, int*, 
                    std::complex<double> *, int *, double *, int *, int *);
}

#ifndef Hamiltonian_TL_class
#define Hamiltonian_TL_class
#define PI acos(-1.0)

class Hamiltonian_TL{

    public:
    //Constructor:
    Hamiltonian_TL(Parameters_TL &Parameters_TL__, Neighbors_TL &Neighbors_TL__, MFParam_TL &MFParam_TL__):
                     Parameters_TL_(Parameters_TL__), Neighbors_TL_(Neighbors_TL__), MFParam_TL_(MFParam_TL__){
        Initialize();
    }
    void Initialize();
    int makeIndex(int spin, int rx, int ry);
    void addHoppingTerms();
    void addInteractionTerms();
    void addPairingTerms();
    double getbondAngleNNPairing(int from_site, int to_site);
    void Diagonalizer();
    void performDiagonalization();
    double fermifunction(double en_);
    pair_doub getEnergies();
    void getNewOPs();
    void updateOrderParamsBroyden(int iter);
    void updateOrderParamsAnderson(int iter);
    double getOPError();
    void performSVD(Mat_2_doub &A, Mat_2_doub &VT, Mat_2_doub &U, Mat_1_doub &Sigma);
    void performComplexSVD(Mat_2_Complex_doub &A, Mat_2_Complex_doub &VT, Mat_2_Complex_doub &U, Mat_1_doub &Sigma);

    double calculateTotalParticles();
    double calculateLDOSatSite(int rx, int ry, double w_val);
    void calculateAvgdLDOS();
    void calculateQPI(double w_val);

    void calculateLDOSplusQPI();
    double Lorentzian(double val);
    void calculateQPIusingTmatrix(double w_val);
    
    
    Parameters_TL &Parameters_TL_;
    Neighbors_TL &Neighbors_TL_;
    MFParam_TL &MFParam_TL_;

    int lx_, ly_, nspin, nsites, Hsize, BdGsize;
    double onsiteE, B_, U0, U1;
    double mu_;

    Mat_1_doub evals_;
    Mat_2_Complex_doub C_mat;
    Mat_2_Complex_doub evecs_;
    Mat_2_Complex_doub uvecs_,vvecs_;
    
    Mat_1_doub hopping_params;
    Mat_1_Complex_doub OPs_, newOPs_;
    Mat_1_Complex_doub HF_OPs_, new_HF_OPs_;
    Mat_1_Complex_doub Pair_OPs_, new_Pair_OPs_;
    int OPs_size;
    int HF_OPs_size, Pair_OPs_size;
    int posN;
    
    double alpha, beta_damping,alpha_AA;

    //For Broyden Mixing:
    Mat_2_doub B_m, B_mp1;      //mixing (inverse Jacobian) matrices (real)
    Mat_1_doub del_X, del_F;    //difference in OPs and residuals
    Mat_1_doub F_m, F_mp1;

    //For Anderson Mixing:
    Mat_1_Complex_doub x_k, x_km1, del_x_k;
    Mat_1_Complex_doub f_k, f_km1, del_f_k;
    Mat_1_Complex_doub new_x_k, new_f_k, gamma_k, x_kp1;
    Mat_2_Complex_doub Xmat, Fmat, Gmat;
    int AM_depth;

    //for lorentzian:
    double eta;
    double wmin,wmax,dw;
    int wsize;
};

/*
double Hamiltonian_TL::getbondAngleNNPairing(int from_site, int to_site){
    auto [xi,yi] = Neighbors_TL_.getCoordinates(from_site);
    auto [xj,yj] = Neighbors_TL_.getCoordinates(to_site);

    int dx = xj - xi;
    int dy = yj - yi;
    if(Parameters_TL_.PBC_X==true){
        if(dx >  lx_/2){dx = dx - lx_;}
        if(dx < -lx_/2){dx = dx + lx_;}
    }
    if(Parameters_TL_.PBC_Y==true){
        if(dy >  ly_/2){dy = dy - ly_;}
        if(dy < -ly_/2){dy = dy + ly_;}
    }

    double dist_x = dx + 0.5*dy;
    double dist_y = 0.5*sqrt(3.0)*dy;
    double bond_angle = atan2(dist_y, dist_x);

    return bond_angle;
}
*/

void Hamiltonian_TL::Initialize(){
    hopping_params = {Parameters_TL_.t1_hop, Parameters_TL_.t2_hop, Parameters_TL_.t3_hop, Parameters_TL_.t4_hop, Parameters_TL_.t5_hop};
    //hopping_params = {Parameters_TL_.t1_hop, 0.0, 0.0, 0.0, 0.0};
    onsiteE = Parameters_TL_.onsiteE;
    B_ = Parameters_TL_.Bfield;
    U0 = Parameters_TL_.U0_;
    U1 = Parameters_TL_.U1_;
    mu_ = Parameters_TL_.mu_fixed;

    lx_ = Parameters_TL_.Lx_;
    ly_ = Parameters_TL_.Ly_;
    nspin = 2;
    nsites = Parameters_TL_.Total_Sites;
    Hsize = Parameters_TL_.Hsize;
    BdGsize = Parameters_TL_.BdG_Size;

    evals_.resize(BdGsize);
    evecs_.resize(BdGsize);     C_mat.resize(BdGsize);
    for(int i=0;i<BdGsize;i++){
        evecs_[i].resize(BdGsize);      C_mat[i].resize(BdGsize);
    }
    uvecs_.resize(Hsize);       vvecs_.resize(Hsize);
    for(int i=0;i<Hsize;i++){
        uvecs_[i].resize(BdGsize);    vvecs_[i].resize(BdGsize);
    }

    HF_OPs_size = MFParam_TL_.initial_HF_OPs_.size();
    HF_OPs_.resize(HF_OPs_size);
    new_HF_OPs_.resize(HF_OPs_size);
    cout<<"Number of HF OPs= "<<HF_OPs_size<<endl;

    Pair_OPs_size = MFParam_TL_.initial_Pair_OPs_.size();
    Pair_OPs_.resize(Pair_OPs_size);
    new_Pair_OPs_.resize(Pair_OPs_size);
    cout<<"Number of Pairing OPs= "<<Pair_OPs_size<<endl;
    
    OPs_size = HF_OPs_size + Pair_OPs_size;
    OPs_.resize(OPs_size);
    newOPs_.resize(OPs_size);
    alpha = Parameters_TL_.alpha_;
    //For Anderson mixing:
    beta_damping = Parameters_TL_.beta_damp;
    alpha_AA=1.0;
    AM_depth = 10;
    gamma_k.resize(AM_depth);
    x_k.resize(OPs_size);       f_k.resize(OPs_size);       x_kp1.resize(OPs_size);
    del_x_k.resize(OPs_size);   del_f_k.resize(OPs_size);   new_x_k.resize(OPs_size);
    Xmat.resize(OPs_size);      Fmat.resize(OPs_size);      new_f_k.resize(OPs_size);
    for(int x=0;x<OPs_size;x++){
        Xmat[x].resize(AM_depth);
        Fmat[x].resize(AM_depth);
    }
    //For QPI:
    wmax = Parameters_TL_.w_max;
    wmin = Parameters_TL_.w_min;
    dw   = Parameters_TL_.dw_;
    eta  = Parameters_TL_.eta_;
    wsize = (int) ( (wmax - wmin) / dw);
}

int Hamiltonian_TL::makeIndex(int spin, int rx, int ry){
    return spin*nsites + rx + lx_*ry;
}

void Hamiltonian_TL::addHoppingTerms(){
    //Hopping connections upto 6th neighbor:
    for(int spin=0;spin<2;spin++){
        for(int rx=0;rx<lx_;rx++){
            for(int ry=0;ry<ly_;ry++){
                int r = Neighbors_TL_.getIndex(rx,ry);
                int electron_index = spin*nsites + r;   //-> up-e: 0 ... lx*ly-1 & dn-e: lx*ly ... 2*lx*ly-1
                int hole_index = electron_index+ Hsize; //-> up-h: 2*lx*ly ... 3*lx*ly-1  & dn-h: 3*lx*ly ... 4*lx*ly-1

                //Hopping connections upto 6th neighbor:
                for(int order=0; order<Neighbors_TL_.dist_to_neigh_order.size(); order++){
                    Mat_1_int neighbors = Neighbors_TL_.findNeighbor(r,order);
                    double hopping = hopping_params[order];
                    
                    for(auto &neighbor : neighbors){
                        int electron_neigh = neighbor + nsites*spin;
                        int hole_neigh = Hsize + (neighbor + nsites*spin);
                        complex<double> hopval = (-1.0)*hopping * One_Complex;
                        //Upper-left block:
                        C_mat[electron_index][electron_neigh] = hopval;
                        C_mat[electron_neigh][electron_index] = conj(hopval);

                        //Lower-right block:
                        C_mat[hole_neigh][hole_index] = (-1.0)*hopval;
                        C_mat[hole_index][hole_neigh] = conj((-1.0)*hopval);
                    }
                }

                //Chemical potential term:
                C_mat[electron_index][electron_index] += (-1.0)*onsiteE*One_Complex;
                C_mat[hole_index][hole_index]         +=  (1.0)*onsiteE*One_Complex;
                //Magnetic field term:
                C_mat[electron_index][electron_index] +=        (pow(-1.0, 1.0 * spin))*B_*One_Complex;
                C_mat[hole_index][hole_index]         += (-1.0)*(pow(-1.0, 1.0 * spin))*B_*One_Complex;

                if(Parameters_TL_.singleImpurity==true){
                    double impVal = Parameters_TL_.impurityval;
                    if(rx==lx_/2 && ry==ly_/2){
                        if(Parameters_TL_.chargeDisorder==true){ //for charge impurity
                            C_mat[electron_index][electron_index] +=  (1.0)*impVal*One_Complex;
                            C_mat[hole_index][hole_index]         += (-1.0)*impVal*One_Complex;
                        }
                        else{ //for magentic impurity
                            C_mat[electron_index][electron_index] +=        (pow(-1.0, 1.0 * spin))*impVal*One_Complex;
                            C_mat[hole_index][hole_index]         += (-1.0)*(pow(-1.0, 1.0 * spin))*impVal*One_Complex;
                        }
                    }
                }
                else{

                }

            }
        }
    }
}

void Hamiltonian_TL::addInteractionTerms(){
    int ham_row, ham_col;

    for(int x=0;x<HF_OPs_size;x++){
        //OP[row_ind][col_ind] where row_ind<=col_ind
        int row_ind = MFParam_TL_.row_vector_HF[x];
        int col_ind = MFParam_TL_.col_vector_HF[x];
        assert(row_ind <= col_ind);
        assert(row_ind<Hsize && col_ind<Hsize);

        //index = spin*nsites + rx + lx_*ry
        int site_row = row_ind % nsites;
        int spin_row = row_ind / nsites;
        int site_col = col_ind % nsites;
        int spin_col = col_ind / nsites;

        //Hartree terms:
        if(row_ind==col_ind){
            //onsite Hartree term: U0 * <n_{i,s}>n_{i,bar(s)}
            ham_row = (1-spin_row)*nsites + site_row;
            ham_col = ham_row;
            C_mat[ham_row][ham_col] += U0 * One_Complex * HF_OPs_[x].real();
            C_mat[ham_row + Hsize][ham_col + Hsize] += -U0 * One_Complex * HF_OPs_[x].real();
            
            //NN Hartree term: U1 * {i in neigh}<n_j> n_i
            for(int site=0; site<nsites; site++){
                if(abs(Neighbors_TL_.interaction_pairs[site_row][site])>1e-5){
                    for(int spin=0;spin<nspin;spin++){
                        ham_row = spin*nsites + site;
                        ham_col = ham_row;
                        C_mat[ham_row][ham_col] += U1 * One_Complex * HF_OPs_[x].real();
                        C_mat[ham_row + Hsize][ham_col + Hsize] += -U1 * One_Complex * HF_OPs_[x].real();
                    }
                }
            }
        }
        else{//Fock terms:
            //onsite Fock term: 
            if(site_row == site_col){
                assert(spin_row==0);
                assert(spin_col==1);

                ham_row = site_row + spin_col*nsites; //i,dn
                ham_col = site_row + spin_row*nsites; //i,up
                C_mat[ham_row][ham_col] += (-1.0) * U0 * One_Complex * HF_OPs_[x];
                C_mat[ham_col][ham_row] += (-1.0) * U0 * One_Complex * conj(HF_OPs_[x]);

                C_mat[ham_col + Hsize][ham_row + Hsize] += (1.0) * U0 * One_Complex * HF_OPs_[x];
                C_mat[ham_row + Hsize][ham_col + Hsize] += (1.0) * U0 * One_Complex * conj(HF_OPs_[x]);
            }
            else{//offsite Fock terms:
                if(spin_row==spin_col){
                    assert(site_col>site_row);
                    if(abs(Neighbors_TL_.interaction_pairs[site_row][site_col])>1e-6){
                        //sr = site_col     sp_r = spin_row
                        //sc = site_row     sp_c = spin_row
                        ham_row = spin_row*nsites + site_col; //j,s'
                        ham_col = spin_row*nsites + site_row; //i,s
                        C_mat[ham_row][ham_col] += (-1.0) * U1 * One_Complex * HF_OPs_[x];
                        C_mat[ham_col][ham_row] += (-1.0) * U1 * One_Complex * conj(HF_OPs_[x]);

                        C_mat[ham_col + Hsize][ham_row + Hsize] += (1.0) * U1 * One_Complex * HF_OPs_[x];
                        C_mat[ham_row + Hsize][ham_col + Hsize] += (1.0) * U1 * One_Complex * conj(HF_OPs_[x]);
                    }
                }
                else{
                    assert(spin_row==0);
                    assert(spin_col==1);
                    if(abs(Neighbors_TL_.interaction_pairs[site_row][site_col])>1e-6){
                        //sr = site_col     sp_r = spin_col
                        //sc = site_row     sp_c = spin_row
                        ham_row = spin_col*nsites + site_col; //j,s'
                        ham_col = spin_row*nsites + site_row; //i,s
                        C_mat[ham_row][ham_col] += (-1.0) * U1 * One_Complex * HF_OPs_[x];
                        C_mat[ham_col][ham_row] += (-1.0) * U1 * One_Complex * conj(HF_OPs_[x]);
                        
                        C_mat[ham_col + Hsize][ham_row + Hsize] += (1.0) * U1 * One_Complex * HF_OPs_[x];
                        C_mat[ham_row + Hsize][ham_col + Hsize] += (1.0) * U1 * One_Complex * conj(HF_OPs_[x]);
                    }
                }
            }   

        }
        
    }
}

void Hamiltonian_TL::addPairingTerms(){
    double Vp_ = Parameters_TL_.pair_pot;

    if(Parameters_TL_.pairingType==0){//onsite s-wave pairing
        
        int e_i_up,e_i_dn;
        for(int x=0;x<Pair_OPs_size;x++){
            int prow_ind = MFParam_TL_.row_vector_Pair[x];
            int pcol_ind = MFParam_TL_.col_vector_Pair[x];
            assert(prow_ind <= pcol_ind);
            assert(prow_ind<Hsize && pcol_ind>=Hsize);
            
            int site_prow = prow_ind % nsites;
            int spin_prow = prow_ind / nsites;
            int site_pcol = (pcol_ind - Hsize) % nsites;
            int spin_pcol = (pcol_ind - Hsize) / nsites;
            assert(site_prow == site_pcol);
            assert(spin_prow==0 && spin_pcol==1);
            
            //Terms we need for onsite s-wave pairing:
            //a. -D_{i}(c^{dag}_{i,dn}c^{dag}_{i,up} - c^{dag}_{i,up}c^{dag}_{i,dn}) + h.c.
            //where D_{i} = <c_{i,up}c_{i,dn}>

            complex<double> D_i = Vp_ * Pair_OPs_[x];
            e_i_up = MFParam_TL_.getRowColIndex(site_prow, 0);
            e_i_dn = MFParam_TL_.getRowColIndex(site_prow, 1);
            
            C_mat[e_i_dn][e_i_up + Hsize] += (-1.0)*One_Complex*D_i;
            C_mat[e_i_up][e_i_dn + Hsize] +=  (1.0)*One_Complex*D_i;
            
            C_mat[e_i_up + Hsize][e_i_dn] += (-1.0)*One_Complex*conj(D_i);
            C_mat[e_i_dn + Hsize][e_i_up] +=  (1.0)*One_Complex*conj(D_i);
        }
    }
    else if(Parameters_TL_.pairingType==1){//NN triplet pairing
        int e_i_up,e_i_dn,e_j_up,e_j_dn;
        int e_i_s,e_j_s;
        for(int x=0;x<Pair_OPs_size;x++){
            int prow_ind = MFParam_TL_.row_vector_Pair[x];
            int pcol_ind = MFParam_TL_.col_vector_Pair[x];
            assert(prow_ind <= pcol_ind);
            assert(prow_ind<Hsize && pcol_ind>=Hsize);
            int spin_prow = prow_ind / nsites;
            int spin_pcol = (pcol_ind - Hsize) / nsites;
            
            int e_i = prow_ind % nsites;
            int e_j = (pcol_ind - Hsize) % nsites;
            //double bond_angle = getbondAngleNNPairing(e_i,e_j);
            
            if(spin_prow==spin_pcol){
                //Terms we need for chiral p-wave: up-up and dn-dn
                //a. -D_{ij}^{up-up}(c^{dag}_{j,up}c^{dag}_{i,up} - c^{dag}_{i,up}c^{dag}_{j,up}) + h.c.
                //b. -D_{ij}^{dn-dn}(c^{dag}_{j,dn}c^{dag}_{i,dn} - c^{dag}_{i,dn}c^{dag}_{j,dn}) + h.c.
                
                //where D_{ij}^{s,s} = <c_{i,s}c_{j,s}> = |D_{ij}| exp(i*theta_{ij})

                //complex<double> Dij_ss = norm(Pair_OPs_[x])*polar(1.0,bond_angle);
                complex<double> Dij_ss = Vp_ * Pair_OPs_[x];
                e_i_s = MFParam_TL_.getRowColIndex(e_i,spin_prow);
                e_j_s = MFParam_TL_.getRowColIndex(e_j,spin_prow);
                
                C_mat[e_j_s][e_i_s + Hsize] += (-1.0)*One_Complex*Dij_ss;
                C_mat[e_i_s][e_j_s + Hsize] +=  (1.0)*One_Complex*Dij_ss;
                
                C_mat[e_i_s + Hsize][e_j_s] += (-1.0)*One_Complex*conj(Dij_ss);
                C_mat[e_j_s + Hsize][e_i_s] +=  (1.0)*One_Complex*conj(Dij_ss);
            }
            else{
                //Terms we need for chiral p-wave: up-dn
                //a. -D_{ij}(c^{dag}_{j,dn}c^{dag}_{i,up} + c^{dag}_{j,up}c^{dag}_{i,dn})
                //b.  D_{ij}(c^{dag}_{i,up}c^{dag}_{j,dn} + c^{dag}_{i,dn}c^{dag}_{j,up})
                //c. -conj(D_{ij})(c_{i,up}c_{j,dn} + c_{i,dn}c_{j,up})
                //d.  conj(D_{ij})(c_{j,dn}c_{i,up} + c_{j,up}c_{i,dn})
                
                //where D_{ij} = (<c_{i,up}c_{j,dn}> - <c_{i,dn}c_{j,up}>) = |D_{ij}| exp(i*theta_{ij})
                
                assert(spin_prow==0 && spin_pcol==1);
                //complex<double> Dij = norm(Pair_OPs_[x])*polar(1.0,bond_angle);
                complex<double> Dij = 0.5 * Vp_ * Pair_OPs_[x];

                e_i_up = MFParam_TL_.getRowColIndex(e_i,0);
                e_i_dn = MFParam_TL_.getRowColIndex(e_i,1);
                e_j_up = MFParam_TL_.getRowColIndex(e_j,0);
                e_j_dn = MFParam_TL_.getRowColIndex(e_j,1);
                
                C_mat[e_j_dn][e_i_up + Hsize] += (-1.0)*One_Complex*Dij;
                C_mat[e_j_up][e_i_dn + Hsize] += (-1.0)*One_Complex*Dij;
                C_mat[e_i_up][e_j_dn + Hsize] +=  (1.0)*One_Complex*Dij;
                C_mat[e_i_dn][e_j_up + Hsize] +=  (1.0)*One_Complex*Dij;
                
                C_mat[e_i_up + Hsize][e_j_dn] += (-1.0)*One_Complex*conj(Dij);
                C_mat[e_i_dn + Hsize][e_j_up] += (-1.0)*One_Complex*conj(Dij);
                C_mat[e_j_dn + Hsize][e_i_up] +=  (1.0)*One_Complex*conj(Dij);
                C_mat[e_j_up + Hsize][e_i_dn] +=  (1.0)*One_Complex*conj(Dij);
            }
        }
    }
    else if(Parameters_TL_.pairingType==2){//NN singlet pairing
        int e_i_up,e_i_dn,e_j_up,e_j_dn;
        for(int x=0;x<Pair_OPs_size;x++){
            int prow_ind = MFParam_TL_.row_vector_Pair[x];
            int pcol_ind = MFParam_TL_.col_vector_Pair[x];
            assert(prow_ind <= pcol_ind);
            assert(prow_ind<Hsize && pcol_ind>=Hsize);
            int spin_prow = prow_ind / nsites;
            int spin_pcol = (pcol_ind - Hsize) / nsites;
            assert(spin_prow==0 && spin_pcol==1);

            
            int e_i = prow_ind % nsites;
            int e_j = (pcol_ind - Hsize) % nsites;
            //double bond_angle = 2.0*getbondAngleNNPairing(e_i, e_j);
            
            //complex<double> Dij = norm(Pair_OPs_[x])*polar(1.0,bond_angle);
            complex<double> Dij = 0.5* Vp_ * Pair_OPs_[x];
            e_i_up = MFParam_TL_.getRowColIndex(e_i,0);
            e_i_dn = MFParam_TL_.getRowColIndex(e_i,1);
            e_j_up = MFParam_TL_.getRowColIndex(e_j,0);
            e_j_dn = MFParam_TL_.getRowColIndex(e_j,1);

            //Terms we need for chiral d-wave:
            //a. -D_{ij}(c^{dag}_{j,dn}c^{dag}_{i,up} - c^{dag}_{j,up}c^{dag}_{i,dn})
            //b.  D_{ij}(c^{dag}_{i,up}c^{dag}_{j,dn} - c^{dag}_{i,dn}c^{dag}_{j,up})
            //c. -conj(D_{ij})(c_{i,up}c_{j,dn} - c_{i,dn}c_{j,up})
            //d.  conj(D_{ij})(c_{j,dn}c_{i,up} - c_{j,up}c_{i,dn})
            
            //where D_{ij} = (<c_{i,up}c_{j,dn}> - <c_{i,dn}c_{j,up}>) = |D_{ij}| exp(2*i*theta_{ij})
            
            C_mat[e_j_dn][e_i_up + Hsize] += (-1.0)*One_Complex*Dij;
            C_mat[e_j_up][e_i_dn + Hsize] +=  (1.0)*One_Complex*Dij;
            C_mat[e_i_up][e_j_dn + Hsize] +=  (1.0)*One_Complex*Dij;
            C_mat[e_i_dn][e_j_up + Hsize] += (-1.0)*One_Complex*Dij;

            C_mat[e_i_up + Hsize][e_j_dn] += (-1.0)*One_Complex*conj(Dij);
            C_mat[e_i_dn + Hsize][e_j_up] +=  (1.0)*One_Complex*conj(Dij);
            C_mat[e_j_dn + Hsize][e_i_up] +=  (1.0)*One_Complex*conj(Dij);
            C_mat[e_j_up + Hsize][e_i_dn] += (-1.0)*One_Complex*conj(Dij);
        }
    }
    else{
        throw std::logic_error("Pairing ansatz ill-defined!");
    }
    
    //Antisymmetry check for pairing OPs: D_{ij} = - D_{ji}
    for(int i=0;i<Hsize;i++){
        for(int j=0;j<Hsize;j++){
            assert( std::abs(C_mat[i][j+Hsize] + C_mat[j][i+Hsize]) < 1e-12);
        }
    }

}

void Hamiltonian_TL::Diagonalizer(){

    //cout<<"Starting the diagonalizer"<<endl;
    std::vector<std::complex<double>> Ham_(BdGsize*BdGsize);
    //std::vector<lapack_complex_double> Ham_(H_size*H_size);

    
    //#pragma omp parallel for default(shared) 
    for (int i=0;i<BdGsize;i++) {
        for (int j=0;j<BdGsize;j++) {
            //Ham_[i*H_size + j] = lapack_make_complex_double(C_mat[i][j].real(), C_mat[i][j].imag());
            Ham_[j*BdGsize + i] = C_mat[i][j];
        }
    }

    //cout<<"C_mat copied in Ham"<<endl;

    // LAPACK routine variables
    char jobz = 'V'; // Computing both eigenvalues and eigenvectors
    char uplo = 'L'; // Using the lower triangular part of the matrix
    int n = BdGsize;
    int lda = BdGsize;
    int info;

    std::vector<double> eigs_(BdGsize);

    std::vector<std::complex<double>> work(1);
    //std::vector<lapack_complex_double> work(1);
    std::vector<double> rwork(3 * n - 2);
    int lwork = -1;

    //querying with lwork=-1
    zheev_(&jobz, &uplo, &n, Ham_.data(), &lda, eigs_.data(), work.data(), &lwork, rwork.data(), &info);

    //cout<<"Diagonalization begins"<<endl;
    //Set optimal workspace size
    lwork = static_cast<int>(work[0].real());
    //lwork = static_cast<int>(lapack_complex_double_real(work[0]));
    work.resize(lwork);

    // Perform the eigenvalue decomposition
    zheev_(&jobz, &uplo, &n, Ham_.data(), &lda, eigs_.data(), work.data(), &lwork, rwork.data(), &info);

    // Check for successful execution
    if (info != 0) {
        std::cerr << "LAPACK zheev_ failed with info=" << info << std::endl;
        return;
    }
    //cout<<"Diagonalization finished"<<endl;

    for(int i=0;i<BdGsize;i++){
        evals_[i] = eigs_[i];
    }
    eigs_.clear();

    for(int i=0;i<BdGsize;i++){//i->component-index
        for(int n=0;n<BdGsize;n++){//n->band-index
            evecs_[i][n] = Ham_[n*BdGsize+i]; //evecs_[i][n] = Psi_n(i)
        }
    }
    Ham_.clear();
    for(int i=0;i<Hsize;i++){ //Psi_n = (u_n , v_n)
        for(int n=0;n<BdGsize;n++){
            uvecs_[i][n] = evecs_[i][n];
            vvecs_[i][n] = evecs_[i+Hsize][n];
        }
    }
    //Sanity check: Sum_{n,i}[|u_n(i)|^2 + |v_n(i)|^2 =1]
    for(int n=0;n<BdGsize;n++){
        double val = 0;
        for(int i=0;i<Hsize;i++){
            val += norm(uvecs_[i][n]) + norm(vvecs_[i][n]);
        }
        assert( abs(val - 1.0) < 1e-6 );
    }

    string Evals_out="Eigenvalues.txt";
    ofstream Evals_file_out(Evals_out.c_str());

    for(int n=0;n<evals_.size();n++){
        Evals_file_out<<n<<"    "<<evals_[n]<<endl;
    }

    //Finding index for positive evals: E[n]>=0
    const double Etol = 1e-12;
    int n0=0;
    while(n0<BdGsize && evals_[n0]<-Etol) n0++;
    posN = n0;
}

void Hamiltonian_TL::performDiagonalization(){
    for(int i=0;i<BdGsize;i++){
        for(int j=0;j<BdGsize;j++){
            C_mat[i][j] = Zero_Complex;
        }
    }

    addHoppingTerms();
    addInteractionTerms();
    addPairingTerms();

    //Checking Hamiltonian Hermiticity:
    for(int i=0;i<BdGsize;i++){
        for(int j=0;j<BdGsize;j++){
            assert(abs(C_mat[i][j] - conj(C_mat[j][i]))<1e-12);
        }
    }

    Diagonalizer();
    getNewOPs();
}

double Hamiltonian_TL::fermifunction(double en_){
    const double temp = Parameters_TL_.temp_;

    double x = en_/temp;
    if(x>50.0) return 0.0; //exp(50)~10^21
    if(x<-50.0) return 1.0;//exp(-50)~10^{-22}

    return (1.0/(1.0 + exp(x)));
}

void Hamiltonian_TL::getNewOPs(){
    Mat_1_doub ffn_,omffn_;
    ffn_.resize(BdGsize);  omffn_.resize(BdGsize);
    for(int n=posN;n<BdGsize;n++){
        ffn_[n] = fermifunction(evals_[n]);
        omffn_[n] = 1.0-ffn_[n];
    }

    for(int x=0;x<HF_OPs_size;x++){     new_HF_OPs_[x] = Zero_Complex;      }
    for(int x=0;x<Pair_OPs_size;x++){   new_Pair_OPs_[x] = Zero_Complex;    }


    //Obtaining new HF OPs:<c^{dag}_{row}c_{col}>
    for(int x=0;x<HF_OPs_size;x++){
        int row = MFParam_TL_.row_vector_HF[x];
        int col = MFParam_TL_.col_vector_HF[x];
        complex<double> val=Zero_Complex;

        for(int n=posN;n<BdGsize;n++){
            val += ( conj(uvecs_[row][n])*uvecs_[col][n] )*ffn_[n] + ( conj(vvecs_[row][n])*vvecs_[col][n] )*omffn_[n];
        }
        new_HF_OPs_[x] = val;
    }

    //Obtaining new Pair OPs:
    if(Parameters_TL_.pairingType==0){
        //OP[x] = <c_{i,up}c_{i,dn}>
        for(int x=0;x<Pair_OPs_size;x++){
            int e_i_up = MFParam_TL_.row_vector_Pair[x];          //i,up
            int e_i_dn = MFParam_TL_.col_vector_Pair[x] - Hsize;  //i,dn
            complex<double> val = Zero_Complex;

            for(int n=posN;n<BdGsize;n++){
                val += ( ( uvecs_[e_i_up][n]*conj(vvecs_[e_i_dn][n]) )*omffn_[n] + ( conj(vvecs_[e_i_up][n])*uvecs_[e_i_dn][n] )*ffn_[n] );
            }
            new_Pair_OPs_[x] = val;
        }
    }
    else if(Parameters_TL_.pairingType==1){
        //a. OP_upup[x] = <c_{i,up}c_{j,up}> , OP_dndn[x] = <c_{i,dn}c_{j,dn}>
        //b. OP_updn[x] = <c_{i,up}c_{j,dn}> + <c_{i,dn}c_{j,up}>
        for(int x=0;x<Pair_OPs_size;x++){
            int prow_ind = MFParam_TL_.row_vector_Pair[x];
            int pcol_ind = MFParam_TL_.col_vector_Pair[x];
            
            int spin_prow = prow_ind / nsites;
            int spin_pcol = (pcol_ind - Hsize) / nsites;
            
            int e_i = prow_ind % nsites;
            int e_j = (pcol_ind - Hsize) % nsites;

            if(spin_prow==spin_pcol){
                int e_i_s = MFParam_TL_.getRowColIndex(e_i,spin_prow);
                int e_j_s = MFParam_TL_.getRowColIndex(e_j,spin_prow);

                complex<double> val = Zero_Complex;
                for(int n=posN;n<BdGsize;n++){
                    val += ( uvecs_[e_i_s][n]*conj(vvecs_[e_j_s][n]) )*omffn_[n] + ( conj(vvecs_[e_i_s][n])*uvecs_[e_j_s][n] )*ffn_[n];
                }
                new_Pair_OPs_[x] = val;
            }
            else{
                assert(spin_prow==0 && spin_pcol==1);
                int e_i_up = MFParam_TL_.getRowColIndex(e_i,spin_prow);
                int e_i_dn = MFParam_TL_.getRowColIndex(e_i,spin_pcol);
                int e_j_up = MFParam_TL_.getRowColIndex(e_j,spin_prow);
                int e_j_dn = MFParam_TL_.getRowColIndex(e_j,spin_pcol);

                complex<double> val = Zero_Complex;
                for(int n=posN;n<BdGsize;n++){
                    val += 0.5*( ( uvecs_[e_i_up][n]*conj(vvecs_[e_j_dn][n]) +  uvecs_[e_i_dn][n]*conj(vvecs_[e_j_up][n]) )*omffn_[n]
                            + ( conj(vvecs_[e_i_up][n])*uvecs_[e_j_dn][n] + conj(vvecs_[e_i_dn][n])*uvecs_[e_j_up][n] )*ffn_[n] );
                }
                new_Pair_OPs_[x] = val;
            }
        }
    }
    else if(Parameters_TL_.pairingType==2){
        //OP_[x] = <c_{i,up}c_{j,dn}> - <c_{i,dn}c_{j,up}>
        for(int x=0;x<Pair_OPs_size;x++){
            int prow_ind = MFParam_TL_.row_vector_Pair[x];
            int pcol_ind = MFParam_TL_.col_vector_Pair[x];

            int spin_prow = prow_ind / nsites;
            int spin_pcol = (pcol_ind - Hsize) / nsites;
            assert(spin_prow==0 && spin_pcol==1);
            
            int e_i = prow_ind % nsites;
            int e_j = (pcol_ind - Hsize) % nsites;

            int e_i_up = MFParam_TL_.getRowColIndex(e_i,spin_prow);
            int e_i_dn = MFParam_TL_.getRowColIndex(e_i,spin_pcol);
            int e_j_up = MFParam_TL_.getRowColIndex(e_j,spin_prow);
            int e_j_dn = MFParam_TL_.getRowColIndex(e_j,spin_pcol);

            complex<double> val = Zero_Complex;
            for(int n=posN;n<BdGsize;n++){
                val += 0.5*( ( uvecs_[e_i_up][n]*conj(vvecs_[e_j_dn][n]) -  uvecs_[e_i_dn][n]*conj(vvecs_[e_j_up][n]))*omffn_[n]
                        + ( conj(vvecs_[e_i_up][n])*uvecs_[e_j_dn][n] - conj(vvecs_[e_i_dn][n])*uvecs_[e_j_up][n] )*ffn_[n] );
            }
            new_Pair_OPs_[x] = val;
        }
    }
    else{
    }
    
    //Generating new OPs:
    for(int i=0;i<HF_OPs_size;i++){
        newOPs_[i] = new_HF_OPs_[i];
    }
    for(int i=0;i<Pair_OPs_size;i++){
        newOPs_[i + HF_OPs_size] = new_Pair_OPs_[i];
    }
    assert(newOPs_.size() == OPs_size);
}

double Hamiltonian_TL::getOPError(){
    double op_error = 0.0;
    for(int x=0;x<OPs_size;x++){
        op_error += norm(OPs_[x]-newOPs_[x]);
    }
    return op_error;
}

void Hamiltonian_TL::updateOrderParamsBroyden(int iter){
    /*For Broyden mixing, we use the multidimensional quasi-Newton-Raphson method:->    
        X^{m+1}_{in} = X^{m}_{in} - B^{m}*F^{m};

        where, based on Sherman-Morrison formula:
        B^{m+1} = B^{m} + [ |d(X_{in})> - B^{m}|dF> ] * <d(X_{in})|B^{m};
                             ----------------------
                              <d(X_{in})|B^{m}|dF>

        with B^{0} = alpha*I;   and, |dF> = F^{m+1} - F^{m};
        |d(X_{in})> = X^{m+1}_{in} - X^{m}_{in};
    */

    //Generating full-rank row and col vectors:
    Mat_1_int row_full, col_full;
    row_full.clear();   col_full.clear();
    row_full.reserve(HF_OPs_size+Pair_OPs_size);
    col_full.reserve(HF_OPs_size+Pair_OPs_size);

    for(int x=0;x<HF_OPs_size;x++){
        row_full.push_back(MFParam_TL_.row_vector_HF[x]);
        col_full.push_back(MFParam_TL_.col_vector_HF[x]);
    }
    for(int x=0;x<Pair_OPs_size;x++){
        row_full.push_back(MFParam_TL_.row_vector_Pair[x]);
        col_full.push_back(MFParam_TL_.col_vector_Pair[x]);
    }
    assert(row_full.size() == OPs_size);


    //Flattening the complex OP vectors into real vectors:->
    int flat_OPs_size;
    Mat_1_int flat_index;
    flat_index.clear();
    //First OPs_size -> real parts; remaining -> imag part of Fock
    for(int x=0;x<OPs_size;x++){
        flat_index.push_back(x);
    }
    for(int x=0;x<OPs_size;x++){
        if(row_full[x] != col_full[x]){
            flat_index.push_back(x);
        }
    }
    flat_OPs_size = flat_index.size();

//    cout<<"flat OPs size = "<<flat_OPs_size<<endl;

    F_m.resize(flat_OPs_size);      F_mp1.resize(flat_OPs_size);
    del_F.resize(flat_OPs_size);    del_X.resize(flat_OPs_size);

    B_m.resize(flat_OPs_size);      B_mp1.resize(flat_OPs_size);
    for(int i=0;i<flat_OPs_size;i++){
        B_m[i].resize(flat_OPs_size);   B_mp1[i].resize(flat_OPs_size);
    }

    Mat_1_doub Left_vec,Right_vec;
    Left_vec.resize(flat_OPs_size);       Right_vec.resize(flat_OPs_size);
    double amplitude;
    if(iter == 0){
        //Initializing the approximated inverse Jacobian matrix (B):->
        for(int i=0;i<flat_OPs_size;i++){
            for(int j=0;j<flat_OPs_size;j++){
                if(i==j){
                    B_mp1[i][j] = -alpha;
                }
                else{
                    B_mp1[i][j] = 0.0;
                }
            }
        }

        //Calculating difference in OPs:->
        for(int i=0;i<flat_OPs_size;i++){
            int idx = flat_index[i];
            if(i<OPs_size){
                F_mp1[i] = (newOPs_[idx] - OPs_[idx]).real();
            }
            else{
                F_mp1[i] = (newOPs_[idx] - OPs_[idx]).imag();
            }
        }

        //Computing delX:-> (delX - Bm*Fm)
        for(int i=0;i<flat_OPs_size;i++){
            del_X[i]=0.0;
            for(int j=0;j<flat_OPs_size;j++){
                del_X[i] += -B_mp1[i][j]*F_mp1[j];
            }
        }

        //Getting new estimates of order parameters:->
        for(int i=0;i<flat_OPs_size;i++){
            int idx = flat_index[i];
            if(i<OPs_size){
                newOPs_[idx].real(OPs_[idx].real() + del_X[i]);
            }
            else{
                newOPs_[idx].imag(OPs_[idx].imag() + del_X[i]);
            }
        }

        //Saving F_m and B_m for next iteration:->
        F_m=F_mp1;
        B_m=B_mp1;
    }
    else{
        //Obtaining new F-vetor:->
        for(int i=0;i<flat_OPs_size;i++){
            int idx = flat_index[i];
            if(i<OPs_size){
                F_mp1[i] = (newOPs_[idx] - OPs_[idx]).real();
            }
            else{
                F_mp1[i] = (newOPs_[idx] - OPs_[idx]).imag();
            }
        }
        //Calculating del_F vector for calculating new B:->
        for(int i=0;i<flat_OPs_size;i++){
            del_F[i] = F_mp1[i] - F_m[i];
        }

        //Calculating left and right vectors for getting the new B:->
        for(int i=0;i<flat_OPs_size;i++){
            Left_vec[i]=0.0;
            Right_vec[i]=0.0;

            for(int j=0;j<flat_OPs_size;j++){
                Left_vec[i] += -B_m[i][j]*del_F[j];
                Right_vec[i] += B_m[j][i]*del_X[j];
            }
        }

        amplitude=0.0;
        for(int i=0;i<flat_OPs_size;i++){
            amplitude +=-del_X[i]*Left_vec[i];
        }
        if(abs(amplitude) < 1e-10){
            amplitude = 1e-10;
        }

        for(int i=0;i<flat_OPs_size;i++){
            Left_vec[i] = del_X[i]+Left_vec[i];
        }

        for(int i=0;i<flat_OPs_size;i++){
            for(int j=0;j<flat_OPs_size;j++){
                B_mp1[i][j] = B_m[i][j] + (Left_vec[i]*Right_vec[j])/amplitude;
            }
        }

        //Calculating difference in OPs:->
        for(int i=0;i<flat_OPs_size;i++){
            del_X[i]=0.0;
            for (int j=0;j<flat_OPs_size;j++){
                del_X[i] += -B_mp1[i][j]*F_mp1[j];
            }
        }

        //Getting new estimates of order parameters:->
        for(int i=0;i<flat_OPs_size;i++){
            int idx = flat_index[i];
            if(i<OPs_size){
                newOPs_[idx].real(OPs_[idx].real() + del_X[i]);
            }
            else{
                newOPs_[idx].imag(OPs_[idx].imag() + del_X[i]);
            }
        }

        //Saving F_m and B_m for next iteration:->
        F_m=F_mp1;
        B_m=B_mp1;
    }

    for(int x=0;x<HF_OPs_size;x++){
        new_HF_OPs_[x] = newOPs_[x];
    }
    for(int x=0;x<Pair_OPs_size;x++){
        new_Pair_OPs_[x] = newOPs_[x + HF_OPs_size];
    }
}

void Hamiltonian_TL::updateOrderParamsAnderson(int iter){
    //Updating OPs based on Anderson mixing, referenced from : [https://doi.org/10.1137/10078356X]
    /*For Anderson mixing, we use the following algorithm:->
        given x_0 and m>=1;
        set x_1 = g(x_0);
        for k=0,1,2,...{
            if(k<m_k){
                x_{k+1} = x_k + alpha*f_k
            }
            else{
                set F_k=(Df_{k-m_k} ,..., Df_{k-1}) where Df_i = f_{i} - f_{i-1};
                set X_k=(Dx_{k-m_k} ,..., Dx_{k-1}) where Dx_i = x_{i} - x_{i-1};
                determine gamma^{k} = (gamma_0^{k}, ..., gamma_{m_k}^{k})^T that solves: min_{gamma} ||f_k - F_k gamma||_2;
                
                set x_{k+1} = x_k - G_k*f_k
                (where G_k = -I + (X_k + F_k)*({F_k}^t{F_k})^{-1}{F_k}^t)
            }
        }
    */

    for(int x=0;x<OPs_size;x++){
        x_k[x] = OPs_[x];
        f_k[x] = newOPs_[x] - OPs_[x];
    }
    //0th iteration
    if(iter == 0){
        //Simple mixing till we build full rank:
        for(int x=0;x<OPs_size;x++){
            newOPs_[x] = x_k[x] + alpha * f_k[x];
        }
        //storing for next iteration
        x_km1 = x_k;    f_km1 = f_k;
    }
    else{
        for(int x=0;x<OPs_size;x++){
            del_x_k[x] = x_k[x] - x_km1[x];
            del_f_k[x] = f_k[x] - f_km1[x];
        }

        //Building the full-rank [OPs_size x AM_depth] matrices: Fmat and Xmat:
        if(iter >=1 && iter <= AM_depth){
            for(int x=0;x<OPs_size;x++){
                Xmat[x][iter - 1] = del_x_k[x];
                Fmat[x][iter - 1] = del_f_k[x];
            }
            //Simple mixing till we build full-rank:
            for(int x=0;x<OPs_size;x++){
                newOPs_[x] = x_k[x] + alpha * f_k[x];
            }

            //storing x_k, f_k for next iteration:
            x_km1 = x_k;    f_km1 = f_k;
        }
        else{
            //Now we will perform the Anderson Acceleration:
            //adding new vectors in Xmat and Fmat: sliding window!
            Mat_2_Complex_doub Xmat_temp = Xmat;
            Mat_2_Complex_doub Fmat_temp = Fmat;
            for(int x=0;x<OPs_size;x++){
                for(int r=0;r<AM_depth-1;r++){
                    Xmat[x][r] = Xmat_temp[x][r+1];
                    Fmat[x][r] = Fmat_temp[x][r+1];
                }
                Xmat[x][AM_depth-1] = del_x_k[x];
                Fmat[x][AM_depth-1] = del_f_k[x];
            }

            //we will always perform SVD here for the least square minimization to get gamma_k's: min_{g}||f_k - Fmat*g||_2
            //F.gammma = f_k
            //F=U*Sigma*V^{dag} => gamma = V*Sigma^{-1}*U^{dag} f_k

            int y_ = min(OPs_size,AM_depth);
            Mat_2_Complex_doub A_=Fmat; //rank = OPs_size x AM_depth
            Mat_2_Complex_doub VT_, U_; 
            Mat_1_doub Sigma_;  //size = y_
            U_.clear();     VT_.clear();    Sigma_.clear();

            performComplexSVD(A_,VT_,U_,Sigma_);

            Mat_1_Complex_doub UT_times_f(y_,Zero_Complex); //U^{dag} f_k
            for(int j=0;j<y_;j++){
                for(int i=0;i<OPs_size;i++){
                    UT_times_f[j] += conj(U_[i][j]) * f_k[i];
                }
            }
            Mat_1_Complex_doub Sinv_times_UTf(y_,Zero_Complex); //Sigma^{-1}*U^{dag} f_k
            for(int j=0;j<y_;j++){
                if(abs(Sigma_[j])>=0.0001){
                    Sinv_times_UTf[j] = (1.0/Sigma_[j])*UT_times_f[j];
                }
                else{
                    Sinv_times_UTf[j] = Zero_Complex;
                }
            }
            for(int r=0;r<AM_depth;r++){
                gamma_k[r]=Zero_Complex;
                for(int j=0;j<y_;j++){ //V*Sigma^{-1}*U^{dag} f_k
                    gamma_k[r] += conj(VT_[j][r])*Sinv_times_UTf[j];
                }
            }

            //obtaining the OP for next iter: x_{k+1} = x_k + f_k - (Xmat+Fmat)*gamma_k
            for(int x=0;x<OPs_size;x++){
                complex<double> temp_x = Zero_Complex;
                complex<double> temp_f = Zero_Complex;
                for(int r=0;r<AM_depth;r++){
                    temp_x += Xmat[x][r] * gamma_k[r];
                    temp_f += Fmat[x][r] * gamma_k[r];
                }
                new_x_k[x] = x_k[x] - 1.0*temp_x;
                new_f_k[x] = f_k[x] - 1.0*temp_f;
            }
            for(int x=0;x<OPs_size;x++){
                x_kp1[x] = new_x_k[x] + alpha_AA * new_f_k[x];
            }

            //obtaining newOPs, keeping flexibility of damping if required: 
            //to reduce oscillations of Anderson by damping it with simple mixing
            for(int x=0;x<OPs_size;x++){
                newOPs_[x] = (1-beta_damping) * x_kp1[x] + beta_damping * (x_k[x] + alpha * f_k[x]);
            }
            //Keeping the Hartree OPs strictly real:
            for(int x=0;x<OPs_size;x++){
                if(x<HF_OPs_size){
                    if(MFParam_TL_.row_vector_HF[x]==MFParam_TL_.col_vector_HF[x]){
                        newOPs_[x].real(newOPs_[x].real());
                        newOPs_[x].imag(0.0);
                    }
                }
            }

            //storing x_k, f_k for next iteration:
            x_km1 = x_k;    f_km1 = f_k;
        }
    }

    for(int i=0;i<HF_OPs_size;i++){
        new_HF_OPs_[i] = newOPs_[i];        
    }
    for(int i=0;i<Pair_OPs_size;i++){
        new_Pair_OPs_[i] = newOPs_[i + HF_OPs_size];
    }

}

void Hamiltonian_TL::performSVD(Mat_2_doub &A, Mat_2_doub &VT, Mat_2_doub &U, Mat_1_doub &Sigma){
    char jobz='S';
    int m = A.size();
    int n = A[0].size();

    int lda=m;
    int ldu=m;
    int ldvt=n;
    int min_mn=min(m,n);

    //Matrix-A is of rank=nxm
    vector<double> A_flat(m * n);
    for(int r=0;r<m;r++){
        for(int c=0;c<n;c++){
            A_flat[c * m + r] = A[r][c];  // column-major
        }
    }
    /*
    for(int r=0;r<m;r++){
        for(int i=0;i<n;i++){
            A_flat[r *n + i] = A[r][i]; // transpose for column-major
        }
    }
    */

    Sigma.clear();
    Sigma.resize(min(m,n));
    vector<double> U_flat(ldu * min_mn);
    vector<double> VT_flat(ldvt * n);

/*
    U.resize(ldu);
    for(int i=0;i<ldu;i++){
        U[i].resize(min_mn);
    }
    VT.resize(ldvt);
    for(int j=0;j<ldvt;j++){
        VT[j].resize(n);
    }
*/

    vector<double> work(3);
    int info;
    int lwork= -1;
    vector<int> iwork(8*min(m,n));
    //query:
    dgesdd_(&jobz, &m, &n, A_flat.data(),&lda, Sigma.data(),U_flat.data(), &ldu, VT_flat.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);

    if(info != 0){
        std::cerr << "SVD query failed with info = " << info << std::endl;
        return;
    }
    lwork = static_cast<int>(work[0]);
    work.resize(lwork);

    // Actual SVD
    dgesdd_(&jobz, &m, &n, A_flat.data(), &lda, Sigma.data(),U_flat.data(), &ldu, VT_flat.data(), &ldvt, work.data(), &lwork, iwork.data(), &info);

    if(info != 0){
        std::cerr << "SVD failed with info = " << info << std::endl;
        return;
    }
    // Convert U and VT back to 2D row-major
    U.resize(m);
    for(int i = 0; i < m; ++i){
        U[i].resize(min_mn);
        for(int j = 0; j < min_mn; ++j){
            U[i][j] = U_flat[j * m + i]; // column-major back to row-major
        }
    }
    VT.resize(ldvt);
    for(int i = 0; i < ldvt; ++i){
        VT[i].resize(n);
        for(int j = 0; j < n; ++j){
            VT[i][j] = VT_flat[j * ldvt + i];
        }
    }

    VT.resize(min_mn);
    for(int i = 0; i < min_mn; ++i){
        VT[i].resize(n);
        for(int j = 0; j < n; ++j){
            VT[i][j] = VT_flat[j * min_mn + i]; // column-major back to row-major
        }
    }

    VT.resize(min_mn);
for(int i=0; i<min_mn; ++i){
  VT[i].resize(n);
  for(int j=0; j<n; ++j){
    // correct column-major unpacking with ldvt = n
    VT[i][j] = VT_flat[i + j*ldvt];  
  }
}

    /*
    U.resize(ldu);
    for(int i=0;i<ldu;i++){
        U[i].resize(min_mn);
        for(int j=0;j<min_mn; ++j){
            U[i][j] = U_flat[i * min_mn + j];
        }
    }

    VT.resize(ldvt);
    for(int i = 0; i < ldvt; ++i){
        VT[i].resize(n);
        for(int j = 0; j < n; ++j){
            VT[i][j] = VT_flat[i * n + j];
        }
    }
    */
}

void Hamiltonian_TL::performComplexSVD(Mat_2_Complex_doub &A, Mat_2_Complex_doub &VT, Mat_2_Complex_doub &U, Mat_1_doub &Sigma){
    char jobz='S';          //[if we need all cols of U and all rows of VT then put 'A' instead of 'S']
    int M = (int) A.size(); //no of rows
    int N = (int) A[0].size(); //no of cols
    int min_MN = min(M,N);

    //Here: Rank(U) = (LDU,UCOLS); Rank(Sigma) = min(M,N); Rank(VT) = (LDVT,N)
    //A = MxN complex matrix
    //U = MxM unitary matrix (for jobz='A') / Mx[min(M,N)] (for jobz='S')
    //Sigma = MxN real matrix which is zero except for its min(M,N) diagonal elements (for 'A' or 'S')
    //VT = NxN unitary matrix (for jobz='A') / [min(M,N)]xN (for jobz='S')

    int lda=M;  
    int ldu=M;
    int ldvt=min_MN;

    U.clear();  U.resize(ldu);
    for(int i=0;i<ldu;i++){     U[i].resize(min_MN);    }
    VT.clear(); VT.resize(ldvt);
    for(int i=0;i<ldvt;i++){    VT[i].resize(N);        }

    std::vector<std::complex<double>> A_flat(M*N);
    for(int j=0;j<N;j++){  //col-index
        for(int i=0;i<M;i++){  //row-index
            A_flat[j*M + i] = A[i][j];
        }
    }

    Sigma.clear();      Sigma.resize(min_MN);
    std::vector<std::complex<double>> U_flat(ldu * min_MN);
    std::vector<std::complex<double>> VT_flat(ldvt * N);

    vector<std::complex<double>> work(3);
    int info;
    int lwork= -1;
    vector<int> iwork(8*min_MN);
    /*
    RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK))
    If JOBZ = 'N', LRWORK >= 5*min(M,N).
    Otherwise,
    LRWORK >= min(M,N)*max(5*min(M,N)+7,2*max(M,N)+2*min(M,N)+1)
    */
    int lrwork = min_MN * max( 5*min_MN + 7, 2*max(M,N) + 2*min_MN + 1 );
    vector<double> rwork(lrwork);

    //zgesdd_ (char *, int *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, std::complex<double> *, int*, 
    //                    std::complex<double> *, int *, double *, int *, int *);

    //query:
    zgesdd_(&jobz, &M, &N, A_flat.data(), &lda, Sigma.data(), U_flat.data(), &ldu, VT_flat.data(), &ldvt, work.data(), &lwork, rwork.data(), iwork.data(), &info);

    if(info != 0){
        std::cerr << "SVD query failed with info = " << info << std::endl;
        return;
    }
    lwork = static_cast<int>(work[0].real());
    work.resize(lwork);

    // Actual SVD
    zgesdd_(&jobz, &M, &N, A_flat.data(), &lda, Sigma.data(),U_flat.data(), &ldu, VT_flat.data(), &ldvt, work.data(), &lwork, rwork.data(), iwork.data(), &info);

    if(info != 0){
        std::cerr << "SVD failed with info = " << info << std::endl;
        return;
    }

    //Convert U and VT back:
    for(int j=0;j<min_MN;j++){
        for(int i=0;i<ldu;i++){
            U[i][j] = U_flat[j*ldu + i];
        }
    }
    for(int j=0;j<N;j++){
        for(int i=0;i<ldvt;i++){
            VT[i][j] = VT_flat[j*ldvt + i];
        }
    }
}

double Hamiltonian_TL::calculateTotalParticles(){
    double Nup_tot,Ndn_tot;
    Nup_tot=0.0;    Ndn_tot=0.0;
    for(int rx=0;rx<lx_;rx++){
        for(int ry=0;ry<ly_;ry++){
            int r = Neighbors_TL_.getIndex(rx, ry);
            int r_up = Neighbors_TL_.getIndex(rx, ry);
            int r_dn = Neighbors_TL_.getIndex(rx, ry) + nsites;

            for(int n=posN;n<BdGsize;n++){
                double ffn = fermifunction(evals_[n]);
                Nup_tot += norm(uvecs_[r_up][n]) * ffn + norm(vvecs_[r_up][n]) * (1.0 - ffn);
                Ndn_tot += norm(uvecs_[r_dn][n]) * ffn + norm(vvecs_[r_dn][n]) * (1.0 - ffn);
            }

        }
    }
    double Ntot = Nup_tot + Ndn_tot;
    double n_fill = Ntot / (1.0*nsites);
    double p_dop  = (1.0 - n_fill);

    cout<<"Total_Up-Dn_particles=("<<Nup_tot<<","<<Ndn_tot<<")"<<endl;
    cout<<"Electron_filling_and_hole_doping=("<<n_fill<<","<<p_dop<<")"<<endl;

    return Ntot;
}

double Hamiltonian_TL::Lorentzian(double val){
    return ( (1.0/PI)*(eta/(val*val + eta*eta)) );
}

double Hamiltonian_TL::calculateLDOSatSite(int rx, int ry, double w_val){
    double rho=0.0;
    int site = Neighbors_TL_.getIndex(rx, ry);
    //calculating LDOS: rho(r,w) = sum_{n,s}[|u_n(r,s)|^2 delta(w-En) + |v_n(r,s)|^2 delta(w+En)]
    for(int spin=0;spin<nspin;spin++){
        int r = site + spin*nsites;
        for(int n=posN;n<BdGsize;n++){
            double En = evals_[n];
            rho += norm(uvecs_[r][n]) * Lorentzian(w_val-En) + norm(vvecs_[r][n]) * Lorentzian(w_val+En);
        }
    }
    return rho;
}

void Hamiltonian_TL::calculateAvgdLDOS(){
    ofstream file_avg_ldos_out("Avg_LDOS_vs_omega.txt");
    //calculating avgd LDOS: rho(w) = (1/nsites) * sum_{n,r,s}[|u_n(r,s)|^2 delta(w-En) + |v_n(r,s)|^2 delta(w+En)]
    for(int om=0;om<wsize;om++){
        double w = wmin + om*dw;
        double sum = 0.0;
        for(int rx=0;rx<lx_;rx++){
            for(int ry=0;ry<ly_;ry++){
                sum += calculateLDOSatSite(rx, ry, w);
            }
        }
        double avg_ldos = sum/double(nsites);
        file_avg_ldos_out<<w<<"     "<<avg_ldos<<endl;
    }
}

void Hamiltonian_TL::calculateQPI(double w_val){
    double wp_val = 0.001*w_val;
    Mat_2_doub ldos_,mod_ldos;
    ldos_.resize(lx_);  mod_ldos.resize(lx_);
    for(int rx=0;rx<lx_;rx++){
        ldos_[rx].resize(ly_);  mod_ldos[rx].resize(ly_);
        for(int ry=0;ry<ly_;ry++){
            ldos_[rx][ry] = 0.0;    mod_ldos[rx][ry] = 0.0;
        }
    }

    double avg_ldos=0.0;
    for(int rx=0;rx<lx_;rx++){
        for(int ry=0;ry<ly_;ry++){
            ldos_[rx][ry] = calculateLDOSatSite(rx,ry,wp_val);
            avg_ldos += ldos_[rx][ry];
        }
    }
    avg_ldos = avg_ldos/double(nsites);

    //Calculating local density modulations: d_rho(r,w) = rho(r,w) - <rho(w)>_r
    for(int rx=0;rx<lx_;rx++){
        for(int ry=0;ry<ly_;ry++){
            mod_ldos[rx][ry] = ldos_[rx][ry] - avg_ldos;
        }
    }

    string head="QPI_for_w";
    string tail="mV.txt";
    ofstream file_qpi_out(head + to_string(w_val) + tail);

    for(int qx_ind=-lx_/2;qx_ind<lx_/2;qx_ind++){
        double qx = (2.0*PI*qx_ind)/(1.0*lx_);
        for(int qy_ind=-ly_/2;qy_ind<ly_/2;qy_ind++){
            double qy = (2.0*PI*qy_ind)/(1.0*ly_);
            
            complex<double> raw_qpi = Zero_Complex;
            for(int rx=0;rx<lx_;rx++){
                for(int ry=0;ry<ly_;ry++){//if x_r, y_r doesn't work; directly use rx, ry instead
                    double x_real = 1.0*rx + 0.5*ry;
                    double y_real = (sqrt(3.0)/2.0)*ry;

                    raw_qpi += exp(-Iota_Complex*(qx*x_real + qy*y_real)) * mod_ldos[rx][ry] * (1.0/double(nsites));
                }
            }
            file_qpi_out<<qx<<"     "<<qy<<"    "<<norm(raw_qpi)<<endl;    
        }
        file_qpi_out<<endl;
    }

}

void Hamiltonian_TL::calculateLDOSplusQPI(){
    
    Mat_1_doub avg_ldos(wsize,0.0);
    Mat_3_doub ldos_,mod_ldos;
    ldos_.resize(lx_);  mod_ldos.resize(lx_);

    for(int rx=0;rx<lx_;rx++){
        ldos_[rx].resize(ly_);      mod_ldos[rx].resize(ly_);
        for(int ry=0;ry<ly_;ry++){
            ldos_[rx][ry].resize(wsize);    mod_ldos[rx][ry].resize(wsize);
            for(int om=0;om<wsize;om++){
                ldos_[rx][ry][om] = 0.0;        mod_ldos[rx][ry][om] = 0.0;
            }
        }
    }
    //ofstream file_ldos_out("LDOS_real_space.txt");
    ofstream file_avg_ldos_out("Avg_LDOS_vs_omega.txt");

    for(int om=0;om<wsize;om++){
        double w = wmin + om*dw;
        for(int rx=0;rx<lx_;rx++){
            for(int ry=0;ry<ly_;ry++){

                double rho = 0.0;
                int site = Neighbors_TL_.getIndex(rx,ry);
                for(int spin=0;spin<nspin;spin++){
                    int r = site + spin*nsites;

                    for(int n=posN;n<BdGsize;n++){
                        double En = evals_[n];
                        rho += norm(uvecs_[r][n]) * Lorentzian(w-En) + norm(vvecs_[r][n]) * Lorentzian(w+En);
                    }
                }
                ldos_[rx][ry][om] = rho;
                avg_ldos[om] += rho;

                //file_ldos_out<<rx<<"    "<<ry<<"    "<<w<<"     "<<ldos_[rx][ry][om]<<endl;
            }
        }
        avg_ldos[om] = avg_ldos[om]/(1.0*nsites);
        file_avg_ldos_out<<w<<" "<<avg_ldos[om]<<endl;
    }

    if(Parameters_TL_.singleImpurity){
        //Calculating local density modulations: d_rho(r,w) = rho(r,w) - <rho(w)>_r
        for(int om=0;om<wsize;om++){
            for(int rx=0;rx<lx_;rx++){
                for(int ry=0;ry<ly_;ry++){
                    mod_ldos[rx][ry][om] = ldos_[rx][ry][om] - avg_ldos[om];
                }
            }
        }

        //Calculating QPI: I(q,w) = |d_rho(q,w)|^2 = |sum_r e^{-iq.r}d_rho(r,w)|^2
/*        Mat_3_doub qpi_;
        qpi_.resize(lx_);
        for(int qx_ind=0;qx_ind<lx_;qx_ind++){
            qpi_[qx_ind].resize(ly_);
            for(int qy_ind=0;qy_ind<ly_;qy_ind++){
                qpi_[qx_ind][qy_ind].resize(wsize);
                for(int om=0;om<wsize;om++){
                    qpi_[qx_ind][qy_ind][om] = 0.0;
                }
            }
        }
*/
        //string file_qpi_pm = "QPI_pm3d.txt";
        //ofstream file_qpi_pm_out(file_qpi_pm.c_str());

        string head="QPI_for_w";
        string tail="mV.txt";
        
        for(int om=0;om<wsize;om++){
            double w = wmin + om*dw;
            ofstream file_qpi_out(head + to_string(w) + tail);

            for(int qy_ind=-ly_/2;qy_ind<=ly_/2;qy_ind++){
                double qy = (2.0*PI*qy_ind)/(1.0*ly_);
                for(int qx_ind=-lx_/2;qx_ind<=lx_/2;qx_ind++){
                    double qx = (2.0*PI*qx_ind)/(1.0*lx_);

                    complex<double> raw_qpi = Zero_Complex;
                    for(int rx=0;rx<lx_;rx++){
                        for(int ry=0;ry<ly_;ry++){
                            
                            raw_qpi += exp(-Iota_Complex*(qx*rx + qy*ry)) * mod_ldos[rx][ry][om] * (1.0/(1.0*nsites));
                        }
                    }
                    file_qpi_out<<qx<<"     "<<qy<<"    "<<norm(raw_qpi)<<endl;
                    //qpi_[qx_ind][qy_ind][om] = norm(raw_qpi);
                    //file_qpi_out<<qx<<"     "<<qy<<"     "<<qpi_[qx_ind][qy_ind][om]<<endl;
                    //file_qpi_pm_out<<qx<<"  "<<qy<<"    "<<w<<"     "<<qpi_[qx_ind][qy_ind][om]<<endl;
                }
                //file_qpi_pm_out<<endl;
            }
        }

    }    

}

void Hamiltonian_TL::calculateQPIusingTmatrix(double w_val){

}

#endif
