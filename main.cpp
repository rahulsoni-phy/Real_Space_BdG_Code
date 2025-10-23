#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <complex>
#include <cmath>
#include <cassert>
#include <utility>

#include "tensors.h"
#include "Parameters_TL.h"
#include "Neighbors_TL.h"
#include "MFParameters_TL.h"
#include "Hamiltonian_TL.h"
//#include "Observables_TL.h"

using namespace std;

int main(int argc, char *argv[]){
    string inputfile = argv[1];
    
    Parameters_TL Parameters_TL_;
    Parameters_TL_.Initialize(inputfile);
    
    Neighbors_TL Neighbors_TL_(Parameters_TL_);
    Neighbors_TL_.generateInteractionPairs();

    mt19937_64 Generator_(Parameters_TL_.Seed_);

    MFParam_TL MFParam_TL_(Parameters_TL_,Neighbors_TL_,Generator_);
    Hamiltonian_TL Hamiltonian_TL_(Parameters_TL_,Neighbors_TL_,MFParam_TL_);

    //----------------------------------------------//
    // Initializing the order parameters:----->
    Mat_1_Complex_doub OPs_HF_Old,OPs_Pair_Old;
    OPs_HF_Old = MFParam_TL_.initial_HF_OPs_;
    OPs_Pair_Old = MFParam_TL_.initial_Pair_OPs_;

    Mat_1_Complex_doub OPs_old;
    OPs_old.resize(Hamiltonian_TL_.OPs_size);
    for(int x=0;x<Hamiltonian_TL_.HF_OPs_size;x++){
        OPs_old[x] = OPs_HF_Old[x];
    }
    for(int x=0;x<Hamiltonian_TL_.Pair_OPs_size;x++){
        OPs_old[x+Hamiltonian_TL_.HF_OPs_size] = OPs_Pair_Old[x];
    }
    
    //----------------------------------------------//
    ofstream file_out_progress("output_selfconsistency.txt");
    if (!file_out_progress){
        cerr << "Failed to open output file." << endl;
        return 1;
    }
    // Self-consistency loop:---->
    cout << "target error = " << Parameters_TL_.Conv_err << endl;
    cout << "Max iterations = " << Parameters_TL_.Max_iters << endl;

    int iter = 0;
    double OP_error;
    OP_error = 10.0;
    int warmup=40;

    while((OP_error > Parameters_TL_.Conv_err) && (iter < Parameters_TL_.Max_iters)){
        OP_error = 0.0;
        
        Hamiltonian_TL_.HF_OPs_ = OPs_HF_Old;
        Hamiltonian_TL_.Pair_OPs_ = OPs_Pair_Old;
        Hamiltonian_TL_.OPs_ = OPs_old;

        Hamiltonian_TL_.performDiagonalization();

        if (Parameters_TL_.Simple_mixing == true){
            double alpha = Parameters_TL_.alpha_;
            
            for(int x = 0; x < Hamiltonian_TL_.OPs_size; x++){
                Hamiltonian_TL_.newOPs_[x] = (1 - alpha) * Hamiltonian_TL_.OPs_[x] + alpha * Hamiltonian_TL_.newOPs_[x];
            }
            for(int x=0; x<Hamiltonian_TL_.HF_OPs_size;x++){
                Hamiltonian_TL_.new_HF_OPs_[x] = Hamiltonian_TL_.newOPs_[x];
            }
            for(int x=0; x<Hamiltonian_TL_.Pair_OPs_size;x++){
                Hamiltonian_TL_.new_Pair_OPs_[x] = Hamiltonian_TL_.newOPs_[x+Hamiltonian_TL_.HF_OPs_size];
            }
        }
        else{
            //Hamiltonian_TL_.updateOrderParamsBroyden(iter);
            Hamiltonian_TL_.updateOrderParamsAnderson(iter);
        }
        
        OP_error = Hamiltonian_TL_.getOPError();
        file_out_progress<<iter<<"  "<<OP_error<<endl;

        OPs_old = Hamiltonian_TL_.newOPs_;
        OPs_HF_Old = Hamiltonian_TL_.new_HF_OPs_;
        OPs_Pair_Old = Hamiltonian_TL_.new_Pair_OPs_;

        iter++;
    }

    string File_out_HFOPs="Final_HF_OPs.txt";
    ofstream file_out_HFOPs(File_out_HFOPs.c_str());
    for(int x=0;x<Hamiltonian_TL_.HF_OPs_size;x++){
        file_out_HFOPs<<x<<"    "<<Hamiltonian_TL_.new_HF_OPs_[x].real()<<"     "<<Hamiltonian_TL_.new_HF_OPs_[x].imag()<<endl;
    }

    string File_out_PairOPs="Final_Pair_OPs.txt";
    ofstream file_out_PairOPs(File_out_PairOPs.c_str());
    for(int x=0;x<Hamiltonian_TL_.Pair_OPs_size;x++){
        file_out_PairOPs<<x<<"    "<<Hamiltonian_TL_.new_Pair_OPs_[x].real()<<"     "<<Hamiltonian_TL_.new_Pair_OPs_[x].imag()<<endl;
    }

    if(Parameters_TL_.get_qpi){
        //Hamiltonian_TL_.calculateLDOSplusQPI();
        Hamiltonian_TL_.calculateQPIonly(0.0);
        Hamiltonian_TL_.calculateQPIonly(1.0);
        Hamiltonian_TL_.calculateQPIonly(-1.0);
    }

    return 0;
}