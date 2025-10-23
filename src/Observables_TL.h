#include "tensors.h"
#include "Parameters_TL.h"
#include "Neighbors_TL.h"
#include "Hamiltonian_TL.h"

#ifndef Observables_TL_class
#define Observables_TL_class

class Observables_TL{
    public:
        //Constructor:
        Observables_TL(Parameters_TL &Parameters_TL__, Neighbors_TL &Neighbors_TL__, Hamiltonian_TL &Hamiltonian_TL__): Parameters_TL_(Parameters_TL__), Neighbors_TL_(Neighbors_TL__), Hamiltonian_TL_(Hamiltonian_TL__){
        initialize();
    }

    Parameters_TL &Parameters_TL_;
    Neighbors_TL &Neighbors_TL_;
    Hamiltonian_TL &Hamiltonian_TL_;

    void initialize();


    int lx_,ly_,block_size_,H_size;

    double one_by_PI_=1.0/(1.0*PI);
    double w_min,w_max,dw,eta;
    int w_size;

    double mu;

    Mat_1_Complex_doub new_OPs_;

    //For Broyden:
    int OP_size;
    Mat_2_doub B_m, B_mp1;
    Mat_1_doub del_X,del_F;
    Mat_1_doub F_m,F_mp1;
};

void Observables_TL::initialize(){
    lx_ = Parameters_TL_.Lx_;
    ly_ = Parameters_TL_.Ly_;
    H_size = Parameters_TL_.BdG_Size;

    w_min = Parameters_TL_.w_min;
    w_max = Parameters_TL_.w_max;
    dw = Parameters_TL_.dw_;
    eta = Parameters_TL_.eta_;

    w_size = (int) ( (w_max - w_min)/dw );

}



#endif