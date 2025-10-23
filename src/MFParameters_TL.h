#include <sstream>
#include "tensors.h"
#include "Parameters_TL.h"
#include "Neighbors_TL.h"
#include "random"

#ifndef MFParam_TL_class
#define MFParam_TL_class

class MFParam_TL{

public:
	// Constructor for the OP class
	MFParam_TL(Parameters_TL &Parameters_TL__, Neighbors_TL &Neighbors_TL__, mt19937_64 &Generator1__)
		: Parameters_TL_(Parameters_TL__), Neighbors_TL_(Neighbors_TL__), Generator1_(Generator1__)	{
		Initialize();
		initialHFOrderParams();
		initialPairingOrderParams();
	}

	Parameters_TL &Parameters_TL_;
	Neighbors_TL &Neighbors_TL_;

	void Initialize();
	double random1(); // returns random in range: [0,1)
	int getRowColIndex(int ind, int spin);
	void initialHFOrderParams();
	void initialPairingOrderParams();
	

	uniform_real_distribution<double> dis1_;
	mt19937_64 &Generator1_;

	int lx_, ly_, nspin;
	int Hsize, nsites, BdGsize;
	
	Mat_1_Complex_doub initial_HF_OPs_,initial_Pair_OPs_;
	Mat_1_int row_vector_HF,col_vector_HF,Super_to_flatind_HF;
	Mat_1_int row_vector_Pair,col_vector_Pair,Super_to_flatind_Pair;
};

void MFParam_TL::Initialize(){
	lx_ = Parameters_TL_.Lx_;
	ly_ = Parameters_TL_.Ly_;

	nspin = 2;
	nsites = Parameters_TL_.Total_Sites;
	Hsize = nspin * nsites;
	BdGsize = Parameters_TL_.BdG_Size;

	Super_to_flatind_HF.resize(Hsize*Hsize,-1);
	Super_to_flatind_Pair.resize(BdGsize*BdGsize,-1);
}

double MFParam_TL::random1(){
	return dis1_(Generator1_);
}

int MFParam_TL::getRowColIndex(int ind, int spin){
	int rx = ind % lx_;
	int ry = ind / lx_;

	int index = spin*nsites + rx + lx_*ry;
	return index;
}

void MFParam_TL::initialHFOrderParams(){
	initial_HF_OPs_.clear();
	row_vector_HF.clear();			col_vector_HF.clear();
	
	if(Parameters_TL_.read_OP_==false){
		string file_OP_="initial_HF_OPs_file.txt";
		ofstream file_out_OP(file_OP_.c_str());
		file_out_OP<<"#row_index    col_index   initial_OP_"<<endl;

		for(int ind1=0;ind1<nsites;ind1++){
			for(int ind2=0;ind2<nsites;ind2++){
				if(Neighbors_TL_.interaction_pairs[ind1][ind2]!=0.0){
					
					//generating onsite Hartree-Fock OPs:
					if(ind2==ind1){
						if(Parameters_TL_.useHartree==true){//onsite Hartree:
							for(int spin=0;spin<nspin;spin++){
								int row = getRowColIndex(ind1,spin);
								int col = row;
								int flat_ind = row + Hsize*col;
								complex<double> val_op(random1(),0.0);
								
								initial_HF_OPs_.push_back(val_op);
								row_vector_HF.push_back(row);
								col_vector_HF.push_back(col);
								Super_to_flatind_HF[flat_ind] = initial_HF_OPs_.size() - 1;
								assert(Super_to_flatind_HF[flat_ind] >= 0);
								file_out_OP<<row<<"             "<<col<<"       "<<initial_HF_OPs_[Super_to_flatind_HF[flat_ind]].real()
												<<"             "<<initial_HF_OPs_[Super_to_flatind_HF[flat_ind]].imag()<<endl;
							}
						}
						if(Parameters_TL_.useFock==true){//onsite Fock:
							int row = getRowColIndex(ind1,0);
							int col = getRowColIndex(ind1,1);
							int flat_ind = row + Hsize*col;
							complex<double> val_op(random1(), random1());
							
							initial_HF_OPs_.push_back(val_op);
							row_vector_HF.push_back(row);
							col_vector_HF.push_back(col);
							Super_to_flatind_HF[flat_ind] = initial_HF_OPs_.size() - 1;
							assert(Super_to_flatind_HF[flat_ind] >= 0);
							file_out_OP<<row<<"             "<<col<<"       "<<initial_HF_OPs_[Super_to_flatind_HF[flat_ind]].real()
											<<"             "<<initial_HF_OPs_[Super_to_flatind_HF[flat_ind]].imag()<<endl;
						}
					}
					else{
						if(Parameters_TL_.useFock==true){//offsite Fock:
							for(int spin=0;spin<nspin;spin++){
								int row = getRowColIndex(ind1,spin);
								
								for(int spin_p=spin;spin_p<nspin;spin_p++){
									int col = getRowColIndex(ind2,spin_p);
									
									if(spin_p==spin){//up-up and dn-dn terms:
										if(ind2>ind1){
											int flat_ind = row + Hsize*col;
											complex<double> val_op(random1(), random1());
											
											initial_HF_OPs_.push_back(val_op);
											row_vector_HF.push_back(row);
											col_vector_HF.push_back(col);
											Super_to_flatind_HF[flat_ind] = initial_HF_OPs_.size() - 1;
											assert(Super_to_flatind_HF[flat_ind] >= 0);
											file_out_OP<<row<<"             "<<col<<"       "<<initial_HF_OPs_[Super_to_flatind_HF[flat_ind]].real()
															<<"             "<<initial_HF_OPs_[Super_to_flatind_HF[flat_ind]].imag()<<endl;
										}
									}
									else{//up-dn terms:
										int flat_ind = row + Hsize*col;
										complex<double> val_op(random1(), random1());
										
										initial_HF_OPs_.push_back(val_op);
										row_vector_HF.push_back(row);
										col_vector_HF.push_back(col);
										Super_to_flatind_HF[flat_ind] = initial_HF_OPs_.size() - 1;
										assert(Super_to_flatind_HF[flat_ind] >= 0);
										file_out_OP<<row<<"             "<<col<<"       "<<initial_HF_OPs_[Super_to_flatind_HF[flat_ind]].real()
														<<"             "<<initial_HF_OPs_[Super_to_flatind_HF[flat_ind]].imag()<<endl;
									}
								}
							}
						}

					}
				}
			}
		}


	}
	else{
		throw std::logic_error("Reading HF OPs from file not yet implemented");
	}

}

void MFParam_TL::initialPairingOrderParams(){
	initial_Pair_OPs_.clear();
	row_vector_Pair.clear();		col_vector_Pair.clear();

	double eps = 0.01;

	if(Parameters_TL_.read_OP_==false){
		if(Parameters_TL_.doPairing==true){
			string file_OP_="initial_pairing_OPs_file.txt";
			ofstream file_out_OP(file_OP_.c_str());
			file_out_OP<<"#row_index    col_index   initial_OP_"<<endl;

			for(int ind1=0;ind1<nsites;ind1++){
				for(int ind2=ind1;ind2<nsites;ind2++){
					if(Neighbors_TL_.pairing_pairs[ind1][ind2]!=0.0){
						//generating onsite pairing OPs:
						if(ind2==ind1){

							if(Parameters_TL_.pairingType==0){//onsite s-wave pairing ansatz
								//only up-dn terms:
								int prow = getRowColIndex(ind1,0);
								int pcol = getRowColIndex(ind1,1) + Hsize;
								int pflat_ind = prow + BdGsize*pcol;
								complex<double> val_op(random1(),random1());

								initial_Pair_OPs_.push_back(val_op*eps);
								row_vector_Pair.push_back(prow);
								col_vector_Pair.push_back(pcol);
								Super_to_flatind_Pair[pflat_ind] = initial_Pair_OPs_.size() - 1;
								assert(Super_to_flatind_Pair[pflat_ind] >= 0);
								file_out_OP<<prow<<"            "<<pcol<<"      "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].real()
												<<"             "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].imag()<<endl;
							}

						}
						else{
							//Note: for NN up-dn OPs we are interested particularly in identifying the NN sites.
							//The <c_{i,up}c_{j,dn}> pair will be modified into something like 0.5*<c_{i,up}c_{j,dn} +/- c_{i,dn}c_{j,up}>
							if(Parameters_TL_.pairingType==1){//offsite chiral p-wave ansatz
								for(int spin=0;spin<nspin;spin++){
									int prow = getRowColIndex(ind1,spin);
								
									for(int spin_p=spin;spin_p<nspin;spin_p++){
										int pcol = getRowColIndex(ind2,spin_p) + Hsize;
									
										if(spin_p==spin){//up-up and dn-dn terms:
											int pflat_ind = prow + BdGsize*pcol;
											complex<double> val_op((2*random1()-1),(2*random1()-1));
											initial_Pair_OPs_.push_back(val_op*eps);

											row_vector_Pair.push_back(prow);
											col_vector_Pair.push_back(pcol);
											Super_to_flatind_Pair[pflat_ind] = initial_Pair_OPs_.size() - 1;
											assert(Super_to_flatind_Pair[pflat_ind] >= 0);
											file_out_OP<<prow<<"            "<<pcol<<"      "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].real()
															<<"             "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].imag()<<endl;
										}
										else{//up-dn terms:
											int pflat_ind = prow + BdGsize*pcol;
											complex<double> val_op((2*random1()-1),(2*random1()-1));
											initial_Pair_OPs_.push_back(val_op*eps);
										
											row_vector_Pair.push_back(prow);
											col_vector_Pair.push_back(pcol);
											Super_to_flatind_Pair[pflat_ind] = initial_Pair_OPs_.size() - 1;
											assert(Super_to_flatind_Pair[pflat_ind] >= 0);
											file_out_OP<<prow<<"            "<<pcol<<"      "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].real()
															<<"             "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].imag()<<endl;
										}
									}
								}
							}
							if(Parameters_TL_.pairingType==2){//offsite chiral d-wave ansatz
								//only up-down terms:
								int prow = getRowColIndex(ind1,0);
								int pcol = getRowColIndex(ind2,1) + Hsize;
								int pflat_ind = prow + BdGsize*pcol;
								complex<double> val_op((2*random1()-1),(2*random1()-1));
							
								initial_Pair_OPs_.push_back(val_op);
								row_vector_Pair.push_back(prow);
								col_vector_Pair.push_back(pcol);
								Super_to_flatind_Pair[pflat_ind] = initial_Pair_OPs_.size() - 1;
								assert(Super_to_flatind_Pair[pflat_ind] >= 0);
								file_out_OP<<prow<<"            "<<pcol<<"      "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].real()
												<<"             "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].imag()<<endl;
							}
						
						}

					}
				}
			}
		}
	}
	else{
		throw std::logic_error("Reading pairing OPs from file not implemented");
	}

}


	/*For generic code purposes
	if(Parameters_TL_.singleImpurity==false){//onsite Pairing: generic
	
	int prow = getRowColIndex(ind1,0);
	int pcol = getRowColIndex(ind1,1) + Hsize;
	int pflat_ind = prow + BdGsize*pcol;
	complex<double> val_op(random1(), random1());
	initial_Pair_OPs_.push_back(val_op);
	row_vector_Pair.push_back(prow);
	col_vector_Pair.push_back(pcol);
	Super_to_flatind_Pair[pflat_ind] = initial_Pair_OPs_.size() - 1;
	assert(Super_to_flatind_Pair[pflat_ind] >= 0);
	file_out_OP<<prow<<"            "<<pcol<<"      "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].real()
					<<"             "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].imag()<<endl;
	}
	
	when ind1!=ind2
	if(Parameters_TL_.singleImpurity==false){//offsite Pairing: generic
	for(int spin=0;spin<nspin;spin++){
	int prow = getRowColIndex(ind1,spin);

	for(int spin_p=spin;spin_p<nspin;spin_p++){
	int pcol = getRowColIndex(ind2,spin_p) + Hsize;
	
	if(spin_p==spin){
	if(ind2>ind1){
	int pflat_ind = prow + BdGsize*pcol;
	complex<double> val_op(random1(), random1());
	
	initial_Pair_OPs_.push_back(val_op);
	row_vector_Pair.push_back(prow);
	col_vector_Pair.push_back(pcol);
	Super_to_flatind_Pair[pflat_ind] = initial_Pair_OPs_.size() - 1;
	assert(Super_to_flatind_Pair[pflat_ind] >= 0);
	file_out_OP<<prow<<"            "<<pcol<<"      "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].real()
					<<"             "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].imag()<<endl;
	}
	}
	else{
	int pflat_ind = prow + BdGsize*pcol;
	complex<double> val_op(random1(), random1());
	
	initial_Pair_OPs_.push_back(val_op);
	row_vector_Pair.push_back(prow);
	col_vector_Pair.push_back(pcol);
	Super_to_flatind_Pair[pflat_ind] = initial_Pair_OPs_.size() - 1;
	assert(Super_to_flatind_Pair[pflat_ind] >= 0);
	file_out_OP<<prow<<"            "<<pcol<<"      "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].real()
					<<"             "<<initial_Pair_OPs_[Super_to_flatind_Pair[pflat_ind]].imag()<<endl;
	}
	}
	}
	}
	*/


#endif
