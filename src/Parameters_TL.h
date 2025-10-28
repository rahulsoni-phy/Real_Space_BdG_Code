#ifndef Parameters_TL_class
#define Parameters_TL_class
#include <iostream>
#include <fstream>
#include "tensors.h"


//this class will read from the **input file**
class Parameters_TL{

public:
	// Declaring system params:---->
	bool PBC_X, PBC_Y;
	int Lx_, Ly_, N_spin_=2;
	int Total_Sites, Hsize, BdG_Size;

	// Declaring model params:----->
	double t1_hop, t2_hop, t3_hop, t4_hop, t5_hop;
	double mu_fixed, temp_, beta;
	double Bfield, onsiteE;

	double U0_, U1_, U1byU0_factor;
	bool useHartree, useFock, doPairing;
	int pairingType;
	double pair_pot;

	//Disorder params:---->
	bool singleImpurity,chargeDisorder;
	double impurityval;
	int disorder_seed;

	// Observables details:--->
	bool get_Akw, get_Akyw;
	double w_min, w_max, dw_, eta_;

	//Declaring self-consistency params:---->
	bool Simple_mixing, Broyden_mixing;
	double Conv_err, alpha_,beta_damp;
	int Seed_, Max_iters;

	bool read_OP_;
	string input_OP_file;
	//----------------------->
	int npthreads_;

	void Initialize(std::string input_file);
	double matchstring(std::string file, std::string match);
	std::string matchstring2(std::string file, std::string match);
};

void Parameters_TL::Initialize(string input_file){

	cout << "-------------------------------------------\n"
		 << "Reading the input file = " << input_file << "\n"
		 << "-------------------------------------------\n"
		 << endl;
	//------------------------------------------------------------------------//
	// Reading and setting boundary conditions
	string PBC_X_string, PBC_Y_string;
	string pbc_x_out, pbc_y_out;

	PBC_X_string = matchstring2(input_file, "PBC_X");
	if (PBC_X_string == "True")
	{	PBC_X = true;		pbc_x_out = "PBC";		}
	else {	PBC_X = false;		pbc_x_out = "OBC";		}

	PBC_Y_string = matchstring2(input_file, "PBC_Y");
	if (PBC_Y_string == "True")
	{	PBC_Y = true;		pbc_y_out = "PBC";		}
	else {	PBC_Y = false;		pbc_y_out = "OBC";		}


	Lx_ = int(matchstring(input_file, "Sites_X"));
	Ly_ = int(matchstring(input_file, "Sites_Y"));

	Total_Sites = Lx_ * Ly_;
	Hsize = N_spin_*Total_Sites;
	BdG_Size = 2 * Hsize;

	cout << "Boundary conditions = " << pbc_x_out << "x" << pbc_y_out << "\n"
		 << "System size = " << Lx_<<"x"<<Ly_<<"\n"
		 << "Total size of the Hamiltonian = " << BdG_Size << "x" << BdG_Size << endl;
	//-----------------------------------------------------------------------//

	t1_hop = matchstring(input_file, "NN_Hopping");
	t1_hop = t1_hop * 0.001;
	t2_hop = -0.3881 * t1_hop;
	t3_hop = 0.1444 * t1_hop;
	t4_hop = -0.0228 * t1_hop;
	t5_hop = -0.0318 * t1_hop;

	U0_ =  matchstring(input_file, "Onsite_Int_U");
	U0_ = U0_ * 0.001;
	U1byU0_factor = matchstring(input_file, "NN_Int_factor");
	U1_ = U0_*U1byU0_factor;
	cout<<"Onsite interaction: U0 = "<<U0_<<", NN interaction: U1 = "<<U1_<<endl;

	onsiteE = matchstring(input_file, "OnsiteE");
	Bfield = matchstring(input_file, "MagneticField");

	mu_fixed = matchstring(input_file, "mu_fixed");
	mu_fixed = mu_fixed * 0.001;

	temp_ = matchstring(input_file, "Temperature");
	beta = 1.0 / (1.0 * temp_);

	string hartree_,fock_;
	if(abs(U0_)>1e-5){
		hartree_ = matchstring2(input_file,"Enable_Hartree");
		if(hartree_=="True"){	useHartree=true;
			cout<<"Hartree is turned ON"<<endl;
		}
		else{	useHartree=false;
			cout<<"Hartree is turned OFF"<<endl;
		}

		fock_ = matchstring2(input_file,"Enable_Fock");
		if(fock_=="True"){	useFock=true;
			cout<<"Fock is turned ON"<<endl;
		}
		else{	useFock=false;
			cout<<"Fock is turned OFF"<<endl;
		}
	}
	else{
		useHartree=false;	useFock=false;
		cout<<"Both Hartree and Fock are DISABLED"<<endl;
	}

	string dopairing_;
	dopairing_ = matchstring2(input_file,"Enable_Pairing");
	if(dopairing_=="True"){		
		doPairing=true;
		cout<<"Pairing is enabled"<<endl;

		pairingType = int(matchstring(input_file, "PairingType"));
		if(pairingType==0){
			cout<<"NN singlet pairing is considered"<<endl;
		}
		else{
			cout<<"NN triplet pairing is considered"<<endl;
		}
		pair_pot = matchstring(input_file, "Phenom_Pair_Pot");
		pair_pot = pair_pot * 0.001;
	}
	else{
		doPairing=false;
		cout<<"Pairing is disabled"<<endl;
	}

	//----------------------------------------------------------------------//
	string impurity_str, type_str;
	impurity_str = matchstring2(input_file,"SingleImpurity");
	if(impurity_str == "True"){
		singleImpurity = true;
		impurityval = matchstring(input_file, "SingleImpurityValue");
		impurityval = impurityval * 0.001;
		cout<<"Single-impurity check! impurity value ="<<impurityval<<endl;
	}
	else{
		disorder_seed = int(matchstring(input_file, "DisorderSeed"));
	}

	type_str = matchstring2(input_file,"DisorderType");
	if(type_str == "Charge"){
		chargeDisorder = true;
	}
	else{
		chargeDisorder =false;
	}
	
	//----------------------------------------------------------------------//
	string mixing_, mixing_out;
	mixing_ = matchstring2(input_file,"Simple_Mixing");
	if(mixing_=="True"){
		Simple_mixing=true;		Broyden_mixing=false;
		mixing_out="Simple_Mixing";
	}
	else{
		Simple_mixing=false;    Broyden_mixing=true;
		mixing_out="Broyden_Mixing";
	}

	Conv_err = matchstring(input_file, "Convergence_Error");
	alpha_ = matchstring(input_file, "alpha_OP");
	beta_damp = matchstring(input_file, "beta_damp");

	Seed_ = int(matchstring(input_file, "Random_Seed"));
	Max_iters = int(matchstring(input_file, "Max_Iterations"));
	cout<<"Self-consistency solver = "<<mixing_out<<endl;
	//----------------------------------------------------------------------//
	
	string qpi_str, Akw_str, Akyw_str;
	/*qpi_str = matchstring2(input_file, "Calculate_QPI");
	if (qpi_str == "True"){
		get_qpi = true;
		cout << "Measuring QPI" << endl;		
	}
	else {	get_qpi = false;	} */

	Akw_str = matchstring2(input_file, "Spectral_Akw");
	if (Akw_str == "True"){
		get_Akw = true;
		cout << "Measuing spectral function along path A(kx,ky;w)" << endl;
	}
	else {	get_Akw = false;	}

	Akyw_str = matchstring2(input_file, "Spectral_Akyw");
	if (Akyw_str == "True"){
		get_Akyw = true;
		cout << "Measuing spectral function A(ky;w)" << endl;
	}
	else {	get_Akyw = false;	}

	w_min = matchstring(input_file, "omega_min");
	w_max = matchstring(input_file, "omega_max");
	dw_ = matchstring(input_file, "d_omega");
	eta_ = matchstring(input_file, "broadening");
	//----------------------------------------------------------------------//
	string read_op_str;
	read_op_str = matchstring2(input_file, "Read_OP_from_file");
	if (read_op_str == "True"){
		read_OP_ = true;
		input_OP_file = matchstring2(input_file, "OP_Filename");
		cout<<"Reading order-params from file:"<<input_OP_file<<endl;
	}
	else{
		read_OP_ = false;
	}
	//----------------------------------------------------------------------//

	npthreads_ = int(matchstring(input_file, "Threads"));
	if (npthreads_ > 1)
	{
		cout << "Total threads= " << npthreads_ << endl;
	}

	cout << "-------------------------------------------\n"
		 << "Finish Reading the input file" << "\n"
		 << "-------------------------------------------\n"
		 << endl;
}

double Parameters_TL::matchstring(string file, string match){

	std::string test, line;
	std::ifstream inputfile(file);
	double amount;

	bool pass = false;

	while (std::getline(inputfile, line)){
		std::istringstream iss(line);

		if (std::getline(iss, test, '=') && !pass){
			if (iss >> amount && test == match){
				pass = true;
			}
			if (pass){
				break;
			}
		}
	}

	if (!pass){
		throw std::invalid_argument("Missing argument in the input file: " + match);
	}

	return amount;
}

string Parameters_TL::matchstring2(string file, string match){
	std::string line;
    std::ifstream inputfile(file);
    std::string amount;
    int offset;

	/* Referenced from https://cplusplus.com/forum/beginner/121556/
			 * https://stackoverflow.com/questions/12463750/ */
	if(inputfile.is_open()) {
        while(std::getline(inputfile, line)) {
            if((offset = line.find(match, 0)) != std::string::npos) {
                amount = line.substr(offset + match.length() + 1);
                break; // Break early if the match is found
            }
        }
        inputfile.close();
    } else {
        std::cerr << "Unable to open input file while in the Parameter class." << std::endl;
        return "";
    }

//    std::cout << match << " = " << amount << std::endl;
    return amount;
}

#endif
