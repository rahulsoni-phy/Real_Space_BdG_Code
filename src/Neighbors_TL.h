#ifndef Neighbors_TL_class
#define Neighbors_TL_class
#include <algorithm>
#include "tensors.h"
#include "Parameters_TL.h"


//this class will find the neighbors based on distance and boundary conditions
class Neighbors_TL{

    public:
    //class constructor
    Neighbors_TL(Parameters_TL &Parameters_TL__): Parameters_TL_(Parameters_TL__){
        Initialize();
    }

    void Initialize();
    int getIndex(int rx, int ry);
    pair_int getCoordinates(int index);
    void generateNeighborOffsets(int max_order);
    Mat_1_int findNeighbor(int index, int order);
    void generateInteractionPairs();

    Parameters_TL &Parameters_TL_;
    Mat_1_int dist_to_neigh_order = {1, 3, 4, 7, 12};
    //dist_to_neigh_order[0]=1 -> corresponds to NN
    //dist_to_neigh_order[1]=3 -> corresponds to NNN 
    //dist_to_neigh_order[2]=4 -> corresponds to NNNN and so on
   
    int lx_,ly_,nsites;
    Mat_2_intpair neighbor_offsets_per_order;
    Mat_2_doub interaction_pairs;
    Mat_2_doub pairing_pairs;
};

void Neighbors_TL::Initialize(){
    lx_ = Parameters_TL_.Lx_;
    ly_ = Parameters_TL_.Ly_;
    nsites = Parameters_TL_.Total_Sites;

    generateNeighborOffsets(dist_to_neigh_order.size());
    interaction_pairs.resize(nsites);
    pairing_pairs.resize(nsites);
    for(int ind1=0; ind1<nsites; ind1++){
        interaction_pairs[ind1].resize(nsites);
        pairing_pairs[ind1].resize(nsites);
    }
}

int Neighbors_TL::getIndex(int rx, int ry){
    return rx + lx_*ry;
}

pair_int Neighbors_TL::getCoordinates(int index){
    if (index < 0 || index >= nsites) {
        throw out_of_range("index out of range in getCoordinates function!");
    }
    int rx = index % lx_;
    int ry = index / lx_;

    return {rx, ry}; 
}

void Neighbors_TL::generateNeighborOffsets(int max_order){
    
    if(max_order < 1 || max_order > dist_to_neigh_order.size()){
        throw std::invalid_argument("invalid max_order in generateNeighborOffsets function!");
    }

    neighbor_offsets_per_order.resize(max_order);

    int range = dist_to_neigh_order[max_order - 1];

    //distance between two points r1 = a1*v1 + b1*v2, r2 = a2*v1 + b2*v2 on triangular lattice will be
    //d = dist^2 = (a1-a2 + (b1-b2)/2)^2 + (3/4)(b1-b2)^2 = (a1-a2)^2 + (b1-b2)^2 + (a1-a2)*(b1-b2)
    //d = dist^2 = dx^2 + dy^2 + dx*dy

    //Identifying the corresponding neighbors based on their order and distance
    for(int dx = -range; dx <= range; dx++){
        for(int dy = -range; dy <= range; dy++){

            if(dx == 0 && dy == 0){
                continue;
            }

            int d = dx*dx + dy*dy + dx*dy;

            for(int order = 0; order < max_order; order++){
                if(d == dist_to_neigh_order[order]){
                    // Adding (dx, dy) pair to the corresponding neighbor order
                    neighbor_offsets_per_order[order].emplace_back(pair_int{dx,dy});
                    break; // Going for the next (dx, dy) pair once a match is found
                }
            }

        }
    }

    //Removing duplicated entries if any
    for(int order = 0; order < max_order; order++){
        Mat_1_intpair &offsets = neighbor_offsets_per_order[order];

        sort(offsets.begin(), offsets.end());
        auto new_end = unique(offsets.begin(), offsets.end());
        offsets.erase(new_end, offsets.end());
    }

}

Mat_1_int Neighbors_TL::findNeighbor(int index, int order){
    if (index < 0 || index >= nsites) {
        throw out_of_range("index out of range in findNeighbor function!");
    }

    pair_int coordinates = getCoordinates(index);
    int rx = coordinates.first;
    int ry = coordinates.second;
    
    //auto [rx,ry] = getCoordinates(index);
    Mat_1_int neighbors;
    Mat_1_intpair &offsets = neighbor_offsets_per_order[order];

    for(auto &offset : offsets){
        int nx = rx + offset.first;
        int ny = ry + offset.second;

        if(Parameters_TL_.PBC_X==true && Parameters_TL_.PBC_Y==true){
            if(nx < 0){nx += lx_;}
            if(nx >= lx_){nx -= lx_;}

            if(ny < 0){ny += ly_;}
            if(ny >= ly_){ny -= ly_;}
        }
        else if (Parameters_TL_.PBC_X==true && Parameters_TL_.PBC_Y==false){
            if(nx < 0){nx += lx_;}
            if(nx >= lx_){nx -= lx_;}
        }
        else if (Parameters_TL_.PBC_X==false && Parameters_TL_.PBC_Y==true){
            if(ny < 0){ny += ly_;}
            if(ny >= ly_){ny -= ly_;}
        }
        else {
            // No wrapping in x or y (both false)
        }
        

        if(nx >= 0 && nx < lx_ && ny >= 0 && ny < ly_){
            
            int neighbor_index = getIndex(nx, ny);
            neighbors.push_back(neighbor_index);
        }
    }

    return neighbors;
}

void Neighbors_TL::generateInteractionPairs(){
    for(int ind1=0; ind1<nsites; ind1++){
        for(int ind2=0; ind2<nsites; ind2++){
            interaction_pairs[ind1][ind2] = 0.0;
            pairing_pairs[ind1][ind2] = 0.0;
        }
    }

    double temp_val=1.0;
    for(int rx=0;rx<lx_;rx++){
        for(int ry=0;ry<ly_;ry++){
            int ind = getIndex(rx, ry);
            
            //HF interaction pairs:
            //onsite pairs:
            if(abs(Parameters_TL_.U0_)>1e-5){
                interaction_pairs[ind][ind] = 1.0*Parameters_TL_.U0_;
            }
            //NN pairs:
            if(abs(Parameters_TL_.U1_)>1e-5){
                Mat_1_int neighbors_ = findNeighbor(ind, 0);
                for(auto &neigh : neighbors_){
                    interaction_pairs[ind][neigh] = 1.0*Parameters_TL_.U1_;
                }
            }
           
            if(Parameters_TL_.doPairing==true){
                //Pairing interaction pairs:
                if(Parameters_TL_.pairingType==0){//onsite pairs
                    pairing_pairs[ind][ind] = 1.0*Parameters_TL_.pair_pot;
                }
                else{
                    Mat_1_int neighbors_ = findNeighbor(ind, 0);
                    for(auto &neigh : neighbors_){
                        pairing_pairs[ind][neigh] = 1.0*Parameters_TL_.pair_pot;
                    }
                }
            }
            
        }
    }
    //verification:
    for(int ind1=0;ind1<nsites;ind1++){
        for(int ind2=0;ind2<nsites;ind2++){
            assert(interaction_pairs[ind1][ind2] == interaction_pairs[ind2][ind1]);
            assert(pairing_pairs[ind1][ind2] == pairing_pairs[ind2][ind1]);
        }
    }
    
    string file_int_pairs = "interaction_pairs.txt";
    ofstream file_out(file_int_pairs.c_str());
    file_out<<"#site1    site2      values"<<endl;
    for(int ind1=0; ind1<nsites; ind1++){
        for(int ind2=ind1; ind2<nsites; ind2++){
            if(interaction_pairs[ind1][ind2] != 0.0){
                file_out<<ind1<<"       "<<ind2<<"      "<<interaction_pairs[ind1][ind2]<<endl;
            }
        }
    }
    string file_pairing_int_pairs = "pairing_interaction_pairs.txt";
    ofstream file_pairing_out(file_pairing_int_pairs.c_str());
    file_pairing_out<<"#site1    site2      values"<<endl;
    for(int ind1=0; ind1<nsites; ind1++){
        for(int ind2=ind1; ind2<nsites; ind2++){
            if(pairing_pairs[ind1][ind2] != 0.0){
                file_pairing_out<<ind1<<"       "<<ind2<<"      "<<pairing_pairs[ind1][ind2]<<endl;
            }
        }
    }
}

#endif
