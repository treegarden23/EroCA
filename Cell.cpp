///////////////////////////////////////////////////////
//                                                   //
//  Cell.cpp                                         //
//  Runoff_Model_Based_On_Cellular_Automata (RunCA)  //
//                                                   //
//  Created by Qi Shao on 26/2/2020.                 //
//  Copyright © 2020 Qi Shao. All rights reserved.   //
//                                                   //
///////////////////////////////////////////////////////



#include "Cell.h"
#include <iostream>

using namespace std;


Cell::Cell()
{

}


Cell::Cell(int x, int y, double elevation, bool valid) 
	: x_(x),
		y_(y),
		land_elevation_(elevation), //m
		water_surface_elevation_(elevation*1000.), //mm
		water_depth_(0.), //mm
		valid_(valid),
		flows_(9, 0.), //mm
		tot_sediment_transport_(9, 0.),
		land_elevation_change_ (0.),
		cumulative_land_elevation_change_ (0.),
		tot_sediment_transport_b(9, 0.)
		//*********
{
	//initialise suspended sediment, bed material and eroded material volumes to zero - Qi
	
	vector<double> zeros(4, 0.);
	SusVol_i_ = zeros;
	BmVol_i_ = zeros;
	EroVol_i_ = zeros;
	vector<vector<double>> zeros_zeros (4, vector<double>(9,0.)); 
	Ssus_i_flow_ = zeros_zeros;
    Sbm_i_flow_ = zeros_zeros;
	Sero_i_flow_ = zeros_zeros;
	/*
	accum_sus_sedi_ = zeros;
	accum_bm_sedi_ = zeros;
	accum_ero_sedi_ = zeros;
	*/
    //Qi added
	sus_sedi_inflow_ = zeros;
	bm_sedi_inflow_ = zeros;
	ero_sedi_inflow_ = zeros;
	sus_sedi_outflow_ = zeros;
	bm_sedi_outflow_ = zeros; 
	ero_sedi_outflow_ = zeros;	


		//bedload
	
	
	BmVol_i_b = zeros;
	EroVol_i_b = zeros;
	
	
    Sbm_i_flow_b = zeros_zeros;
	Sero_i_flow_b = zeros_zeros;

	
	bm_sedi_inflow_b = zeros;
	ero_sedi_inflow_b = zeros;
	
	bm_sedi_outflow_b = zeros; 
	ero_sedi_outflow_b = zeros;
	

}




Cell::~Cell()
{
    
}


void Cell::SetFlows(vector<double> flows){

    if(flows.size() != 9){
        cerr<<"Error: in Cell::SetFlows: wrong input flow vector size = "<<flows.size()<<", should be 9";
        terminate();
    }

    flows_ = flows;
}
//****************************************
//void Cell::SettotalSedimentTransport(vector<double> stot) {
	//if (stot.size() != 9) {
		//cerr << "Error:in cell::SetSedimentTransport:wrong input sed vector size=" << stot.size() << ", should be 9";
		//terminate();
	//}
	//tot_sediment_transport_ = stot;
	

//}


