///////////////////////////////////////////////////////////
//                                                       
//  RunCA.h                                             
//  Runoff_Model_Based_On_Cellular_Automata (RunCA)     
//  Erosion_Model_Based_On_Cellular_Automata (EroCA)    
//                                                      
//  Created by Shahla yavari and Qi Shao on 26/2/2021.  
//   All rights reserved.                               
//  Copyright © 2020 Qi Shao.
//  Copyright © 2021 Shahla yavari.
//                                                   
//////////////////////////////////////////////////////////


#ifndef CELL_H
#define CELL_H

#include <vector>

using namespace std;

class Cell {
  public:
    Cell();
    Cell(int x, int y, double elevation, bool valid);
    ~Cell();

    void SetLandElevation(double elevation) {land_elevation_ = elevation;};
    double GetLandElevation() {return land_elevation_;};

    void SetWaterSurfaceElevation(double height) {water_surface_elevation_ = height;};
    double GetWaterSurfaceElevation() {return water_surface_elevation_;};   

    void SetWaterDepth(double depth) {water_depth_ = depth;};
    double GetWaterDepth() {return water_depth_;};    

    void SetSoilType(int type) {soil_type_ = type;};
    int GetSoilType() {return soil_type_;};      

    void SetStoragePotential(double ws) {water_storage_potential_ = ws;};
    double GetStoragePotential() {return water_storage_potential_;};    

    void SetValid(bool valid) {valid_ = valid;};
    bool GetValid() {return valid_;};

    void SetX(int x) {x_ = x;};
    int GetX() {return x_;};  
    void SetY(int y) {y_ = y;};
    int GetY() {return y_;}; 

    void SetInfiltrationRate(double rate) {infiltration_rate_ = rate;};
    double GetInfiltrationRate() {return infiltration_rate_;};

    void SetInfiltrationCapacity(double capacity) {infiltration_capacity_ = capacity;};
    double GetInfiltrationCapacity() {return infiltration_capacity_;};

    void SetDrainageRate(double rate) {drainage_rate_ = rate;};
    double GetDrainageRate() {return drainage_rate_;};       
     
    void SetFlows(vector<double> flows);
    vector<double> GetFlows() {return flows_;}; 

	//*****************************************************************
	
	void SetSuspendedSedimentVolumes(vector<double> SusVol) { SusVol_i_ = SusVol; }; //i
	vector<double> GetSuspendedSedimentVolumes() { return SusVol_i_; }; //i	

	void SetBedMaterialVolumes(vector<double> BmVol) { BmVol_i_ = BmVol; }; //i
	vector<double> GetBedMaterialVolumes() { return BmVol_i_; }; //i	

	void SetErodedMaterialVolumes(vector<double> EroVol) { EroVol_i_ = EroVol; }; //i
	vector<double> GetErodedMaterialVolumes() { return EroVol_i_; }; //i

	void SetLandElavationChange(double change) {land_elevation_change_ = change;};
    double GetLandElavationChange() {return land_elevation_change_;};	

	void SetCumulativeLandElavationChange(double change) {cumulative_land_elevation_change_ = change;};
    double GetCumulativLandElavationChange() {return cumulative_land_elevation_change_;};		

	void SetSsus_i_flow(vector < vector<double> > Ssus_i_flow) { Ssus_i_flow_ = Ssus_i_flow; };
	vector<vector<double>>GetSsus_i_flow() {return Ssus_i_flow_;};

	void SetSbm_i_flow(vector < vector<double> > Sbm_i_flow) { Sbm_i_flow_ = Sbm_i_flow; };
	vector<vector<double>> GetSbm_i_flow() { return Sbm_i_flow_; };

	void SetSero_i_flow(vector < vector<double> > Sero_i_flow) { Sero_i_flow_ = Sero_i_flow; };
	vector<vector<double>> GetSero_i_flow() { return Sero_i_flow_; };
//**test
	void Setsero_i(double Sero){ Sero_i_=Sero;};
	double Getsero_i(){ return Sero_i_;};

	void Setsbm_i(double sbm_i){ sbm_i_=sbm_i;};
	double Getsbm_i() { return sbm_i_; };
	///

	void Setaccum_totsed(double com_accum_totsed) { accum_totsed_ = com_accum_totsed; };
	double Getaccum_totsed() { return accum_totsed_; };

    /*
	void Setaccum_sus_sedi(vector<double> com_accum_sus_sedi) { accum_sus_sedi_ = com_accum_sus_sedi; };
	vector<double> Getaccum_sus_sedi() {return  accum_sus_sedi_;};

	void Setaccum_bm_sedi(vector<double> com_accum_bm_sedi) { accum_bm_sedi_ = com_accum_bm_sedi; };
	vector<double> Getaccum_bm_sedi() { return  accum_bm_sedi_;};

	void Setaccum_ero_sedi(vector<double> com_accum_ero_sedi) { accum_ero_sedi_ = com_accum_ero_sedi; };
	vector<double> Getaccum_ero_sedi() { return accum_ero_sedi_;};
	*/
    
	//Qi added
	void SetSusSediInflow(vector<double> sus_sedi_inflow) { sus_sedi_inflow_ = sus_sedi_inflow; };
	vector<double> GetSusSediInflow() {return sus_sedi_inflow_;};

	void SetBmSediInflow(vector<double> bm_sedi_inflow) { bm_sedi_inflow_ = bm_sedi_inflow; };
	vector<double> GetBmSediInflow() { return  bm_sedi_inflow_;};

	void SetEroSediInflow(vector<double> ero_sedi_inflow) { ero_sedi_inflow_ = ero_sedi_inflow; };
	vector<double> GetEroSediInflow() { return ero_sedi_inflow_;};


	void SetSusSediOutflow(vector<double> sus_sedi_outflow) { sus_sedi_outflow_ = sus_sedi_outflow; };
	vector<double> GetSusSediOutflow() {return sus_sedi_outflow_;};

	void SetBmSediOutflow(vector<double> bm_sedi_outflow) { bm_sedi_outflow_ = bm_sedi_outflow; };
	vector<double> GetBmSediOutflow() { return  bm_sedi_outflow_;};

	void SetEroSediOutflow(vector<double> ero_sedi_outflow) { ero_sedi_outflow_ = ero_sedi_outflow; };
	vector<double> GetEroSediOutflow() { return ero_sedi_outflow_;};


		
	void SetFlowVelocity(vector<double> fv) { flow_velocities_ = fv; };
	vector<double> GetFlowVelocity() { return flow_velocities_; };
	
	void SetSlope(vector<double> s) { slopes_ = s; };
	vector<double> GetSlope() { return slopes_; };
	
	void SettotalSedimentTransport(vector<double> Stot) { tot_sediment_transport_ = Stot; };
	vector<double> GettotalSedimentTransport() { return tot_sediment_transport_; };

 	void SetBackupStates(vector<double> states) { backup_states_ = states; };
	vector<double> GetBackupStates() { return backup_states_; };   


    // bedload transportation parameters
	

	void SetBedMaterialVolumesb(vector<double> BmVolb) { BmVol_i_b = BmVolb; }; //i            1b
	vector<double> GetBedMaterialVolumesb() { return BmVol_i_b; }; //i	

	void SetErodedMaterialVolumesb(vector<double> EroVolb) { EroVol_i_b = EroVolb; }; //i           2b
	vector<double> GetErodedMaterialVolumesb() { return EroVol_i_b; }; //i

	void SetLandElavationChangeb(double changeb) {land_elevation_change_b = changeb;};
    double GetLandElavationChangeb() {return land_elevation_change_b;};	

	void SetCumulativeLandElavationChangeb(double changeb) {cumulative_land_elevation_change_b = changeb;};
    double GetCumulativLandElavationChangeb() {return cumulative_land_elevation_change_b;};		

   void SetCumulativeLandElavationChangesus(double changes) {cumulative_land_elevation_change_s = changes;};
    double GetCumulativLandElavationChangesus() {return cumulative_land_elevation_change_s;};	
	
	void SetSbm_i_flowb(vector < vector<double> > Sbm_i_flowb) { Sbm_i_flow_b = Sbm_i_flowb; };                //3b
	vector<vector<double>> GetSbm_i_flowb() { return Sbm_i_flow_b; };

	void SetSero_i_flowb(vector < vector<double> > Sero_i_flowb) { Sero_i_flow_b = Sero_i_flowb; };            //4b
	vector<vector<double>> GetSero_i_flowb() { return Sero_i_flow_b; };

	void Setsero_ib(double Serob){ Sero_i_b=Serob;}; //5b
	double Getsero_ib(){ return Sero_i_b;};

	void Setsbm_ib(double sbm_ib){ sbm_i_=sbm_ib;};   //6b
	double Getsbm_ib() { return sbm_i_b; };


	void Setaccum_totsedb(double com_accum_totsedb) { accum_totsed_b = com_accum_totsedb; };       //7b
	double Getaccum_totsedb() { return accum_totsed_b; };


	void SetBmSediInflowb(vector<double> bm_sedi_inflowb) { bm_sedi_inflow_ = bm_sedi_inflowb; };         //8b
	vector<double> GetBmSediInflowb() { return  bm_sedi_inflow_b;};

	void SetEroSediInflowb(vector<double> ero_sedi_inflowb) { ero_sedi_inflow_b = ero_sedi_inflowb; };           //9b
	vector<double> GetEroSediInflowb() { return ero_sedi_inflow_b;};


	void SetBmSediOutflowb(vector<double> bm_sedi_outflowb) { bm_sedi_outflow_b = bm_sedi_outflowb; };         //10b
	vector<double> GetBmSediOutflowb() { return  bm_sedi_outflow_b;};

	void SetEroSediOutflowb(vector<double> ero_sedi_outflowb) { ero_sedi_outflow_b = ero_sedi_outflowb; };        //11b
	vector<double> GetEroSediOutflowb() { return ero_sedi_outflow_b;};

	
	void SettotalSedimentTransportb(vector<double> Stotb) { tot_sediment_transport_b = Stotb; };          //12b
	vector<double> GettotalSedimentTransportb() { return tot_sediment_transport_b; };	


  private:
    int x_, y_;
    double land_elevation_; //m
    double water_surface_elevation_; //mm
    double water_depth_; //mm
    int soil_type_;
    bool valid_;
    double infiltration_rate_; //mm/h
    double infiltration_capacity_; //mm/h
    double drainage_rate_; //mm/h
    double water_storage_potential_; //mm
	
	
									 
	//****************************************************************************************
	
	vector<double> SusVol_i_; //volume (m3) of suspended sediment in differnt size fractions - 
	vector<double> BmVol_i_; //volume (m3) of bed material in differnt size fractions - 
	vector<double> EroVol_i_; //volume (m3) of eroded material in differnt size fractions - 
	double land_elevation_change_;
	double cumulative_land_elevation_change_;
	double Sero_i_;
	double sbm_i_;
	


	



    // ul  u  ur       0  1  2
    // l   c  r   <=>  3  4  5
    // dl  d  dr       6  7  8
    vector<double> flows_; //flow to neighbour cells (mm)
	
	vector<double> flow_velocities_;
	vector<double> slopes_;
	//**********************
    
	//total sediment tranpsorted to neighbour cells ([9])
	vector<double> tot_sediment_transport_; //stot (m3)
	//suspended sediment transported to neighbour cells from different fractions ([4][9])
	vector<vector<double>> Ssus_i_flow_; 
	//bed material/deposited sediment transported to neighbour cells from different fractions ([4][9])
	vector<vector<double>> Sbm_i_flow_; 
	//parent/eroded mateiral transported to neighbour cells from different fractions ([4][9])
    vector<vector<double>> Sero_i_flow_; 
    
	

    //incoming sediment amount at each fraction size from different sources 
	vector<double> sus_sedi_inflow_, bm_sedi_inflow_, ero_sedi_inflow_; //[4]
	vector<double> sus_sedi_outflow_, bm_sedi_outflow_, ero_sedi_outflow_; //[4]

	double accum_totsed_;

	vector<double> backup_states_; //added


		// Bed load parameter in private
	
	vector<double> BmVol_i_b; //volume (m3) of bed material in differnt size fractions - 
	vector<double> EroVol_i_b; //volume (m3) of eroded material in differnt size fractions - 
	double land_elevation_change_b;
	double cumulative_land_elevation_change_b;
	double Sero_i_b;
	double sbm_i_b;
	double land_elevation_change_s;
	double cumulative_land_elevation_change_s;
	
 
	//vector<double> slopes_;
	//**********************
    
	//total sediment tranpsorted to neighbour cells ([9])
	vector<double> tot_sediment_transport_b; //stot (m3)
	
	//bed material/deposited sediment transported to neighbour cells from different fractions ([4][9])
	vector<vector<double>> Sbm_i_flow_b; 
	//parent/eroded mateiral transported to neighbour cells from different fractions ([4][9])
    vector<vector<double>> Sero_i_flow_b; 
    	
    //incoming sediment amount at each fraction size from different sources 
	vector<double>  bm_sedi_inflow_b, ero_sedi_inflow_b; //[4]
	vector<double>  bm_sedi_outflow_b, ero_sedi_outflow_b; //[4]

	double accum_totsed_b;
};

#endif
