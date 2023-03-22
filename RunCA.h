///////////////////////////////////////////////////////
//                                                   //
//  RunCA.h                                          //
//  Runoff_Model_Based_On_Cellular_Automata (RunCA)  //
//  Erosion_Model_Based_On_Cellular_Automata (EroCA)  //
//                                                   //
//  Created by Shahla yavari and Qi Shao on 26/2/2021. //
//  Copyright © 2020 Qi Shao. All rights reserved.   //
//                                                   //
///////////////////////////////////////////////////////


#ifndef RUNCA_H
#define RUNCA_H

#include "Cell.h"
#include <vector>

using namespace std;

int main();

class RunCA {
  public:
    RunCA( const char* model_name );
    ~RunCA();

    int Run();
    void Initialisation();
    void TimeIntegration();
    
    //reading input files
    void ReadElevations();
    void SetupBordersAndOutlet();
    void ReadSoilTypes();
    void ReadInfiltrationParameters();
    void ReadRainfallRates();
  	
    void ReadErosionParameters(); //Qi: to be consistent with names of other functions    

    //computations
    double ComputeTimeToPonding(Cell* cell);
    void UpdateStates(Cell* cell, double& sum_r, double& sum_i, double sum_d); 
  	void UpdateErosionStates(Cell* cell);
	  double ComputeFlows(Cell* cell);
    double ComputeFlows2(Cell* cell); //added
  	double ComputeSedimentTransport(Cell* cell);
    void AccumulateFlows(Cell* cell);
  	void AccumulateSedimentTransport(Cell*cell);
    void AccumulateBedloadTransport(Cell* cell); //
    double ComputeBedloadTransport(Cell* cell);  //
    void  UpdateErosionStatesByBedload(Cell* cell);  //

    double ComputeTimeStepSize(Cell *cell);
    double ComputeGlobalTimeStepSize(bool initial_guess);
    void ComputeFlowVelocities(Cell* cell); 
    void ComputeSlopes(Cell*cell); 

    
    //outputs
    void SaveInputDataToVti(); //for visualisation in ParaView
    void SaveResultsToVti();
    void SaveWaterDepthToText(); //for visualisation in ArcGIS
    void WriteOutputFileHeaders();
    void WriteResults(); //for plotting in Excel
    void WriteSummary();
    void SaveElevationUpdateToText(); //  for visualisation change elevation in ArcGIS
    void SaveTotalAccumulationSedimentVolumeToText();
    
    void SaveCumulativeLandElavationChangeToText(); //test for effect of change of elevation on infiltration
    


  private:
    const char* model_name_;
    bool debug_;
    int nx_, ny_, effective_cells_;
    double cell_size_;
    double duration_;
    double time_increment_;
    double adjust_coefficient_ = 1.;
    double model_time_ = 0.;
    int iteration_ = 0;
    int threads_;
    double tot_runoff_=0., tot_rain_=0., tot_infiltration_=0., tot_drainage_=0., tot_ponding_=0.;
    double run_time_;

    vector< vector<Cell*> > cells_;    
    vector<double> rainfall_;

    //infiltration parameters vector for Holtan equation
    //[0] soil type
    //[1] If (steady infiltration rate)
    //[2] I0 (initial infiltration rate)
    //[3] P (infiltration decay coefficient)
    //[4] TP (total porosity)
    //[5] WC0 (initial soil water content)
    //[6] FC (filed capacity)
    //[7] D (control zone depth: mm)
    //[8] N (manning's coefficent)
    vector< vector<double> > iparams_; 

  	//vector<double>  BM_i , SusVol_i ; //: these should be member variables of Cell?
	  //double BMtot = 0., SusVol = 0.; //: these should be member variable of Cell?
    
    //erosion parameters vector
    //: why do we need the input percentages of size fractions?? 
	  //[0] soil type
	  //[1] clay
	  //[2] silt  
	  //[3] sand
	  //[4] gravel
    //: fall velocities should not be included here
    /* 
	  //[5] wc  (fall velocity clay)
	  //[6] ws  (fall velocity silt)
	  //[7] wsa (fall velocity sand)
	  //[8] wg  (fall velocity gravel)
	  //[9] k   (erodibility factor in the USLE equation)
	  //[10] c  (cropping management factor)
	  //[11] p  (conservation practice factor) 
    */
    //
  	//[5] k (erodibility factor in the USLE equation)
	  //[6] c (cropping management factor)
	  //[7] p (conservation practice factor)    
	vector< vector<double> > erosion_para_;
  double tot_transportedsediment_=0.;
   bool with_bedload_transport_ = false;
		
};

#endif
