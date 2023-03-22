
///////////////////////////////////////////////////////
//                                                   //
//  Cell.cpp                                         //
//  Runoff_Model_Based_On_Cellular_Automata (RunCA)  //
//                                                   //
//  
//EroCA_Model_Based_On_Cellular_Automata             //
//                                                   //
//  Created by Shahla yavari and  Qi Shao on 26/2/2021.
//                                                   //
//  Created by Qi Shao on 26/2/2021. 
//  Copyright © 2020 Qi Shao. All rights reserved.   //
//                                                   //
///////////////////////////////////////////////////////


#include "RunCA.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <sstream>
#include <math.h>
#include <algorithm> 
#include <iomanip>
#include <omp.h>

using namespace std;

int main()
{
    cout << "________________________________________________________________\n\n";
    cout << "Runoff Model Based On Cellular Automata (RunCA), vs. 2020 beta.\n";
    cout << "________________________________________________________________\n\n";

    string model_name("undefined");
    cout <<"\nEnter name of simulation model: ";
    cin >> model_name; 

    RunCA simulator(model_name.c_str());
    simulator.Run();
    return 0;
}


RunCA::RunCA(const char* model_name)
  : model_name_(model_name),
    debug_(true)
{

}


RunCA::~RunCA()
{
    delete model_name_;
}


int RunCA::Run()
{
    cout<<"\n*****Simulation Started*****"<<endl;
    Initialisation();
    TimeIntegration();
    cout<<"\n*****Simulation Completed*****"<<endl;
    return 0;
}


void RunCA::Initialisation()
{
    cout<<"\n------Start model initialisation------"<<endl;
    ReadElevations();
    SetupBordersAndOutlet();
    ReadInfiltrationParameters();
    ReadSoilTypes();
    ReadRainfallRates();
    ReadErosionParameters(); 
    SaveInputDataToVti();
    cout<<"------Model initialisation finished------"<<endl;
}


void RunCA::ReadElevations()
{
    cout<<"\nStart reading elevation file..."<<endl;
    string file_name(model_name_);
    file_name += "-input-elevations.txt";
    ifstream ifs(file_name.c_str());
    if( !ifs.is_open() ) {
        cerr<<"input elevation file ( " << file_name << " ) does not exist, please check"<<endl;
        terminate();
    } else{
        cout<<"  elevation file opened successfully"<<endl;
    }

    string text_line, token;
    //read in 1st line, should be something like: ncols	99
    getline(ifs,text_line);  
    istringstream iss (text_line); 
    iss >> token;
    iss >> token;
    nx_ = stoi(token);
    if(debug_)  cout<<"  number of cells in x direction (each row) = "<<nx_<<endl;
    
    //read in 2nd line, should be something like: nrows	86
    getline(ifs,text_line);
    while( text_line.length()==0 )  getline(ifs,text_line); //skip blank lines if there is any
    iss.str(text_line);
    iss >> token;
    iss >> token;
    ny_ = stoi(token);
    if(debug_)  cout<<"  number of cells in y direction (each column) = "<<ny_<<endl;

    cells_.resize(ny_+2); //include two border rows at top and bottom each
    for(int y=0; y<ny_+2; y++) cells_[y].resize(nx_+2); //include two side borders

    //skip 3rd and 4th lines
    getline(ifs,text_line);
    while( text_line.length()==0 )  getline(ifs,text_line); //skip blank lines if there is any   
    getline(ifs,text_line);
    while( text_line.length()==0 )  getline(ifs,text_line); //skip blank lines if there is any  

    //read in 5th line, should be something like: cellsize	15
    getline(ifs,text_line);
    while( text_line.length()==0 )  getline(ifs,text_line); //skip blank lines if there is any   
    iss.str(text_line);
    iss >> token;
    iss >> token;
    cell_size_ = stod(token);
    if(debug_)  cout<<"  cell size = "<<cell_size_<<" m"<<endl;
    double lx = nx_ * cell_size_ * 1.0;
    double ly = ny_ * cell_size_ * 1.0; 
    if(debug_)  cout<<"  model dimensions: "<<lx<<"m(X) x "<<ly<<"m(Y)"<<endl;

    //skip 6th line
    getline(ifs,text_line);
    while( text_line.length()==0 )  getline(ifs,text_line); //skip blank lines if there is any

    //start reading elevation values row by row
    int y_count(0), x_count(0), cell_count(0), valid_cells(0);
    double elevation;
    bool valid;
    getline(ifs,text_line);
    while( text_line.length()==0 )  getline(ifs,text_line); //skip blank lines if there is any
    do{
        y_count++;
        istringstream ss(text_line); 
        x_count=0;
        while(ss >> token){
            x_count++;
            elevation = stod(token);
            valid = (elevation == -9999) ? false : true;
            Cell* cell = new Cell(x_count, y_count, elevation, valid);
            if(valid) valid_cells++;
            cells_[y_count][x_count]=cell;
            cell_count++;
        }

        if(x_count != nx_) {
            cerr<<"  ERROR: in row "<<y_count<<" counted values("<<x_count<<") do not match specified number of columns ("<<nx_<<")"<<endl;
            cerr<<"  Please check your input elevation file"<<endl;
            terminate();        
        }

    } while(getline(ifs,text_line) && text_line.length()!=0);

    effective_cells_ = valid_cells;

    if(y_count != ny_) {
        cerr<<"  ERROR: counted rows("<<y_count<<") do not match specified number of rows ("<<ny_<<")"<<endl;
        cerr<<"  Please check your input elevation file"<<endl;
        terminate();        
    }

    if(cell_count != nx_ * ny_) {
        cerr<<"  ERROR: inserted number of cells("<<cell_count<<") do not match values in input files("<<nx_*ny_<<")"<<endl;
        cerr<<"  Please check your input elevation file"<<endl;
        terminate();        
    }    

    if(debug_) cout<<"  read in elevation value for "<<cell_count<<" cells, including "<<effective_cells_<<"  valid cells"<<endl;

    ifs.close();

    cout<<"Finish reading elevation file."<<endl; 
}


void RunCA::SetupBordersAndOutlet()
{
    cout<<"\nStart setting up borders and outlet..."<<endl;
    for (int y=0,x=0;x<nx_+2;x++){ //top border
        Cell* cell = new Cell(x,y,-9999,false);
        cells_[y][x] = cell;
    };
    for (int y=ny_+1,x=0;x<nx_+2;x++){ //bottom border
        Cell* cell = new Cell(x,y,-9999,false);
        cells_[y][x] = cell;
    };
    for (int x=0,y=0;y<ny_+2;y++){ //left border
        Cell* cell = new Cell(x,y,-9999,false);
        cells_[y][x] = cell;
    };
     for (int x=nx_+1,y=0;y<ny_+2;y++){ //right border
        Cell* cell = new Cell(x,y,-9999,false);
        cells_[y][x] = cell;
    };   
    int outlet_count(0);
    for (int y=ny_,x=1;x<nx_+1;x++) { //bottom outlet
        if (!cells_[y][x]->GetValid()) { //invalid cell
            continue;
        } else { //valid cell
            double h_above = cells_[y-1][x]->GetLandElevation();
            double h = cells_[y][x]->GetLandElevation();
            double h_below = h - (h_above - h);
            //double h_below = h - (h_above - h) * 0.1; //added
            //double h_below = h - (h_above - h) * 5.; //Qi added
            //double h_below = h - (h_above - h);
            //double h_below = h - (h_above - h) * 0.1;

            Cell* cell = new Cell(x, y+1, h_below, true);
            cells_[y+1][x] = cell;
            outlet_count++;
        };
    };  
    if(outlet_count == 0){
        cerr<<"  ERROR: there seems to be no outlet cells, please check your input elevation file"<<endl;
        terminate();
    }
    if(debug_)  cout<<"  found "<<outlet_count<<" outlet cells"<<endl;

    cout<<"Finish setting up borders and outlet."<<endl;
}




void RunCA::ReadSoilTypes()
{
    cout<<"\nStart reading soil types..."<<endl;
    string file_name(model_name_);
    file_name += "-input-soiltypes.txt";
    ifstream ifs2(file_name.c_str());
    if( !ifs2.is_open() ) {
        cerr<<"input soil type file ( " << file_name << " ) does not exist, please check"<<endl;
        terminate();
    } else{
        cout<<"  soil type file opened successfully"<<endl;
    }

    //read soil types file
    string text_line, token;
    //1st line, should be something like: ncols	99
    getline(ifs2,text_line);  
    istringstream iss2 (text_line); 
    iss2 >> token;
    iss2 >> token;
    if (stoi(token) != nx_) {
        cerr<<"  ERROR: ncols = "<<stoi(token)<<" does not match number of cells in x direction = "<<nx_<<endl;
        cerr<<"  Please check your input soil type file"<<endl;
        terminate();
    }
    //2nd line, should be something like: nrows	86
    getline(ifs2,text_line);
    while( text_line.length()==0 )  getline(ifs2,text_line); //skip blank lines if there is any
    iss2.str(text_line);
    iss2 >> token;
    iss2 >> token;
    if (stoi(token) != ny_) {
        cerr<<"  ERROR: nrows = "<<stoi(token)<<" does not match number of cells in y direction = "<<ny_<<endl;
        cerr<<"  Please check your input soil type file"<<endl;
        terminate();
    }
    //skip 3rd and 4th lines
    getline(ifs2,text_line);
    while( text_line.length()==0 )  getline(ifs2,text_line); //skip blank lines if there is any   
    getline(ifs2,text_line);
    while( text_line.length()==0 )  getline(ifs2,text_line); //skip blank lines if there is any  
    //read in 5th line, should be something like: cellsize	15
    getline(ifs2,text_line);
    while( text_line.length()==0 )  getline(ifs2,text_line); //skip blank lines if there is any   
    iss2.str(text_line);
    iss2 >> token;
    iss2 >> token;
    if (stod(token) != cell_size_) {
        cerr<<"  ERROR: cellsize = "<<stoi(token)<<" m does not match size of cell in the model = "<<cell_size_<<" m"<<endl;
        cerr<<"  Please check your input soil type file"<<endl;
        terminate();
    }
    //skip 6th line
    getline(ifs2,text_line);
    while( text_line.length()==0 )  getline(ifs2,text_line); //skip blank lines if there is any
    //start reading soil type values row by row
    int y_count(0), x_count(0), cell_count(0), valid_cells(0);
    int soiltype;
    getline(ifs2,text_line);
    while( text_line.length()==0 )  getline(ifs2,text_line); //skip blank lines if there is any
    do{
        y_count++;
        if(y_count>ny_){
            cerr<<"  ERROR: row number = "<<y_count<<" exceeds maximun row number =  "<<ny_<<endl;
            cerr<<"  Please check your input soil type file"<<endl;
            terminate();                  
        }
    
        istringstream ss(text_line); 
        x_count=0;
        while(ss >> token){
            x_count++;
            if(x_count>nx_){
                cerr<<"  ERROR: column number = "<<x_count<<" exceeds maximun column number =  "<<nx_<<endl;
                cerr<<"  Please check your input soil type file"<<endl;
                terminate();                  
            }            
            soiltype = stoi(token); 
            if(cells_[y_count][x_count]->GetValid()){
                if(soiltype==-9999){
                    cerr<<"  ERROR: nodata soiltype value is assigned for a valid cell located at row: "<<y_count<<" and column: "<<x_count<<endl;
                    cerr<<"  Please check your input soil type file"<<endl;
                    terminate();                  
                }
                if(soiltype<1 || soiltype>iparams_.size()){
                    cerr<<"  ERROR: input soiltype = "<<soiltype<<" is out of range between 0 and "<<iparams_.size()<<endl;
                    cerr<<"  Please check your input soil type file"<<endl;
                    terminate();         
                }
                cells_[y_count][x_count]->SetSoilType(soiltype);
                //initialise water storage potential
                vector<double> params = iparams_[soiltype-1];
                double porosity = params[4];
                double initial_water_content = params[5];
                double control_zone_depth = params[7];
                double water_storage_potential = (porosity - initial_water_content) * control_zone_depth; //mm
                water_storage_potential = max(water_storage_potential, 0.);
                cells_[y_count][x_count]->SetStoragePotential(water_storage_potential);
                cells_[y_count][x_count]->SetInfiltrationRate(params[2]); //initialised to I0
                cells_[y_count][x_count]->SetInfiltrationCapacity(params[2]); //initialised to I0
                cells_[y_count][x_count]->SetDrainageRate(0.); //initialised to 0



              
                valid_cells++;
            }
            cell_count++;
        }

        if(x_count != nx_) {
            cerr<<"  ERROR: in row "<<y_count<<" counted values("<<x_count<<") do not match specified number of columns ("<<nx_<<")"<<endl;
            cerr<<"  Please check your input soil type file"<<endl;
            terminate();        
        }

    } while(getline(ifs2,text_line) && text_line.length()!=0);

    if(y_count != ny_) {
        cerr<<"  ERROR: counted rows("<<y_count<<") do not match specified number of rows ("<<ny_<<")"<<endl;
        cerr<<"  Please check your input soil type file"<<endl;
        terminate();        
    }

    if(cell_count != nx_ * ny_) {
        cerr<<"  ERROR: assigned number of cells("<<cell_count<<") do not match values in input files("<<nx_*ny_<<")"<<endl;
        cerr<<"  Please check your input soil type file"<<endl;
        terminate();        
    }   

     if(valid_cells != effective_cells_) {
        cerr<<"  ERROR: assigned valid cells("<<valid_cells<<") do not match effective cells in the model("<<effective_cells_<<")"<<endl;
        cerr<<"  Please check your input soil type file"<<endl;
        terminate();        
    }       

    if(debug_) cout<<"  read in soil type value for "<<cell_count<<" cells, including "<<valid_cells<<" valid cells"<<endl;

    ifs2.close();
    cout<<"Finish reading soil types"<<endl;
}



void RunCA::ReadInfiltrationParameters()
{
    cout<<"\nStart reading infiltration parameters..."<<endl;
    string file_name1(model_name_);
    file_name1 += "-input-infiltration.txt";
    ifstream ifs1(file_name1.c_str());
    if( !ifs1.is_open() ) {
        cerr<<"input infiltration file ( " << file_name1 << " ) does not exist, please check"<<endl;
        terminate();
    } else{
        cout<<"  infiltration file opened successfully"<<endl;
    }

    //read infiltration file
    string text_line, token;
    getline(ifs1,text_line); //1st title line
    if(debug_) {
        cout<<"  reading following parameters:"<<endl; 
        cout<<"  "<<text_line<<endl;
    }
    
    //start reading parameters line by line
    int row(0);
    while(getline(ifs1,text_line) && text_line.length()!=0){
        row++;
        vector<double> row_para;
        istringstream iss1 (text_line);
        while(iss1 >> token)  row_para.push_back(stod(token));

        if(row_para.size() != 9) {
            cerr<<"  ERROR: wrong number of infiltration parameters: "<<row_para.size()<<" , should be 9"<<endl;
            terminate();
        }
        if((int)row_para[0] != row) {
            cerr<<"  ERROR: soil type ("<<(int)row_para[0]<<" } does not match row number ("<<row<<")"<<endl;
            cerr<<"  soil types need to be ordered row by row, please check input infiltration file"<<endl;
            terminate();            
        }
        iparams_.push_back(row_para);
    }
    ifs1.close();
    if(debug_)  cout<<"  read in infiltration parameters for "<<iparams_.size()<<" soil types"<<endl;
    cout<<"Finish reading infiltration parameters."<<endl;
}


//**********************
//void RunCA::Readerosionpara(){
void RunCA::ReadErosionParameters(){ 
	cout << "\nStart reading erosion parameters..." << endl;
	string file_name1(model_name_);
	file_name1 += "-input-erosion.txt"; 
	ifstream ifs4(file_name1.c_str());
	if (!ifs4.is_open()) {
		cerr << "input erosion file ( " << file_name1 << " ) does not exist, please check" << endl;
		terminate();
	}
	else {
		cout << "  input erosion file opened successfully" << endl;
	}
	//read sizefraction file
	string text_line, token;
	getline(ifs4, text_line); //1st title line
	if (debug_) {
		//cout << "  reading following size:" << endl;
        cout << "  reading following size fractions:" << endl; 
		cout << "  " << text_line << endl;
	}
	//start reading size line by line
	int row(0);
	while (getline(ifs4, text_line) && text_line.length() != 0) {
		row++;
		vector<double> row_i;
		istringstream iss1(text_line);
		while (iss1 >> token)  row_i.push_back(stod(token));
		if (row_i.size() != 11.) {
			cerr << "  ERROR: wrong number of erosion parameters : " << row_i.size() << " , should be 11." << endl;
			terminate();
		}
		if ((int)row_i[0] != row) {
			cerr << "  ERROR: soil type (" << (int)row_i[0] << " } does not match row number (" << row << ")" << endl;
			
            cerr << "  soil types need to be ordered row by row, please check input erosion file" << endl; 
			terminate();
		}
		erosion_para_.push_back(row_i);
	}

    for(int i=0;i<erosion_para_.size();i++){
        for(int j=0;j<erosion_para_[i].size();j++)
            cout<<erosion_para_[i][j]<<" ";
        cout<<endl;
    }    
    cout<<endl;   
      	
	ifs4.close();
	if (debug_)  cout << "  read in erosion parameters for " << erosion_para_.size() << " soil types" << endl;
	cout << "Finish reading erosion parameters." << endl;


}

//*************************
void RunCA::ReadRainfallRates()
{
    cout<<"\nStart reading input rainfall rates..."<<endl;
    string file_name(model_name_);
    file_name += "-input-rainfall.txt";
    ifstream ifs(file_name.c_str());
    if( !ifs.is_open() ) {
        cerr<<"input rainfall file ( " << file_name << " ) does not exist, please check"<<endl;
        terminate();
    } else{
        cout<<"  rainfall file opened successfully"<<endl;
    }   
    //read rainfall file
    string text_line, token;
    getline(ifs,text_line); //1st title line
    if(debug_) {
        cout<<"    rainfall file header:"<<endl; 
        cout<<"    "<<text_line<<endl;
    }
    //start reading rainfall rate line by line
    int minute(0);
    double sum(0.);
    rainfall_.push_back(0.); //Qi added - 6/11: 0 at 0 min
    while(getline(ifs,text_line) && text_line.length()!=0){
        minute++;
        istringstream iss (text_line);
        iss >> token;
        if(stoi(token) != minute){
            cerr<<"    ERROR: there seems to be missing data in rainfall input file, rainfall rates must be specified in consecutive minutes"<<endl;
            terminate();
        }
        iss >>token;
        rainfall_.push_back(stod(token));
        sum += stod(token);
    }
    ifs.close();
    if(rainfall_.size() != minute+1){
        cerr<<"    ERROR: rainfall vector size = "<<rainfall_.size()<<" does not entries in input rainfall file = "<<minute<<endl;
        cerr<<"    Please check your input rainfall file"<<endl;
        terminate();
    }
    duration_ = minute*1.0;
    if(debug_)  cout<<"    read in rainfall rates for "<<minute<<" minutes, mean rainfall rate = "<<sum/minute<<" mm/h"<<endl;
    cout<<"Finish reading input rainfall rates."<<endl;
}


void RunCA::TimeIntegration()
{
    cout<<"\n------Start time integration------"<<endl;

    
    int choice(0);
    cout<<"\nDo you want to include bendload transport? (1 = yes, 0 = no):"<<endl;
    cin>>choice;
    if(choice==1) with_bedload_transport_ = true;
    cout<<"with_bedload_transport = "<<with_bedload_transport_<<endl;

    //number of threads for computations
    threads_ = omp_get_max_threads();
    cout<<"\nPlease enter number of threads used for simulation:"<<endl;
    cout<<"(maximum available threads = "<<omp_get_max_threads()<<")"<<endl;
    cin>>threads_;
    threads_=max(1, threads_);
    threads_=min(omp_get_max_threads(), threads_);
    cout<<"Using "<<threads_<<" out of "<<omp_get_max_threads()<<" available threads"<<endl;
    
    //number of savings for results visualisation
    int num_savings(15);
    cout<<"\nPlease enter number of savings for results visualisation (duration="<<duration_<<" min): "<<endl;
    cout<<"(suggested number: 10~30, default = 15)"<<endl;
    cin>>num_savings;
    double save_interval = duration_*60/num_savings;
    cout<<"resulted saving interval = "<<save_interval<<" sec = "<<save_interval/60.<<" min"<<endl;

    //number of savings for results plotting
    int num_outputs(100);
    cout<<"\nPlease enter number of savings for results plotting (duration="<<duration_<<" min): "<<endl;
    cout<<"(suggested number: 50~500, default = 200)"<<endl;
    cin>>num_outputs;
    double output_interval = duration_*60/num_outputs;
    cout<<"resulted outputting interval = "<<output_interval<<" sec = "<<output_interval/60.<<" min"<<endl;    

    //estimatation of time increment based on 1 mm water depth
    double dt_5_percent = ComputeGlobalTimeStepSize(true);
    /*
    vector<double> min_time_steps(threads_, duration_*60.);  
    #pragma omp parallel num_threads(threads_)
    {
        int id = omp_get_thread_num();
        double min_step(duration_*60.);
        #pragma omp for collapse(2)        
        for(int y=1;y<ny_+1;y++){
            for(int x=1;x<nx_+1;x++){
                Cell* cell = cells_[y][x];               
                if(cell->GetValid()) {
                    double ori_depth = cell->GetWaterDepth();
                    cell->SetWaterDepth(1.); //1mm
                    double time_step = ComputeTimeStepSize(cell);
                    min_step = min(min_step, time_step);
                    cell->SetWaterDepth(ori_depth); //restore original depth
                }
            }
        }
        min_time_steps[id] = min_step;
    }
    sort(min_time_steps.begin(), min_time_steps.end());
    */
    
    
    //ajusting coefficient for time step size
    adjust_coefficient_ = 1.;
    cout<<"\nPlease enter an ajusting coefficient for time step size:"<<endl//
        //<<"(estimated time step size based on 1mm water depth is: "<<min_time_steps[0]<<" sec)"<<endl;
        <<"(estimated time step size based on 1mm water depth is: "<<dt_5_percent<<" sec)"<<endl;
    cin>>adjust_coefficient_;

    //time_increment_=min_time_steps[0] * adjust_coefficient_; //added
    time_increment_= dt_5_percent * adjust_coefficient_; //added

    //time integration
    cout<<"starting time integration..."<<endl;    
    double T_begin = omp_get_wtime();

    iteration_ = 0;
    model_time_ = 0.;
    double save_time(0.), output_time(0.);
    double previous_time(model_time_); //0 sec
    //double max_time_increment(60.); //60 sec
    //double max_time_increment(min_time_steps[0]);
    //double min_time_increment(1.0e-6); //1e-6 sec //Qi added: 17/12/2020
    //double max_time_increment = min_time_steps[0]*10.;
    //double min_time_increment = min_time_steps[0] * 0.001;
    double max_time_increment = dt_5_percent * 10.;
    double min_time_increment = dt_5_percent * 0.001;    
    int percentage(0);
    bool finished(false);
    WriteOutputFileHeaders();
    
    while(!finished) 
    { 
       // cout<<"passed 0"<<endl;
        //if(abs(model_time_-duration_*60.)<1.0e-4) finished = true;
        
        //0. check if there is ponded cell
        bool no_ponded_cell(true);
        for(int y=1;y<ny_+1 && no_ponded_cell;y++){
            for(int x=1;x<nx_+1 && no_ponded_cell;x++){
                Cell* cell = cells_[y][x];               
                if(cell->GetValid() && cell->GetWaterDepth()>0.) 
                    no_ponded_cell = false;
            }
        }

       // cout<<"passed 0.1"<<endl;
        
        //if no ponded cell, then compute time to ponding
        if(no_ponded_cell){ 
            cout<<"  no water on land surface, skipping dry period..."<<endl;
            cout<<"    computing time to ponding..."<<endl;
            vector<double> min_times(threads_, duration_*60.); 
            #pragma omp parallel num_threads(threads_)
            {
                double min_time(duration_*60.);
                int id = omp_get_thread_num();
                #pragma omp for collapse(2)
                for(int y=1;y<ny_+1;y++){
                    for(int x=1;x<nx_+1;x++){
                        Cell* cell = cells_[y][x];
                        if(cell->GetValid())
                            min_time = min(min_time, ComputeTimeToPonding(cell));
                    }
                }
                min_times[id]=min_time;
            }
            
            sort(min_times.begin(), min_times.end());
            double min_ponding_time = min_times[0]; 
         
            if(min_ponding_time == duration_*60.){
                cout<<"    no ponding will occur before the end of simulation, no need to continue"<<endl;
                //added
                #pragma omp parallel num_threads(threads_)
                {
                    #pragma omp for collapse(2)
                    for(int y=1;y<ny_+1;y++){
                        for(int x=1;x<nx_+1;x++){
                            Cell* cell = cells_[y][x];
                            //set flows to zeros
                            vector<double> flows(9,0.);
                            cell->SetFlows(flows);
                            cell->SettotalSedimentTransport(flows);
                            if(with_bedload_transport_) cell->SettotalSedimentTransportb(flows);
                        }
                    }
                } 
                          

                SaveWaterDepthToText();
                SaveResultsToVti();
                cout<<"    bringing all cell states to simulation end time..."<<endl;
                time_increment_ = output_interval; //added
                while(model_time_<duration_*60.) { //added
                //time_increment_ = duration_*60.-model_time_;
                vector<double> sum_rain(threads_, 0.), sum_infiltration(threads_,0.), sum_drainage(threads_,0.);
                #pragma omp parallel num_threads(threads_)
                {
                    int id = omp_get_thread_num();
                    double sum_r(0.), sum_i(0.), sum_d(0.);
                    #pragma omp for collapse(2)
                    for(int y=1;y<ny_+1;y++){
                        for(int x=1;x<nx_+1;x++){
                            Cell* cell = cells_[y][x];
                            if(cell->GetValid()) UpdateStates(cell,sum_r,sum_i,sum_d);  
                        }     
                    }
                    sum_rain[id]=sum_r; sum_infiltration[id]=sum_i; sum_drainage[id]=sum_d;
                }
                

                for(int n=0;n<threads_;n++){
                    tot_rain_+=sum_rain[n]; tot_infiltration_+=sum_infiltration[n]; tot_drainage_+=sum_drainage[n];
                }
                
                //model_time_ = duration_*60.;   
                WriteResults(); //Qi added - 6/11     
                model_time_ += output_interval; //added         
                }
                break;
            }
            

            cout<<"    computed time to ponding = "<<min_ponding_time<<" sec, bringing all cell states to this time level..."<<endl;
            time_increment_ = min_ponding_time;
            
            vector<int> updated_cells(threads_, 0), ponded_cells(threads_, 0); 
            vector<double> sum_rain(threads_, 0.), sum_infiltration(threads_,0.), sum_drainage(threads_,0.);
            #pragma omp parallel num_threads(threads_)
            {
                int id = omp_get_thread_num();
                int updated(0), ponded(0);
                double sum_r(0.), sum_i(0.), sum_d(0.);
                #pragma omp for collapse(2)
                for(int y=1;y<ny_+1;y++){
                    for(int x=1;x<nx_+1;x++){
                        Cell* cell = cells_[y][x];
                        if(cell->GetValid()) {
                            UpdateStates(cell,sum_r,sum_i,sum_d);
                            if(cell->GetWaterDepth()>0.) ponded++;
                            updated++;
                        }
                    }     
                }
                updated_cells[id]=updated; ponded_cells[id]=ponded;
                sum_rain[id]=sum_r; sum_infiltration[id]=sum_i; sum_drainage[id]=sum_d;
            }
            

            int cell_updated(0), cell_ponded(0);
            for(int n=0;n<threads_;n++){
                cell_updated += updated_cells[n]; cell_ponded += ponded_cells[n];
                tot_rain_+=sum_rain[n]; tot_infiltration_+=sum_infiltration[n]; tot_drainage_+=sum_drainage[n];
            }
           
            cout<<"    updated "<<cell_updated<<" cells at "<<min_ponding_time<<" sec, with "<<cell_ponded<<" ponded cells"<<endl; 
            //previous_time=model_time_;
            model_time_ = min_ponding_time; 
            cout<<"  skipped "<<model_time_-previous_time<<" sec dry period."<<endl;   
            previous_time=model_time_;  
        }  

       // cout<<"passed 0.7"<<endl;     


        //Qi: 11/12/2020
        //added - begin
        vector<double> sum_rain(threads_, 0.), sum_infiltration(threads_,0.), sum_drainage(threads_,0.); //added


        
        //1. compute a new time increment, dt
        /*
        vector<double> min_steps(threads_, duration_*60.-model_time_);  
        #pragma omp parallel num_threads(threads_)
        {
            int id = omp_get_thread_num();
            double min_step(duration_*60.-model_time_);
            #pragma omp for collapse(2)        
            for(int y=1;y<ny_+1;y++){
                for(int x=1;x<nx_+1;x++){
                    Cell* cell = cells_[y][x];               
                    if(cell->GetValid() && cell->GetWaterDepth()>0.) {
                        double time_step = ComputeTimeStepSize(cell);
                        min_step = min(min_step, time_step);
                    }
                }
            }
            min_steps[id] = min_step;
        }
        sort(min_steps.begin(), min_steps.end());
        double min_time_step = min_steps[0];
        */
        
        double min_time_step = ComputeGlobalTimeStepSize(false);
        
        //Qi added: 11/12/2020
        min_time_step *= adjust_coefficient_;  
        time_increment_ = min_time_step;
        time_increment_ = min(time_increment_, max_time_increment);
        time_increment_ = max(time_increment_, min_time_increment); //Qi added: 17/12/2020 
        if((model_time_ + time_increment_)> duration_*60.) 
                time_increment_=duration_*60.-model_time_;     
        
        //cout<<"passed 1"<<endl;
        

        //2. update cell states based on dt
        #pragma omp parallel num_threads(threads_)
        {
            int id = omp_get_thread_num();
            double sum_r(0.), sum_i(0.), sum_d(0.);
            #pragma omp for collapse(2) 
            for(int y=1;y<ny_+1;y++){
                for(int x=1;x<nx_+1;x++){
                    Cell* cell = cells_[y][x];
                    if(cell->GetValid()) UpdateStates(cell,sum_r,sum_i,sum_d);
                }
            }
            sum_rain[id]=sum_r; sum_infiltration[id]=sum_i; sum_drainage[id]=sum_d;
        }    

        for(int n=0;n<threads_;n++){
            tot_rain_+=sum_rain[n]; tot_infiltration_+=sum_infiltration[n]; tot_drainage_+=sum_drainage[n];
        }  
        
        model_time_ += time_increment_;
        //previous_time=model_time_;
        if(abs(model_time_-duration_*60.)<1.0e-4) {cout<<"finished = true"<<endl; finished = true;}  

       // cout<<"passed 2"<<endl;
        //3. compute water flows during dt for each valid cell, based on updated water depth
        vector<double> sum_runoff(threads_, 0.);
        #pragma omp parallel num_threads(threads_)
        {
            int id = omp_get_thread_num();
            double runoff(0.);
            #pragma omp for collapse(2)
            for(int y=1;y<ny_+1;y++){
                for(int x=1;x<nx_+1;x++){
                    Cell* cell = cells_[y][x];
                    //reset flows to zeros
                    vector<double> flows(9,0.);
                    cell->SetFlows(flows);
                    if(cell->GetValid()) runoff+=ComputeFlows2(cell);
                    //if(cell->GetValid()) runoff+=ComputeFlows(cell);
                }
            }
            sum_runoff[id] = runoff;
        }
        for(int n=0;n<threads_;n++) tot_runoff_+=sum_runoff[n];
        //cout<<"passed 3"<<endl;
        
		
        //Qi: 11/12/2020
        //this needs to be computed before water depth is updated in AccumulateFlows() function
        //4. compute transported sediments during dt for each valid cell
		vector<double> sum_transportedsediment(threads_, 0.);
        #pragma omp parallel num_threads(threads_)
		{
			int id = omp_get_thread_num();
			double transportedsediment(0.);
            #pragma omp for collapse(2)
			for (int y = 1; y < ny_ + 1; y++) {
				for (int x = 1; x < nx_ + 1; x++) {
					Cell* cell = cells_[y][x];
					//reset Stot to zeros
					vector<double> Stot(9, 0.);
					cell->SettotalSedimentTransport(Stot);
					if (cell->GetValid()) {
                        ComputeFlowVelocities(cell); //i
                        ComputeSlopes(cell); //i
                        transportedsediment += ComputeSedimentTransport(cell);
                        if(with_bedload_transport_) ComputeBedloadTransport(cell);
                    }
				}
			}
					
			sum_transportedsediment[id] = transportedsediment;
		}
		for (int n = 0; n < threads_; n++) tot_transportedsediment_ += sum_transportedsediment[n];
        //cout<<"passed 4"<<endl;
        

        //5. accumulate water flows for each valid cell, water depth is updated after redistribution
        #pragma omp parallel num_threads(threads_)
        {
            #pragma omp for collapse(2)
            for(int y=1;y<ny_+1;y++){
                for(int x=1;x<nx_+1;x++){
                    Cell* cell = cells_[y][x];
                    if(cell->GetValid()) AccumulateFlows(cell);
                }
            }
        } 
       // cout<<"passed 5"<<endl;
        /*
        //4. compute transported sediments during dt for each valid cell
		vector<double> sum_transportedsediment(threads_, 0.);
        #pragma omp parallel num_threads(threads_)
		{
			int id = omp_get_thread_num();
			double transportedsediment(0.);
            #pragma omp for collapse(2)
			for (int y = 1; y < ny_ + 1; y++) {
				for (int x = 1; x < nx_ + 1; x++) {
					Cell* cell = cells_[y][x];
					//reset Stot to zeros
					vector<double> Stot(9, 0.);
					cell->SettotalSedimentTransport(Stot);
					if (cell->GetValid()) {
                        ComputeFlowVelocities(cell); //i
                        ComputeSlopes(cell); //i
                        transportedsediment += ComputeSedimentTransport(cell);
                    }
				}
			}
					
			sum_transportedsediment[id] = transportedsediment;
		}
		for (int n = 0; n < threads_; n++) tot_transportedsediment_ += sum_transportedsediment[n];        
        */

		//**************************************************
        /*
		//compute transported sediments for each valid cell****
		vector<double> sum_transportedsediment(threads_, 0.);
        #pragma omp parallel num_threads(threads_)
		{
			int id = omp_get_thread_num();
			double transportedsediment(0.);
            #pragma omp for collapse(2)
			for (int y = 1; y < ny_ + 1; y++) {
				for (int x = 1; x < nx_ + 1; x++) {
					Cell* cell = cells_[y][x];
					//reset Stot to zeros
					vector<double> Stot(9, 0.);
					cell->SettotalSedimentTransport(Stot);
					if (cell->GetValid()) {
                        ComputeFlowVelocities(cell); //i
                        ComputeSlopes(cell); //i
                        transportedsediment += ComputeSedimentTransport(cell);
                    }
				}
			}
					
			sum_transportedsediment[id] = transportedsediment;
		}
		for (int n = 0; n < threads_; n++) tot_transportedsediment_ += sum_transportedsediment[n];
        */
        
        
		//6. accumulate transported sediment for each valid cell
        #pragma omp parallel num_threads(threads_)
		{
         #pragma omp for collapse(2)
			for (int y = 1; y < ny_ + 1; y++) {
				for (int x = 1; x < nx_ + 1; x++) {
					Cell* cell = cells_[y][x];
					if(cell->GetValid()) {
                        AccumulateSedimentTransport(cell);
                        if(with_bedload_transport_) AccumulateBedloadTransport(cell);
                    }
                    
				}
			}
		}
       // cout<<"passed 6"<<endl;
        

        //Qi: 11/12/2020
        //7. update erosion states (settling, deposition, etc.) for each valid cell, after water depth is updated in step 5
        #pragma omp parallel num_threads(threads_)
        {
            #pragma omp for collapse(2) 
            for(int y=1;y<ny_+1;y++){
                for(int x=1;x<nx_+1;x++){
                    Cell* cell = cells_[y][x];
                    if(cell->GetValid()) {
                        UpdateErosionStates(cell);
                        if(with_bedload_transport_) UpdateErosionStatesByBedload(cell);
                    }
                    
                }
            }
        }
        
        // cout<<"passed 7"<<endl;
         

        
        //8. output results if required
        if(model_time_ >= save_time){
            SaveWaterDepthToText();
            SaveResultsToVti();
            save_time += save_interval;
        }
        
        if(model_time_ >= output_time){
            WriteResults();
            output_time += output_interval;
        }
        //cout<<"passed 8"<<endl;

        iteration_++;

        //number of ponded cells and mean water depth
        int wet_cells(0);
        double mean_depth(0.), depth_sum(0.);
        for(int y=1;y<ny_+1;y++){
            for(int x=1;x<nx_+1;x++){
                Cell* cell = cells_[y][x];
                if(cell->GetValid() && cell->GetWaterDepth()>0.){
                    wet_cells++;
                    depth_sum += cell->GetWaterDepth();
                }
            }
        }  
        if(wet_cells>0) mean_depth=depth_sum/wet_cells;  

        //output to screen if required
        if(model_time_/(duration_*60.) >= percentage/100.) {
            cout<<"  iteration = "<<iteration_<<", time = "<<model_time_/60.<<" min, time step = "<<model_time_-previous_time<<" sec,"
            <<" ponded cells = "<<wet_cells<<", mean depth = "<<mean_depth<<" mm, "
            <<model_time_/(duration_*60.)*100.<<"% completed"<<endl;
            percentage = ceil(model_time_/(duration_*60.)*100.);
        }     

        /*
        //compute new time increment
        vector<double> min_steps(threads_, duration_*60.-model_time_);  
        #pragma omp parallel num_threads(threads_)
        {
            int id = omp_get_thread_num();
            double min_step(duration_*60.-model_time_);
            #pragma omp for collapse(2)        
            for(int y=1;y<ny_+1;y++){
                for(int x=1;x<nx_+1;x++){
                    Cell* cell = cells_[y][x];               
                    if(cell->GetValid() && cell->GetWaterDepth()>0.) {
                        double time_step = ComputeTimeStepSize(cell);
                        min_step = min(min_step, time_step);
                    }
                }
            }
            min_steps[id] = min_step;
        }
        sort(min_steps.begin(), min_steps.end());
        double min_time_step = min_steps[0];

        time_increment_ = min_time_step * adjust_coefficient_;
        time_increment_ = min(time_increment_, max_time_increment);

        if((model_time_ + time_increment_)> duration_*60.) 
            time_increment_=duration_*60.-model_time_;
        */
        
        /*
        //Qi added: update erosion states first
        //update erosion states for each valid cell - i
        #pragma omp parallel num_threads(threads_)
        {
            #pragma omp for collapse(2) 
            for(int y=1;y<ny_+1;y++){
                for(int x=1;x<nx_+1;x++){
                    Cell* cell = cells_[y][x];
                    if(cell->GetValid()) UpdateErosionStates(cell);
                }
            }
        } 
        */
        
        /*
        //update cell states for each valid cell
        vector<double> sum_rain(threads_, 0.), sum_infiltration(threads_,0.), sum_drainage(threads_,0.);
        #pragma omp parallel num_threads(threads_)
        {
            int id = omp_get_thread_num();
            double sum_r(0.), sum_i(0.), sum_d(0.);
            #pragma omp for collapse(2) 
            for(int y=1;y<ny_+1;y++){
                for(int x=1;x<nx_+1;x++){
                    Cell* cell = cells_[y][x];
                    if(cell->GetValid()) UpdateStates(cell,sum_r,sum_i,sum_d);
                }
            }
            sum_rain[id]=sum_r; sum_infiltration[id]=sum_i; sum_drainage[id]=sum_d;
        }    
        for(int n=0;n<threads_;n++){
            tot_rain_+=sum_rain[n]; tot_infiltration_+=sum_infiltration[n]; tot_drainage_+=sum_drainage[n];
        } 
        */  

        /*
        //update erosion states for each valid cell - i
        #pragma omp parallel num_threads(threads_)
        {
            #pragma omp for collapse(2) 
            for(int y=1;y<ny_+1;y++){
                for(int x=1;x<nx_+1;x++){
                    Cell* cell = cells_[y][x];
                    if(cell->GetValid()) UpdateErosionStates(cell);
                }
            }
        } 
        */
		
        previous_time = model_time_;
        //model_time_ += time_increment_;


    } //end while loop

    //output final results
    SaveWaterDepthToText();
    SaveResultsToVti();
    //WriteResults();

    run_time_ = omp_get_wtime()-T_begin; //sec
    WriteSummary();

    cout<<"------Time integration finished in "<<run_time_/60.<<" min, with "<<iteration_<<" iteration steps------"<<endl;
}


double RunCA::ComputeTimeToPonding(Cell* cell)
{
    //retrieve infiltration parameters
    int soiltype = cell->GetSoilType();
    vector<double> params = iparams_[soiltype-1];
    double If = params[1]; //steady infiltration rate
    double I0 = params[2]; //initial infiltration rate
    double P  = params[3]; //infiltration decay coefficient
    double TP = params[4]; //total porosity
    double WC0 = params[5]; //initial soil water content
    double FC = params[6]; //filed capacity
    double D = params[7]; //control zone depth: mm
    double N = params[8]; //manning's coefficent

    //determine time to ponding
    double water_depth = cell->GetWaterDepth(); //orign depth (mm) = 0
    double WS = cell->GetStoragePotential(); //origin soil storage potential
    double time(model_time_);
    while(water_depth == 0. && time<=duration_*60.) {
    //while(water_depth < 1.0 && time<=duration_*60.) {
        time += 1.0;
        //time += 30.0; //added
        //rainfall rate
        //int minute = round(time/60.);
        int minute = ceil(time/60.); //added
        //int minute = round(time/60.)+1; //added
        //double rain = rainfall_[minute-1]; //mm/h
        double rain = rainfall_[minute]; //mm/h
        //add rain to depth
        water_depth += rain/3600.; //mm
        //potential infiltration rate (infiltration capacity)
        double Ip=I0;
        if (WS==0) {Ip=If;}
        else {Ip=If+(I0-If)*pow(WS/(TP*D),P);}; //mm/h
        //available infiltration rate
        double Ia=water_depth/1.*3600.0; //mm/h
        //effective infiltration rate
        double Ie = (Ia < Ip) ? Ia : Ip; //mm/h
        //substract infiltration from depth
        water_depth -= Ie/3600.; //mm
        if(water_depth < 0.) water_depth = 0.;
        //drainage rate
        double Dr(0.);
        if(WS<(TP-FC)*D) { //water storage capacity is exceeded
            Dr=If*pow((1-WS/((TP-FC)*D)),3.); //mm/h
        }    
        //update soil water storage potential
        WS += (Dr-Ie)/3600.;
        if (WS>TP*D) {WS=TP*D;};
        if (WS<0.) {WS=0;};  

    }

    return time;
}

//***

void RunCA::UpdateErosionStates(Cell* cell) {
	
    //Qi: rule of thumb: use underscore(_) for class member variables only
    

	//retrieve erosion parameters
	int soiltype = cell->GetSoilType();
	vector<double> ero_prams = erosion_para_[soiltype - 1];
	double clay = ero_prams[1];
	double silt = ero_prams[2];
	double sand = ero_prams[3];
	double gravel = ero_prams[4];
	double k = ero_prams[5];
	double c = ero_prams[6];
	double p = ero_prams[7];
	//double ws = ero_prams[8]; //

    vector<double> SusVol = cell->GetSuspendedSedimentVolumes(); //^
    vector<double> BmVol = cell->GetBedMaterialVolumes(); //^
    vector<double> EroVol = cell->GetErodedMaterialVolumes(); //^
    
    //Qi added
    vector<double> SusInflow = cell->GetSusSediInflow();
    vector<double> BedInflow = cell->GetBmSediInflow();
    vector<double> EroInflow = cell->GetEroSediInflow();

    vector<double> SusOutflow = cell->GetSusSediOutflow();
    vector<double> BedOutflow = cell->GetBmSediOutflow();
    vector<double> EroOutflow = cell->GetEroSediOutflow();   

    
    //Qi added
    double eroded_volume(0.);
    for (int i = 0; i < 4; i++){
        //incoming suspended, bed and eroded sediments (from other cells) need to be added/transferred to suspended sediment volume
        SusVol[i] += (SusInflow[i] + BedInflow[i] + EroInflow[i]); 
        //outgoing suspended sediment needs to be substracted from suspended sediment volume 
        SusVol[i] += SusOutflow[i]; //note: SusOutflow is already negative
        //outgoing bed material needs to be substracted from bed sediment volume
        BmVol[i] += BedOutflow[i]; //note: BedOutflow is already negative 
        //outgoing bed and eroded sediment need to be substracted from land elevation later
        eroded_volume +=  (BedOutflow[i] + EroOutflow[i]); //note: negative value  
    }
    


	//compute settling and adjust suspended sediment and deposited/bed material
	double water_depth = cell->GetWaterDepth();
    water_depth /= 1000.; //Qi added: convert to m
   
    vector<double> ws(4, 0.); //fall velocities 
    //Table 2-1, Page24
    //ωs=1/18 [(d^2*g)/µ *(z-1)] 
/*
    vector<double> d(4, 0.);
    d[0] = 0.001/10.; //clay mm->cm 0.0001
    d[1] = 0.016/10.; //silt mm->cm  0.0016
    d[2] = 0.35/10.; //sand  mm->cm  0.035
    d[3] = 3/10.; //graval mm->cm  0.3    
    double g=980 ;//cm/s2;
    double o=0.01003;// cm2/s;
    double z=2.65;//2.65 gr/cm3     
    
    for (int i = 0; i < 4; i++){
        ws[i]=((1./18.)*(((pow(d[i],2)*g)/o)*(z-1)))*0.01;  //cm/s to m/s
    }
*/


     
   // ws[0] = 1.5; //clay
   // ws[1] = 2; //silt
    ws[0] = 8.6e-7; //clay (M/S)
    ws[1] = 2.2e-4; //silt
    ws[2] = 0.036; //sand
   ws[3] = 0.36; //graval (?)
     
    
    double deposited_volume(0.); //Qi added
    for (int i = 0; i < 4; i++) {
     
        //// h = h * (2./3.); //added                     
        double settlperc_i(0.); 
        if (water_depth > (ws[i] * time_increment_)) {   
            settlperc_i = ws[i] * time_increment_ / water_depth;
            //settlperc_i = ws[i] * time_increment_ / h;//sh added to use modify water depth
		} else {
            settlperc_i= 1.; 
        }
        
        double settled_vol = SusVol[i] * settlperc_i;
        
        //settled sediment needs to be substracted from suspended sediment volume
        SusVol[i] -= settled_vol; 
        //settled sediment needs to be added to bed sediment volume
        //AccumSus[i] -= settled_vol; 
        BmVol[i] += settled_vol; 
        //if(isnan(AccumBm[i]) ) cout<<"before update, AccumBm[i] = "<<AccumBm[i]<<" i = "<<i<<endl;
        //AccumBm[i] += settled_vol; 
        //BmVol[i] += AccumBm[i]; //Qi added: update bed material volumes
        //settled sediment needs to be added to land elevation later
        deposited_volume += settled_vol; //note: positive value
        //if(isnan(AccumBm[i]) ) cout<<"after update, AccumBm[i] = "<<AccumBm[i]<<" i = "<<i<<" settled_vol = "<<settled_vol<<endl;

        //if(cell->GetY() == 5 && cell->GetX() == 5) //Qi added
          //cout<<"  after settlement, SusVol[i] = "<<SusVol[i]<<", AccumSus[i] = "<<AccumSus[i]<<" BmVol[i] = "<<BmVol[i]<<" AccumBm[i] = "<<AccumBm[i]<<", i = "<<i<<endl<<endl;        
       
    } 

   

    //Store updated suspended sediment, bed material and eroded material
    cell->SetSuspendedSedimentVolumes(SusVol); 
    cell->SetBedMaterialVolumes(BmVol); 
    cell->SetErodedMaterialVolumes(EroVol);    

	//erosion or deposition ->update elevation
	//the balance between insed and outsed
	//balance between inflows and outflows: + for deposition, - for erosion
    
	/*
    double deposition (0.), erosion(0.); 
    for (int i = 0; i < 4; i++) { 
        //if(isnan(deposition)) cout<<"before update, deposition is nan, i = "<<i<<endl;
        deposition += AccumBm[i];
        //if(isnan(deposition)) cout<<"after update, deposition is nan, AccumBm[i] = "<<AccumBm[i]<<" i = "<<i<<endl;
        erosion += AccumEro[i]; 
    } 
    */
    //where y is the rate of erosion if (y<0) or deposition if (y>0) per unit area->hight=volume/area of cell
    //QS: this is wrong, erosion positive means deposition, negative means erosion
    //double y = (deposition - erosion) / (cell_size_ *cell_size_); //i
    //double y = (deposition + erosion) / (cell_size_ *cell_size_); //QS: fixed

    double y = (deposited_volume + eroded_volume) / (cell_size_ *cell_size_); //Qi added
    //if(y>1) cout<<"y = "<<y<<" deposition = "<<deposition<<" erosion = "<<erosion<<" cell_size_ = "<<cell_size_<<endl;
	double elevation = cell->GetLandElevation();
	elevation += y;
    cell-> SetLandElavationChange(y);
    //if(isnan(y)) cout<<"y = "<<y<<" deposition = "<<deposition<<" erosion = "<<erosion<<endl;

    double cumulative_elevation_change = cell->GetCumulativLandElavationChange();
    cumulative_elevation_change += y;
    cell->SetCumulativeLandElavationChange(cumulative_elevation_change);

     double cumulative_elevation_changesus =cell->GetCumulativLandElavationChangesus();
    
     cumulative_elevation_changesus += y;
     cell->SetCumulativeLandElavationChangesus(cumulative_elevation_changesus);
				
	cell->SetLandElevation(elevation); //update land elevation						
	cell->SetWaterSurfaceElevation(elevation*1000. + water_depth); //updated water surface elevation (mm)); sh
    //cell->SetWaterSurfaceElevation(elevation*1000. +( water_depth* (2./3.))); 

}


//updateErosion state for bedload

void RunCA::UpdateErosionStatesByBedload(Cell* cell) {
	
    //Qi: rule of thumb: use underscore(_) for class member variables only
    

	//retrieve erosion parameters
	int soiltype = cell->GetSoilType();
	vector<double> ero_prams = erosion_para_[soiltype - 1];
	double clay = ero_prams[1];
	double silt = ero_prams[2];
	double sand = ero_prams[3];
	double gravel = ero_prams[4];
	double k = ero_prams[5];
	double c = ero_prams[6];
	double p = ero_prams[7];
	

    
    vector<double> BmVolb = cell->GetBedMaterialVolumesb(); //^
    vector<double> EroVolb = cell->GetErodedMaterialVolumesb(); //^
    
    //Qi added
   
    vector<double> BedInflowb = cell->GetBmSediInflowb();
    vector<double> EroInflowb = cell->GetEroSediInflowb();

  
    vector<double> BedOutflowb = cell->GetBmSediOutflowb();
    vector<double> EroOutflowb = cell->GetEroSediOutflowb();   

    
    //Qi added
    double eroded_volumeb(0.), deposited_volumeb(0.);
    for (int i = 2; i < 4; i++){
        //incoming suspended, bed and eroded sediments (from other cells) need to be added/transferred to suspended sediment volume
        BmVolb[i] += ( BedInflowb[i] + EroInflowb[i]); 
        //outgoing suspended sediment needs to be substracted from suspended sediment volume 
      
        //outgoing bed material needs to be substracted from bed sediment volume
        BmVolb[i] += BedOutflowb[i]; //note: BedOutflow is already negative 
        //outgoing bed and eroded sediment need to be substracted from land elevation later
        EroVolb[i] += EroOutflowb[i];
        deposited_volumeb += (BedInflowb[i] + EroInflowb[i]); 
        eroded_volumeb +=  (BedOutflowb[i] + EroOutflowb[i]); //note: negative value  
    }
    double water_depth = cell->GetWaterDepth();
     water_depth /= 1000.; //Qi added: convert to m
  
   
/*
	//compute settling and adjust suspended sediment and deposited/bed material
	double water_depth = cell->GetWaterDepth();
    water_depth /= 1000.; //Qi added: convert to m
  
        
    double deposited_volumeb(0.); //Qi added
    for (int i = 2; i < 4; i++) {
        
                               
        double settlperc_i(0.); 
        settlperc_i= BmVolb[i];
        
        double settled_vol =  settlperc_i;
        
        //settled sediment needs to be substracted from suspended sediment volume
        //SusVol[i] -= settled_vol;

        //settled sediment needs to be added to bed sediment volume
        BmVolb[i] += settled_vol; 
        //if(isnan(AccumBm[i]) ) cout<<"before update, AccumBm[i] = "<<AccumBm[i]<<" i = "<<i<<endl;
        
        //settled sediment needs to be added to land elevation later
        deposited_volumeb += settled_vol; //note: positive value
        //if(isnan(AccumBm[i]) ) cout<<"after update, AccumBm[i] = "<<AccumBm[i]<<" i = "<<i<<" settled_vol = "<<settled_vol<<endl;
            
       
    } 
*/

    //Store updated suspended sediment, bed material and eroded material
    cell->SetBedMaterialVolumesb(BmVolb); 
    cell->SetErodedMaterialVolumesb(EroVolb);    

	//erosion or deposition ->update elevation
	//the balance between insed and outsed
	//balance between inflows and outflows: + for deposition, - for erosion
   
    double yb = (deposited_volumeb + eroded_volumeb) / (cell_size_ *cell_size_); //Qi added
    //if(y>1) cout<<"y = "<<y<<" deposition = "<<deposition<<" erosion = "<<erosion<<" cell_size_ = "<<cell_size_<<endl;
	double elevation = cell->GetLandElevation();
	elevation += yb;
    cell-> SetLandElavationChange(yb);
    //if(isnan(yb)) cout<<"yb = "<<yb<<" deposition = "<<deposition<<" erosion = "<<erosion<<endl;

    double cumulative_elevation_change = cell->GetCumulativLandElavationChange();
     double cumulative_elevation_changeb =cell->GetCumulativLandElavationChangeb();
    

    cumulative_elevation_change += yb;
    cumulative_elevation_changeb+=yb;
    cell->SetCumulativeLandElavationChange(cumulative_elevation_change);
    cell->SetCumulativeLandElavationChangeb(cumulative_elevation_changeb);
				
	cell->SetLandElevation(elevation); //update land elevation						
	cell->SetWaterSurfaceElevation(elevation*1000. + water_depth); //updated water surface elevation (mm)); sh

}

//***
void RunCA::UpdateStates(Cell* cell, double& sum_r, double& sum_i, double sum_d) {
    //retrieve infiltration parameters
    int soiltype = cell->GetSoilType();
    vector<double> params = iparams_[soiltype-1];
    double If = params[1]; //steady infiltration rate
    double I0 = params[2]; //initial infiltration rate
    double P  = params[3]; //infiltration decay coefficient
    double TP = params[4]; //total porosity
    double WC0 = params[5]; //initial soil water content
    double FC = params[6]; //filed capacity
    double D = params[7]; //control zone depth: mm
    double N = params[8]; //manning's coefficent
	

    double water_depth = cell->GetWaterDepth(); //orign dept (mm)
    double WS = cell->GetStoragePotential(); //origin soil storage potential
    double Ip = cell->GetInfiltrationCapacity();
    double Ie = cell->GetInfiltrationRate();
    double Dr = cell->GetDrainageRate();

    double time(0.);
    while(time<time_increment_) {
        double dt(1.); //1 sec
        if(time+dt > time_increment_) dt=time_increment_-time;
        time += dt;             
        //rainfall rate
        //int minute = round((model_time_+time)/60.);
        int minute = ceil((model_time_+time)/60.); //added
        //int minute = round((model_time_+time)/60.) + 1; //added
        //double rain = rainfall_[minute-1]; //mm/h
        double rain = rainfall_[minute]; //mm/h
        sum_r += rain/3600.*dt;
        //add rain to depth
        water_depth += rain/3600.*dt; //mm
        //potential infiltration rate (infiltration capacity)
        if (WS==0) {Ip=If;}
        else {Ip=If+(I0-If)*pow(WS/(TP*D),P);}; //mm/h
        //available infiltration rate
        double Ia=water_depth/dt*3600.0; //mm/h
        //effective infiltration rate
        Ie = (Ia < Ip) ? Ia : Ip; //mm/h
        sum_i += Ie/3600.*dt;
        //substract infiltration from depth
        water_depth -= Ie/3600.*dt; //mm
        if(water_depth < 0.) water_depth = 0.;
        //drainage rate
        if(WS<(TP-FC)*D) { //water storage capacity is exceeded
            Dr=If*pow((1-WS/((TP-FC)*D)),3.); //mm/h
        } else{
            Dr=0.;
        }   
        sum_d += Dr/3600.*dt;  
        //update soil water storage potential
        WS += (Dr-Ie)/3600.*dt;
        if (WS>TP*D) {WS=TP*D;};
        if (WS<0.) {WS=0;};  
    }

    cell->SetWaterDepth(water_depth); //updated water depth
    //if(cell->GetY() == ny_) cell->SetWaterDepth(0.); //added
    cell->SetWaterSurfaceElevation(cell->GetLandElevation()*1000. + water_depth); //updated water surface elevation  
    cell->SetStoragePotential(WS); //updated soil water storage potential
    cell->SetInfiltrationCapacity(Ip);
    cell->SetInfiltrationRate(Ie);
    cell->SetDrainageRate(Dr);
}





double RunCA::ComputeFlows(Cell* cell) 
{        
    if(cell->GetWaterDepth() == 0.) { //there is no water to flow
        return 0.;
    }
    
    vector<double> flows(9,0.);

    //find neighbour cells
    vector<Cell*> neighbours;
    int y_c = cell->GetY();
    int x_c = cell->GetX();
    for(int y=y_c-1;y<=y_c+1;y++){
        for(int x=x_c-1;x<=x_c+1;x++){
            Cell* neighbour = cells_[y][x];
            neighbours.push_back(neighbour);
        }
    }
	
    //calculate cardinal average height:
    int count1(0); 
    double sum(0.), ave(0.);
    int soiltype = cell->GetSoilType();
    double manning_n = iparams_[soiltype-1][8];
    double water_depth = cell->GetWaterDepth(); //mm    
    double wh_c = cell->GetWaterSurfaceElevation(); //surface water elevation of central cell 
    for (int n=1;n<9;n=n+2) {
        if(neighbours[n]->GetValid()) {
            double wh_n = neighbours[n]->GetWaterSurfaceElevation(); //surface water elevation of neighbour cell
            if(wh_n < wh_c) {sum+=wh_n; count1++;}
        }
    }
    ave=(sum+wh_c)/(count1+1.);

    //if there is no cardinal neighbour to flow, then check diagonal neighbours
    if(count1==0){
        count1=0; sum=0;
        for (int n=0;n<9;n=n+2) {
            if(n!=4 && neighbours[n]->GetValid()) {
                double wh_n = neighbours[n]->GetWaterSurfaceElevation(); //surface water elevation of neighbour cell
                if(wh_n < wh_c) {sum+=wh_n; count1++;}
            }
        }
        ave=(sum+wh_c)/(count1+1.);
        //if there is no diagonal neighbour to flow either
        if(count1==0){
            return 0.;           
        } else { //there is only diagonal neighbour to flow
            //identify diagonal flow neighbours:
            int count(0);
            sum=0;
            for (int n=0;n<9;n=n+2) {
                if(n!=4 && neighbours[n]->GetValid()) {
                    double wh_n = neighbours[n]->GetWaterSurfaceElevation();
                    if(wh_n<ave) {sum+=wh_n; count++;}
                }
            }
            ave=(sum+wh_c)/(count+1.);
            while (count < count1) {
                count1=count;
                count=0; sum=0;
                for (int n=0;n<9;n=n+2) {
                    if(n!=4 && neighbours[n]->GetValid()) {
                        double wh_n = neighbours[n]->GetWaterSurfaceElevation();
                        if(wh_n<ave) {sum+=wh_n; count++;}
                    }
                }                
                if (count < count1 ) {ave=(sum+wh_c)/(count+1);};  
            }
            //calculate diagonal water flows
            for (int n=0;n<9;n=n+2) {
                if(n!=4 && neighbours[n]->GetValid()){
                    double wh_n = neighbours[n]->GetWaterSurfaceElevation();
                    if(wh_n<ave && wh_n<wh_c){
                        double flow = ave-wh_n; //mm
                        double time_step=cell_size_*1.414*manning_n/pow(water_depth/1000.,0.667)/pow((wh_c-wh_n)/1000./(cell_size_*1.414),0.5);
                        if(time_step<time_increment_) {flow=flow*time_step/time_increment_;};
                        flows[n] = flow;
                        flows[4] -= flow;
                    }
                }
            }
            // if there is not enough water to flow
            if ((-flows[4]) > cell->GetWaterDepth()) {
                for (int n=0;n<9;n=n+2) {
                    if (n!=4) {
                        double flow = flows[n];
                        flow *= (cell->GetWaterDepth() / (-flows[4]));
                        flows[n] = flow; 
                    }
                }
                flows[4] = -cell->GetWaterDepth();
            }                     
        }
    } 
    //there is cardinal neighbour(s) to flow, diagonal neighbours are excluded from computations 
    else { 
        //identify cardinal flow neighbours:
        int count(0);
        sum=0;
        for (int n=1;n<9;n=n+2) {
            if(neighbours[n]->GetValid()) {
                double wh_n = neighbours[n]->GetWaterSurfaceElevation();
                if(wh_n<ave) {sum+=wh_n; count++;}
            }
        }
        ave=(sum+wh_c)/(count+1.);
        while (count < count1) {
            count1=count;
            count=0; sum=0;
            for (int n=1;n<9;n=n+2) {
                if(neighbours[n]->GetValid()) {
                    double wh_n = neighbours[n]->GetWaterSurfaceElevation();
                    if(wh_n<ave) {sum+=wh_n; count++;}
                }
            }                
            if (count < count1 ) {ave=(sum+wh_c)/(count+1);};  
        }
        //calculate cardinal water flows
        for (int n=1;n<9;n=n+2) {
            if(neighbours[n]->GetValid()){
                double wh_n = neighbours[n]->GetWaterSurfaceElevation();
                if(wh_n<ave && wh_n<wh_c){
                    double flow = ave-wh_n; //mm
                    double time_step=cell_size_*manning_n/pow(water_depth/1000.,0.667)/pow((wh_c-wh_n)/1000./cell_size_,0.5);
                    if(time_step<time_increment_) {flow=flow*time_step/time_increment_;};
                    flows[n] = flow;
                    flows[4] -= flow;
                }
                    //double flow = ave-wh_n; //mm
            }
        }
        // if there is not enough water to flow
        if ((-flows[4]) > cell->GetWaterDepth()) {
            for (int n=1;n<9;n=n+2) {
                double flow = flows[n];
                flow *= (cell->GetWaterDepth() / (-flows[4]));
                flows[n] = flow;     
            }
            flows[4] = -cell->GetWaterDepth();
        }      
    }
    /*
    for (int y=ny_,x=1;x<nx_+1;x++) { //bottom outlet
        if (!cells_[y][x]->GetValid()) { //invalid cell
            continue;
        } else { //valid cell
            double h_above = cells_[y-1][x]->GetWaterSurfaceElevation();
            double h = cells_[y][x]->GetWaterSurfaceElevation();
            double h_below = h - (h_above - h);
            cells_[y+1][x]->SetWaterSurfaceElevation(h_below);
        };
    }; 
    */
    cell->SetFlows(flows);
    
    double runoff(0.);
    if(y_c==ny_) runoff += (flows[6]+flows[7]+flows[8]); //runoff amount (mm)

    return runoff;
}


//added
double RunCA::ComputeFlows2(Cell* cell) 
{        
    if(cell->GetWaterDepth() == 0.) { //there is no water to flow
        return 0.;
    }
    
    vector<double> flows(9,0.);

    //find neighbour cells
    vector<Cell*> neighbours;
    int y_c = cell->GetY();
    int x_c = cell->GetX();
    for(int y=y_c-1;y<=y_c+1;y++){
        for(int x=x_c-1;x<=x_c+1;x++){
            Cell* neighbour = cells_[y][x];
            neighbours.push_back(neighbour);
        }
    }
	
    //check lower cardinal neighbours:
    int count1(0); 
    int soiltype = cell->GetSoilType();
    double manning_n = iparams_[soiltype-1][8];
    double water_depth = cell->GetWaterDepth(); //mm  
    double wh_c = cell->GetWaterSurfaceElevation(); //surface water elevation of central cell 
    double h_c = cell->GetLandElevation(); //m
    for (int n=1;n<9;n=n+2) {
        if(neighbours[n]->GetValid()) {
            double wh_n = neighbours[n]->GetWaterSurfaceElevation(); //surface water elevation of neighbour cell
            if(wh_n < wh_c) {count1++;}
        }
    }

    //if there is no cardinal neighbour to flow, then check diagonal neighbours
    if(count1==0){
        count1=0; 
        for (int n=0;n<9;n=n+2) {
            if(n!=4 && neighbours[n]->GetValid()) {
                double wh_n = neighbours[n]->GetWaterSurfaceElevation(); //surface water elevation of neighbour cell
                if(wh_n < wh_c) {count1++;}
            }
        }
        //if there is no diagonal neighbour to flow either
        if(count1==0){
            return 0.;           
        } else { //there is only diagonal neighbour to flow
            //calculate diagonal water flows
            for (int n=0;n<9;n=n+2) {
                if(n!=4 && neighbours[n]->GetValid()){
                    double wh_n = neighbours[n]->GetWaterSurfaceElevation();
                    double h_n = neighbours[n]->GetLandElevation(); //m
                    //if(h_n<h_c){
                       if(wh_n<wh_c){ 
                        
                        double slope = (wh_c-wh_n)/1000./(cell_size_*1.414);
                       // double slope = (h_c-h_n)/(cell_size_*1.414);
                        double velocity = pow(water_depth/1000.,0.667) * pow(slope,0.5) / manning_n; //m/s
                        //double adjusted_n = manning_n * pow(water_depth/10., -1./3.); 
                        //double velocity = pow(water_depth/1000.,0.667) * pow(slope,0.5) / adjusted_n; //m/s
                        double q = velocity * water_depth/1000.; //m2/s
                        double flow = q * time_increment_ * cell_size_*1.414;; //m3
                        flows[n] = flow;
                        flows[4] -= flow;
                    }                    
                }
            }
            
            // if there is not enough water to flow
            double water_available = cell->GetWaterDepth()/1000.*cell_size_*cell_size_; //m3
            if ((-flows[4]) > water_available) { //m3
                for (int n=0;n<9;n=n+2) {
                    if (n!=4) {
                        double flow = flows[n]; //m3
                        flow *= (water_available / (-flows[4]));
                        flows[n] = flow; 
                    }
                }
                flows[4] = -water_available;
            }  
                               
        }
    } 
    //there is cardinal neighbour(s) to flow, diagonal neighbours are excluded from computations 
    else { 
        //calculate cardinal water flows
        for (int n=1;n<9;n=n+2) {
            if(neighbours[n]->GetValid()){
                double wh_n = neighbours[n]->GetWaterSurfaceElevation();
                double h_n = neighbours[n]->GetLandElevation(); //m
                //if(h_n<h_c){
                   if(wh_n<wh_c){ 
                    double slope = (wh_c-wh_n)/1000./cell_size_;
                    //double slope = (h_c-h_n)/cell_size_;
                    double velocity = pow(water_depth/1000.,0.667) * pow(slope,0.5) / manning_n; //m/s
                    //double adjusted_n = manning_n * pow(water_depth/10., -1./3.); 
                    //double velocity = pow(water_depth/1000.,0.667) * pow(slope,0.5) / adjusted_n; //m/s
                    double q = velocity * water_depth/1000.; //m2/s
                    double flow = q * time_increment_ * cell_size_; //m3
                    flows[n] = flow;
                    flows[4] -= flow;
                }
            }
        }
           
        // if there is not enough water to flow
        double water_available = cell->GetWaterDepth()/1000.*cell_size_*cell_size_; //m3
        if ((-flows[4]) > water_available) { //m3
            for (int n=1;n<9;n=n+2) {
                double flow = flows[n]; //m3
                flow *= (water_available / (-flows[4]));
                flows[n] = flow; 
            }
            flows[4] = -water_available;
        }          
    }

    cell->SetFlows(flows);
    
    double runoff(0.);
    //if(y_c==ny_) runoff += (flows[6]+flows[7]+flows[8]) / cell_size_ / cell_size_*effective_cells_* 1000; //runoff amount (mm)
     if(y_c==ny_) runoff += (flows[6]+flows[7]+flows[8]) / cell_size_ / cell_size_/effective_cells_* 1000; //runoff amount (mm)
    //if(y_c==ny_) runoff += (flows[6]+flows[7]+flows[8]) * 1000.; //runoff amount (L)

    return runoff;
}


void RunCA::AccumulateFlows(Cell* cell)
{
    int y = cell->GetY();
    int x = cell->GetX();

    //neighbourhood
    //(y-1, x-1)  (y-1, x)  (y-1, x+1)
    //( y , x-1)  ( y , x)  ( y , x+1)
    //(y+1, x-1)  (y+1, x)  (y+1, x+1)
 
    //flows
    //0 1 2
    //3 4 5
    //6 7 8

    // calculate accumulated flow (m3)
    double accum_flow = cells_[y+1][x+1]->GetFlows()[0] //down right neighbour, up left flow 
                      + cells_[y+1][x]->GetFlows()[1] //down neighbour, up flow
                      + cells_[y+1][x-1]->GetFlows()[2] //down left neighbour, up right flow
                      + cells_[y][x+1]->GetFlows()[3] //right neighbour, left flow
                      + cells_[y][x]->GetFlows()[4] //central cell, outflow //added - commented out
                      + cells_[y][x-1]->GetFlows()[5] //left neighbour, right flow
                      + cells_[y-1][x+1]->GetFlows()[6] //up right neighbour, down left flow
                      + cells_[y-1][x]->GetFlows()[7] //up neighbour, down flow
                      + cells_[y-1][x-1]->GetFlows()[8]; //up left neighbour, down right flow
    
    //if(y!=1) accum_flow += cells_[y][x]->GetFlows()[4]; //central cell, outflow //added
    
    double water_depth = cell->GetWaterDepth(); //mm
    //water_depth += accum_flow;
    water_depth += accum_flow / cell_size_ / cell_size_ * 1000.; //mm
    //water_depth += accum_flow / cell_size_ /100. / cell_size_ /100. * 10.; //mm
    //water_depth += accum_flow / cell_size_ /100.* 10.; //mm
    if(water_depth < 0.) water_depth = 0.;
    cell->SetWaterDepth(water_depth); //updated water depth (mm)
    cell->SetWaterSurfaceElevation(cell->GetLandElevation()*1000. + water_depth); //updated water surface elevation (mm)
}


//*********************************************************************************************************************************************************

//functions like flow velocity and slope ( digonal and cardinall)


void RunCA::ComputeFlowVelocities(Cell* cell) { 
    
    vector<double> fv(9, 0.); 
	//find neighbour cells
	vector<Cell*> neighbours;
	int y_c = cell->GetY();
	int x_c = cell->GetX();
	for (int y = y_c - 1; y <= y_c + 1; y++) {
		for (int x = x_c - 1; x <= x_c + 1; x++) {
			Cell* neighbour = cells_[y][x];
			neighbours.push_back(neighbour);
		}
	}
    
    double wh_c = cell->GetWaterSurfaceElevation(); 
	for (int n = 0; n < 9; n++) {
		double water_depth = cell->GetWaterDepth();
		int soiltype = cell->GetSoilType();
		double manning_n = iparams_[soiltype - 1][8];
	
		double wh_n = neighbours[n]->GetWaterSurfaceElevation();

		if (neighbours[n]->GetValid()) {
            double d_h = fabs((wh_c - wh_n)/1000.); //ensure positive 
            double distance = cell_size_; 
            if(n==0 || n==2 || n==6 || n==8) distance *= 1.414; ///digonal neighbours 
            fv[n] = (1./manning_n) * pow(water_depth/1000., 0.667) * pow(d_h/distance, 0.5); 
            
        }
        
    
	}
    cell->SetFlowVelocity(fv);  
}
	
//**********slope

void RunCA::ComputeSlopes(Cell*cell) { 
    
	vector<double> slope(9,0.);
	//find neighbour cells
	vector<Cell*> neighbours;
	int y_c = cell->GetY();
	int x_c = cell->GetX();
	for (int y = y_c - 1; y <= y_c + 1; y++) {
		for (int x = x_c - 1; x <= x_c + 1; x++) {
			Cell* neighbour = cells_[y][x];
			neighbours.push_back(neighbour);
		}
	}
    
    double z_c = cell->GetLandElevation(); 
	for(int n=0;n<9;n++){
	    if (neighbours[n]->GetValid()) {
		    
		    double z_n = neighbours[n]->GetLandElevation();
            double d_h = fabs(z_n - z_c); //ensure positive 
            double distance = cell_size_; 
            if(n==0 || n==2 || n==6 || n==8) distance *= 1.414; ///digonal neighbours 
            slope[n] = d_h / distance; 
            //if(slope[n]>1.) cout<<"slope[n] = "<<slope[n]<<" d_h = "<<d_h<<" z_n = "<<z_n<<" z_c = "<<z_c<<" distance = "<<distance<<" n = "<<n<<endl;
           
	    }
        
    }
    cell->SetSlope(slope);
  
}





//*********************************************
//the unit discharge rate of cell unit is equal to the product of groundwater flow velocity and the depth of water.
//q=V*h  //Where q is unit discharge rate (m2/min), V is the runoffvelocity(m / s), h is the water depth of each cell(mm)
//velocity=(1/n)*(h^2/3)*(s^1/2) //=qskr=58390*s^1.664*q(unit discharge[m2 s-1]^2.035*k*C*P*x*t

	
    //Bed load transport Yalin' equation
    //G=2.65, h=water depth, s=slope, y=water density, p=soil density   , g = 9.8 (m/s^2)

    //shear stress=y*h*s
    //shear velocity =(g*S * h)^0.5    unit-> meter/s
  //  Y: Dimensionless shear stress =  (sherstress/y)/((G-1)*g*particle size)
  //Ycr: dimensionless critical shear from Shield diagram
  //∆=(Y/Ycr)-1  (when Y<Ycr ,∆=0)
  //ბ = 2.45 * ∆ *(Ycr^0.5) * (G^-0.4)
    
    /*
    particlesize
    vector<double> d(4, 0.);
    d[0] = 0.001/10.; //clay mm->cm 0.0001        Ycr=0  
    d[1] = 0.016/10.; //silt mm->cm  0.0016       Ycr=10^-1 or 0.1
    d[2] = 0.35/10.; //sand  mm->cm  0.035        Ycr=0.15 (N/M^2)
    d[3] = 3/10.; //graval mm->cm  0.3            Ycr=2 (N/M^2)
    */ 
    //   Bedload transport capacity-> TC=0.635*10^6*G*p*y*particle diameter*shear velocity*∆*[1-ln((1+ბ)/ბ)]


    

//***********************************

//Function ComputeSedimentTransport(cell)
double RunCA::ComputeSedimentTransport(Cell* cell) {
	
    double Ssus=0,Sbm=0,Sero=0; //tota outflowing suspended sediment, bed material and eroded material

	//retrieve erosion parameters
	int soiltype = cell->GetSoilType();
	vector<double> ero_params = erosion_para_[soiltype - 1];
	double clay = ero_params[1];
	double silt = ero_params[2];
	double sand = ero_params[3];
	double gravel = ero_params[4];
	double cf = ero_params[5];
	double cr = ero_params[6];
	double p = ero_params[7];//solid density
    double y = ero_params[8];//i ??density of water constan
 
    double k= ero_params[9];  //rill erodibility(kg s/m-4)
    double c= ero_params[10];  //interrill erodibility (s m-1)
    double pr= ero_params[11];//pa
   
	//double ws = ero_params[8]; //i ??
  
	//find neighbour cells
	vector<Cell*> neighbours;
	int y_c = cell->GetY();
	int x_c = cell->GetX();
	for (int y = y_c - 1; y <= y_c + 1; y++) {
		for (int x = x_c - 1; x <= x_c + 1; x++) {
			Cell* neighbour = cells_[y][x];
			neighbours.push_back(neighbour);
		}
	}

    
    vector<double> Stot(9, 0.); //total sediment (from all sources and all fractions) transported to neighbours - i

	//1. compute total volume of outgoing sediment SKR
	//for each outgoing cell n :
	double SKR = 0; //total volume of outgoing sediment
	
    vector<double> skr_n(9, 0.); //volumes of outgoing sediment from different fractions - i ( different nighbour)
	double water_depth = cell->GetWaterDepth(); 

	for (int n = 0; n < 9; n++) {
		
		if (cell->GetValid()  && cell->GetFlows()[n] != 0.) {//there is flow to transport sediment
          if(n!=4){
            

			//compute total volume of outgoing sediment skr_n(eq 15)
            //i: make sure to use '.' for a double number, otherwise computed result may be an intergal
            //i: double check the equation
			//skr_n[n] = 58390 * pow(cell->GetSlope()[n], 1.664)*((water_depth / 1000) * cell->GetFlowVelocity()[n])*k*c*pr*cell_size_* time_increment_;
            //i: slopes need to be computed previously
            double S = cell->GetSlope()[n]; // m/m
            //QS: this is wrong. To convert flow from mm to m, it should be devided by 1000
            //double q = cell->GetFlows()[n] * 1000. * cell_size_ / time_increment_; // m2/s		
            //double h = cell->GetWaterDepth()*100;//shala addeded unit=cm
            double h = cell->GetWaterDepth()/10.;//Q: mm->cm
            //h = pow(h, 6./5.); //added
            //h = h * (2./3.); //added
            
            //int minute = round((model_time_)/60.); 
            int minute = ceil((model_time_)/60.); //added
            //int minute = round((model_time_)/60.) + 1; //added
            //double rain = rainfall_[minute-1]/10./60.; //  added mm/h to cm/min
            double rain = rainfall_[minute]/10./60.; //  added mm/h to cm/min
            double q = cell->GetFlows()[n] / 1000. * cell_size_ / time_increment_; // m2/s QS: fixed
            //skr_n[n] =((cf/p)*y* S * (h/10.)) + ((cr/p)*pow(rain,2)); //•	g/min/cm2 
            //analitical one
           
           



            //skr_n[n] = (cf / p) * y * S * h + (cr / p) * pow(rain,2.); //g/min/cm2
             skr_n[n] = (cf / p) * y * S * h + (cr / p) * pow(rain,2.); //g/min/cm2 

            //skr_n[n] *= ( (cell_size_*100.) * (time_increment_/60.)/p);//convert this g/min/cm2 to cm3	//(m3)try
            //double q = cell->GetFlows()[n] / 10. * (cell_size_*100.); // cm2 //added
            //skr_n[n] *= ( q * (time_increment_/60.) );//g/min/cm2->g //added
            //skr_n[n] /= p; //g->cm3 //added
           //original skr_n[n] *= ( (cell_size_*100.)* (cell_size_*100.) * (time_increment_/60.) / p );//Q: g/min/cm->cm3, a cell_size*100 was missing //commented out
           double distance = cell_size_; 
            if(n==0 || n==2 || n==6 || n==8) distance *= 1.414; ///digonal neighbours 
           
           //test  skr_n[n] *= ( (cell_size_*100.)* (time_increment_/60.) / p );//from gr/min/cm ->cm3
           skr_n[n] *= ( (distance*100.)* (time_increment_/60.) / p );//from gr/min/cm ->cm3
           
            skr_n[n] /= 1.0e6; //cm3 -> m3 //commented out	
          
       
          //
           // skr_n[n] = 58390. * pow(S, 1.664) * pow(q, 2.035) * k * c * pr * cell_size_ * time_increment_;	//(m3)try
            //if(skr_n[n]>1.) cout<<"skr_n[n] = "<<skr_n[n]<<" S = "<<S<<" q = "<<q<<" k = "<<k<<" c = "<<c<<" p = "<<p<<" cell_size_ = "<<cell_size_<<" time_increment_ = "<<time_increment_<<endl;
                  
			SKR += skr_n[n];//SKR = sum(skr_n) unit:(m3)
          }
		}
	}
       

	//2. compute transport of suspended sediment
  
    vector<double> Ssus_n(9,0.); //transported suspended sediment (in all fractions) to neighbour cells 
    
    //suspended sediment flow to neighbours from different fractions
    vector<vector<double>> Ssus_i_flow(4, vector<double>(9,0.)); 
     
	//compute transported suspended sediment Ssus_i_n(eq 18)
    //total suspended sediment volume 
    vector<double> SusVol = cell->GetSuspendedSedimentVolumes();
    vector<double>consen ;
         
    double tot_sus_vol = SusVol[0]+SusVol[1]+SusVol[2]+SusVol[3];
    
    if(tot_sus_vol > 0.) { 
	    
        for (int i = 0; i < 4; i++) { 
		    
            double Ssus_i(0.); //total outgoing suspended volume from the size fraction     
            for (int n = 0; n < 9; n++) { 

			  //if(n!=4){
              if(n!=4 && cell->GetFlows()[n]>0.){ //QS: I think there should be a check on whether there is a flow
               
               //QS: there is a missing condition to check whether skr_n[n] is less than tot_sus_vol, see eqn (18)
               if(skr_n[n] < tot_sus_vol) { //QS: added
                
                //Q: flow velocity needs to be computed previously
                double num_1 = SusVol[i] * cell->GetFlowVelocity()[n] * time_increment_ / cell_size_; 
               
                double num_2 = skr_n[n] * SusVol[i] / tot_sus_vol; 
                 
               
                Ssus_i_flow[i][n]= max(num_1, num_2);
               
               

               } else { //QS: added
                Ssus_i_flow[i][n] = SusVol[i]; //QS: added
               }
               
                //Ssus_i_flow[i][n] =  skr_n[n] * ero_params[i+1];
                Ssus_i += Ssus_i_flow[i][n];
                
              }
            }
           
            //if(Ssus_i > SusVol[i]) { //when outgoing volume is larger than available volume 
            //QS: performs only when Ssus_i>0
            if(Ssus_i > 0. && SusVol[i] > 0. && Ssus_i > SusVol[i]) { //when outgoing volume is larger than available volume
                for (int n = 0; n < 9 ; n++) //i
                if(n!=4){
                    //if(isnan(Ssus_i_flow[i][n])) cout<<"before update, Ssus_i_flow[i][n] = "<<Ssus_i_flow[i][n]<<" i = "<<i<<" n = "<<n<<endl;
                    Ssus_i_flow[i][n] *= (SusVol[i]/Ssus_i); //make according adjustment 
                    //if(isnan(Ssus_i_flow[i][n])) cout<<"after update, Ssus_i_flow[i][n] = "<<Ssus_i_flow[i][n]<<" SusVol[i] = "<<SusVol[i]<<" Ssus_i = "<<Ssus_i<<" i = "<<i<<" n = "<<n<<endl;
                }            
            } 
            

            for (int n = 0; n < 9; n++) {
              if(n!=4){
                Ssus_n[n] += Ssus_i_flow[i][n]; 
                
                Stot[n] += Ssus_i_flow[i][n]; 
           
                Ssus += Ssus_i_flow[i][n];
                
              }      
            } 
           
        } //end size fraction loop
    } //end if
   		   	
	//3. compute transport of bed material (deposited sediment)
	
    vector<double> Sbm_n(9, 0.); //transported bed material (in all fractions) to neighbour cells 
  //bedload sediment flow to each neighbour in each fraction
    vector<vector<double>> Sbm_i_flow(4, vector<double>(9,0.)); 
    
       //compute total excess capacity totXSScap(eq 19)
	double TotXssCap= max(0.0, (SKR- Ssus));
             
    vector<double> BmVol = cell->GetBedMaterialVolumes(); //bed material volumes in different fractions 
       
    double tot_bm_vol = BmVol[0]+BmVol[1]+BmVol[2]+BmVol[3]; //total material volume in all fractions 

    if (TotXssCap > 0 && tot_bm_vol > 0){ 		
		
        for (int i = 0; i < 4; i++) { 
            //compute total transported bed material (eq 20)
            double sbm_i(0.); 
            if(TotXssCap < tot_bm_vol) sbm_i = TotXssCap * BmVol[i] / tot_bm_vol; 
                           
            else sbm_i = BmVol[i];
                      
            for(int n=0;n<9 ;n++){ 
              //if(n!=4){  
              if(n!=4 && cell->GetFlows()[n]>0.){ //QS: I think there should be a check on whether there is a flow
                 
                Sbm_i_flow[i][n] = sbm_i * skr_n[n] / SKR; 
                
                Sbm_n[n] += Sbm_i_flow[i][n];
                         
                Stot[n] += Sbm_i_flow[i][n];
            
                Sbm += Sbm_i_flow[i][n]; 
      
              }
            } 
            
        } //end size fractions loop
      
    } //end if
    
    
    
	//4. compute transport of parent material (erosion)
	//update total excess capacity totXSScap
	
    vector<double> Sero_n(9, 0.); //transported eroded material (in all fractions) to neighbour cells 
    vector<vector<double>> Sero_i_flow(4, vector<double>(9,0.)); 
	 
    double Remaincap = max(0.0,(TotXssCap - Sbm));
    	
    if(Remaincap>0) { 
        
        for (int i = 0; i < 4; i++) { 
            double sero_i = Remaincap * erosion_para_[soiltype - 1][i+1]; //eqn 21 
            
                            
            for(int n=0;n<9;n++){ 
            
              //if(n!=4) {
              if(n!=4 && cell->GetFlows()[n]>0.){ //QS: I think there should be a check on whether there is a flow    
                              
                Sero_i_flow[i][n] = sero_i * skr_n[n] / SKR;
                
                Sero_n[n] += Sero_i_flow[i][n]; 
                    
                Stot[n] += Sero_i_flow[i][n];
                
                Sero += Sero_i_flow[i][n]; 
           
              }
            } //end neighbours loop
           
		} //end fractions loop                                                       
	} //end if                                                                          
 //if(Sero>0) cout<<"Sero = "<<Sero<<endl;
         
	//computing total sediment outflows for central cell (n=4)/ to each neighbour in each fraction
	
    for (int n = 0; n < 9 ; n++) {
      if(n!=4){
        Stot[4] -= Stot[n]; 
        for (int i = 0; i < 4; i++) {     
	        Ssus_i_flow[i][4] -= Ssus_i_flow[i][n];
			Sbm_i_flow[i][4] -= Sbm_i_flow[i][n];
		    Sero_i_flow[i][4] -= Sero_i_flow[i][n];
            //if(Sero_i_flow[i][n]!=0.) cout<<"Sero_i_flow[i][n] = "<<Sero_i_flow[i][n]<<" Sero_i_flow[i][4] = "<<Sero_i_flow[i][4]<<" i = "<<i<<" n = "<<n<<endl;
      }
		}
	}
      
	cell->SettotalSedimentTransport(Stot);   //totaltransportedmaterial
    //addd
    
   // if(Stot[4]==0) cout<<"Stot = "<<endl;
	cell->SetSsus_i_flow(Ssus_i_flow);
	cell->SetSbm_i_flow(Sbm_i_flow);
	cell->SetSero_i_flow(Sero_i_flow);
    
    
	double transportedsediment(0.);
    if(y_c==ny_) transportedsediment+=(Stot[6]+Stot[7]+Stot[8]); //m3
        //transportedsediment *= (p*1.0e6); //m3 -> g //added      
        //transportedsediment+=(Stot[6]+Stot[7]+Stot[8]); //g/min/cm2 
    
     
    
   // transportedsediment *= (time_increment_/60.); //g/cm2/min -> g/cm2 //added
    
    return transportedsediment;


}
  

		
//compute accumulate transported sediment (m3)
void RunCA::AccumulateSedimentTransport(Cell* cell) {

	int y = cell->GetY();
	int x = cell->GetX();

	//neighbourhood
	//(y-1, x-1)  (y-1, x)  (y-1, x+1)
	//( y , x-1)  ( y , x)  ( y , x+1)
	//(y+1, x-1)  (y+1, x)  (y+1, x+1)

	//flows
	//0 1 2
	//3 4 5
	//6 7 8

   
    //total accumualted sediment amount
    double com_accum_totsed = cells_[y + 1][x + 1]->GettotalSedimentTransport()[0] //down right neighbour, up left flow 
			                + cells_[y + 1][x]->GettotalSedimentTransport()[1] //down neighbour, up flow
			                + cells_[y + 1][x - 1]->GettotalSedimentTransport()[2] //down left neighbour, up right flow
			                + cells_[y][x + 1]->GettotalSedimentTransport()[3] //right neighbour, left flow
			                + cells_[y][x]->GettotalSedimentTransport()[4] //central cell, outflow
			                + cells_[y][x - 1]->GettotalSedimentTransport()[5] //left neighbour, right flow
			                + cells_[y - 1][x + 1]->GettotalSedimentTransport()[6] //up right neighbour, down left flow
			                + cells_[y - 1][x]->GettotalSedimentTransport()[7] //up neighbour, down flow
			                + cells_[y - 1][x - 1]->GettotalSedimentTransport()[8]; //up left neighbour, down right flow
    

	//accumulate sediment flow for each fraction size
    //vector<double> com_accum_sus_sedi(4, 0.), com_accum_bm_sedi(4, 0.), com_accum_ero_sedi(4, 0.); 

    vector<double> sus_sedi_inflow(4, 0.), bm_sedi_inflow(4, 0.), ero_sedi_inflow(4, 0.); //Qi added
	vector<double> sus_sedi_outflow(4, 0.), bm_sedi_outflow(4, 0.), ero_sedi_outflow(4, 0.); //Qi added
	
    for (int i = 0; i < 4; i++) { //i: use a single loop over fractions   
        ///// calculate accumulated suspended sedimenttransport for each fraction size ()
		//com_accum_sus_sedi[i] = cells_[y + 1][x + 1]->GetSsus_i_flow()[i][0] //down right neighbour, up left flow 
        sus_sedi_inflow[i]  = cells_[y + 1][x + 1]->GetSsus_i_flow()[i][0] //down right neighbour, up left flow //Qi added
				            + cells_[y + 1][x]->GetSsus_i_flow()[i][1] //down neighbour, up flow
				            + cells_[y + 1][x - 1]->GetSsus_i_flow()[i][2] //down left neighbour, up right flow
				            + cells_[y][x + 1]->GetSsus_i_flow()[i][3] //right neighbour, left flow
				            //+ cells_[y][x]->GetSsus_i_flow()[i][4] //central cell, outflow
				            + cells_[y][x - 1]->GetSsus_i_flow()[i][5] //left neighbour, right flow
				            + cells_[y - 1][x + 1]->GetSsus_i_flow()[i][6] //up right neighbour, down left flow
				            + cells_[y - 1][x]->GetSsus_i_flow()[i][7] //up neighbour, down flow
				            + cells_[y - 1][x - 1]->GetSsus_i_flow()[i][8];

        sus_sedi_outflow[i] = cells_[y][x]->GetSsus_i_flow()[i][4]; //central cell, outflow //Qi added                   
        /*
        //Qi added
        if(x==5) {
            if(cells_[y + 1][x + 1]->GetSsus_i_flow()[i][0] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", down right neighbour, up left suspended flow = "<<cells_[y + 1][x + 1]->GetSsus_i_flow()[i][0]<<endl;
            if(cells_[y + 1][x]->GetSsus_i_flow()[i][1] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", down neighbour, up suspended flow = "<<cells_[y + 1][x]->GetSsus_i_flow()[i][1]<<endl;
            if(cells_[y + 1][x - 1]->GetSsus_i_flow()[i][2] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", down left neighbour, up right suspended flow = "<<cells_[y + 1][x - 1]->GetSsus_i_flow()[i][2]<<endl;
            if(cells_[y][x + 1]->GetSsus_i_flow()[i][3] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", right neighbour, left suspended flow = "<<cells_[y][x + 1]->GetSsus_i_flow()[i][3]<<endl;
            if(cells_[y][x]->GetSsus_i_flow()[i][4] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", central cell, suspended outflow = "<<cells_[y][x]->GetSsus_i_flow()[i][4]<<endl;
            if(cells_[y][x - 1]->GetSsus_i_flow()[i][5] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", left neighbour, right suspended flow = "<<cells_[y][x - 1]->GetSsus_i_flow()[i][5]<<endl;
            if(cells_[y - 1][x + 1]->GetSsus_i_flow()[i][6] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", up right neighbour, down left suspended flow = "<<cells_[y - 1][x + 1]->GetSsus_i_flow()[i][6]<<endl;
            if(cells_[y - 1][x]->GetSsus_i_flow()[i][7] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", up neighbour, down suspended flow = "<<cells_[y - 1][x]->GetSsus_i_flow()[i][7]<<endl;
            if(cells_[y - 1][x - 1]->GetSsus_i_flow()[i][8] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", up left neighbour, down right suspended flow = "<<cells_[y - 1][x - 1]->GetSsus_i_flow()[i][8]<<endl;
        }  
        */       
        /*
        if(isnan(com_accum_sus_sedi[i])) {
        //if(isnan(com_accum_ero_sedi[i]) || com_accum_ero_sedi[i]>0.) {
            cout<<"i = "<<i<<" cells_[y + 1][x + 1]->GetSsus_i_flow()[i][0] = "<<cells_[y + 1][x + 1]->GetSsus_i_flow()[i][0]<<endl;
            cout<<"i = "<<i<<" cells_[y + 1][x]->GetSsus_i_flow()[i][1] = "<<cells_[y + 1][x]->GetSsus_i_flow()[i][1]<<endl;
            cout<<"i = "<<i<<" cells_[y + 1][x - 1]->GetSsus_i_flow()[i][2] = "<<cells_[y + 1][x - 1]->GetSsus_i_flow()[i][2]<<endl;
            cout<<"i = "<<i<<" cells_[y][x + 1]->GetSsus_i_flow()[i][3] = "<<cells_[y][x + 1]->GetSsus_i_flow()[i][3]<<endl;
            cout<<"i = "<<i<<" cells_[y][x]->GetSsus_i_flow()[i][4] = "<<cells_[y][x]->GetSsus_i_flow()[i][4]<<endl;
            cout<<"i = "<<i<<" cells_[y][x - 1]->GetSsus_i_flow()[i][5] = "<<cells_[y][x - 1]->GetSsus_i_flow()[i][5]<<endl;
            cout<<"i = "<<i<<" cells_[y - 1][x + 1]->GetSsus_i_flow()[i][6] = "<<cells_[y - 1][x + 1]->GetSsus_i_flow()[i][6]<<endl;
            cout<<"i = "<<i<<" cells_[y - 1][x]->GetSsus_i_flow()[i][7] = "<<cells_[y - 1][x]->GetSsus_i_flow()[i][7]<<endl;
            cout<<"i = "<<i<<" cells_[y - 1][x - 1]->GetSsus_i_flow()[i][8] = "<<cells_[y - 1][x - 1]->GetSsus_i_flow()[i][8]<<endl;
        }                            
	    */
		//////// calculate accumulated bed load sedimenttransport for each fraction size ()
		
		//com_accum_bm_sedi[i] = cells_[y + 1][x + 1]->GetSbm_i_flow()[i][0] //down right neighbour, up left flow 
        bm_sedi_inflow[i]    = cells_[y + 1][x + 1]->GetSbm_i_flow()[i][0] //down right neighbour, up left flow //Qi added
				             + cells_[y + 1][x]->GetSbm_i_flow()[i][1] //down neighbour, up flow
				             + cells_[y + 1][x - 1]->GetSbm_i_flow()[i][2] //down left neighbour, up right flow
				             + cells_[y][x + 1]->GetSbm_i_flow()[i][3] //right neighbour, left flow
				             //+ cells_[y][x]->GetSbm_i_flow()[i][4] //central cell, outflow
				             + cells_[y][x - 1]->GetSbm_i_flow()[i][5] //left neighbour, right flow
				             + cells_[y - 1][x + 1]->GetSbm_i_flow()[i][6] //up right neighbour, down left flow
				             + cells_[y - 1][x]->GetSbm_i_flow()[i][7] //up neighbour, down flow
				             + cells_[y - 1][x - 1]->GetSbm_i_flow()[i][8];
        
        bm_sedi_outflow[i]   = cells_[y][x]->GetSbm_i_flow()[i][4]; //central cell, outflow //Qi added
		/*
        //Qi added
        if(x==5) {
            if(cells_[y + 1][x + 1]->GetSbm_i_flow()[i][0] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", down right neighbour, up left bed material flow = "<<cells_[y + 1][x + 1]->GetSbm_i_flow()[i][0]<<endl;
            if(cells_[y + 1][x]->GetSbm_i_flow()[i][1] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", down neighbour, up bed material flow = "<<cells_[y + 1][x]->GetSbm_i_flow()[i][1]<<endl;
            if(cells_[y + 1][x - 1]->GetSbm_i_flow()[i][2] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", down left neighbour, up right bed material flow = "<<cells_[y + 1][x - 1]->GetSbm_i_flow()[i][2]<<endl;
            if(cells_[y][x + 1]->GetSbm_i_flow()[i][3] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", right neighbour, left bed material flow = "<<cells_[y][x + 1]->GetSbm_i_flow()[i][3]<<endl;
            if(cells_[y][x]->GetSbm_i_flow()[i][4] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", central cell, bed material outflow = "<<cells_[y][x]->GetSbm_i_flow()[i][4]<<endl;
            if(cells_[y][x - 1]->GetSbm_i_flow()[i][5] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", left neighbour, right bed material flow = "<<cells_[y][x - 1]->GetSbm_i_flow()[i][5]<<endl;
            if(cells_[y - 1][x + 1]->GetSbm_i_flow()[i][6] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", up right neighbour, down left bed material flow = "<<cells_[y - 1][x + 1]->GetSbm_i_flow()[i][6]<<endl;
            if(cells_[y - 1][x]->GetSbm_i_flow()[i][7] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", up neighbour, down bed material flow = "<<cells_[y - 1][x]->GetSbm_i_flow()[i][7]<<endl;
            if(cells_[y - 1][x - 1]->GetSbm_i_flow()[i][8] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", up left neighbour, down right bed material flow = "<<cells_[y - 1][x - 1]->GetSbm_i_flow()[i][8]<<endl;
        } 
        */		
        
        //////// calculate accumulated eroded sedimenttransport for each fraction size ()
		
		//com_accum_ero_sedi[i] = cells_[y + 1][x + 1]->GetSero_i_flow()[i][0] //down right neighbour, up left flow 
        ero_sedi_inflow[i]      = cells_[y + 1][x + 1]->GetSero_i_flow()[i][0] //down right neighbour, up left flow //Qi added
		                 		+ cells_[y + 1][x]->GetSero_i_flow()[i][1] //down neighbour, up flow
				                + cells_[y + 1][x - 1]->GetSero_i_flow()[i][2] //down left neighbour, up right flow
				                + cells_[y][x + 1]->GetSero_i_flow()[i][3] //right neighbour, left flow
				                //+ cells_[y][x]->GetSero_i_flow()[i][4] //central cell, outflow
				                + cells_[y][x - 1]->GetSero_i_flow()[i][5] //left neighbour, right flow
				                + cells_[y - 1][x + 1]->GetSero_i_flow()[i][6] //up right neighbour, down left flow
				                + cells_[y - 1][x]->GetSero_i_flow()[i][7] //up neighbour, down flow
				                + cells_[y - 1][x - 1]->GetSero_i_flow()[i][8];

        ero_sedi_outflow[i]     = cells_[y][x]->GetSero_i_flow()[i][4]; //central cell, outflow //Qi added
        /*
        //Qi added
        if(x==5) {
            if(cells_[y + 1][x + 1]->GetSero_i_flow()[i][0] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", down right neighbour, up left eroded material flow = "<<cells_[y + 1][x + 1]->GetSero_i_flow()[i][0]<<endl;
            if(cells_[y + 1][x]->GetSero_i_flow()[i][1] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", down neighbour, up eroded material flow = "<<cells_[y + 1][x]->GetSero_i_flow()[i][1]<<endl;
            if(cells_[y + 1][x - 1]->GetSero_i_flow()[i][2] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", down left neighbour, up right eroded material flow = "<<cells_[y + 1][x - 1]->GetSero_i_flow()[i][2]<<endl;
            if(cells_[y][x + 1]->GetSero_i_flow()[i][3] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", right neighbour, left eroded material flow = "<<cells_[y][x + 1]->GetSero_i_flow()[i][3]<<endl;
            if(cells_[y][x]->GetSero_i_flow()[i][4] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", central cell, eroded material outflow = "<<cells_[y][x]->GetSero_i_flow()[i][4]<<endl;
            if(cells_[y][x - 1]->GetSero_i_flow()[i][5] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", left neighbour, right eroded material flow = "<<cells_[y][x - 1]->GetSero_i_flow()[i][5]<<endl;
            if(cells_[y - 1][x + 1]->GetSero_i_flow()[i][6] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", up right neighbour, down left eroded material flow = "<<cells_[y - 1][x + 1]->GetSero_i_flow()[i][6]<<endl;
            if(cells_[y - 1][x]->GetSero_i_flow()[i][7] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", up neighbour, down eroded material flow = "<<cells_[y - 1][x]->GetSero_i_flow()[i][7]<<endl;
            if(cells_[y - 1][x - 1]->GetSero_i_flow()[i][8] != 0.)
                cout<<"i = "<<i<<", y = "<<y<<", up left neighbour, down right eroded material flow = "<<cells_[y - 1][x - 1]->GetSero_i_flow()[i][8]<<endl;
        } 
        */

        /*
        //if(isnan(com_accum_ero_sedi[i])) {
        if(com_accum_ero_sedi[i]<0.) {
        //if(isnan(com_accum_ero_sedi[i]) || com_accum_ero_sedi[i]>0.) {
            cout<<"i = "<<i<<" cells_[y + 1][x + 1]->GetSero_i_flow()[i][0] = "<<cells_[y + 1][x + 1]->GetSero_i_flow()[i][0]<<endl;
            cout<<"i = "<<i<<" cells_[y + 1][x]->GetSero_i_flow()[i][1] = "<<cells_[y + 1][x]->GetSero_i_flow()[i][1]<<endl;
            cout<<"i = "<<i<<" cells_[y + 1][x - 1]->GetSero_i_flow()[i][2] = "<<cells_[y + 1][x - 1]->GetSero_i_flow()[i][2]<<endl;
            cout<<"i = "<<i<<" cells_[y][x + 1]->GetSero_i_flow()[i][3] = "<<cells_[y][x + 1]->GetSero_i_flow()[i][3]<<endl;
            cout<<"i = "<<i<<" cells_[y][x]->GetSero_i_flow()[i][4] = "<<cells_[y][x]->GetSero_i_flow()[i][4]<<endl;
            cout<<"i = "<<i<<" cells_[y][x - 1]->GetSero_i_flow()[i][5] = "<<cells_[y][x - 1]->GetSero_i_flow()[i][5]<<endl;
            cout<<"i = "<<i<<" cells_[y - 1][x + 1]->GetSero_i_flow()[i][6] = "<<cells_[y - 1][x + 1]->GetSero_i_flow()[i][6]<<endl;
            cout<<"i = "<<i<<" cells_[y - 1][x]->GetSero_i_flow()[i][7] = "<<cells_[y - 1][x]->GetSero_i_flow()[i][7]<<endl;
            cout<<"i = "<<i<<" cells_[y - 1][x - 1]->GetSero_i_flow()[i][8] = "<<cells_[y - 1][x - 1]->GetSero_i_flow()[i][8]<<endl;
        }
        */
	    
    }

    /*   
	cell->Setaccum_totsed(com_accum_totsed);//updated sediment transport(m)
	//cell->Setaccum_bm_sedi(com_accum_sus_sed; //QS: this is wrong
	//cell->Setaccum_sus_sedi(com_accum_bm_sedi); //QS: this is wrong
	cell->Setaccum_bm_sedi(com_accum_bm_sedi); //QS: fixed
	cell->Setaccum_sus_sedi(com_accum_sus_sedi); //QS: fixed   
	cell->Setaccum_ero_sedi(com_accum_ero_sedi);
    */

    //Qi added
    cell->Setaccum_totsed(com_accum_totsed);//updated sediment transport(m)
    cell->SetSusSediInflow(sus_sedi_inflow);
    cell->SetSusSediOutflow(sus_sedi_outflow);
    cell->SetBmSediInflow(bm_sedi_inflow);
    cell->SetBmSediOutflow(bm_sedi_outflow);
    cell->SetEroSediInflow(ero_sedi_inflow);
    cell->SetEroSediOutflow(ero_sedi_outflow);
       	
}

//
//Function ComputeSedimentTransport for bedlad(cell)


double RunCA::ComputeBedloadTransport(Cell* cell) {
     vector<double> skr_nb(9, 0.); //volumes of outgoing sediment from different fractions - i ( different nighbour)
    
    //###########################3
	
    double Sbmb=0,Serob=0; //tota outflowing suspended sediment, bed material and eroded material
    double G=2.65; // unit (gr/cm3)
    //double g = 9.8;    //(m/s^2)-> 35280(m/min^2)   ->3531600(cm/min^2)
    double g = 3531600.;
	//retrieve erosion parameters
	int soiltype = cell->GetSoilType();
	vector<double> ero_params = erosion_para_[soiltype - 1];
	double clay = ero_params[1];
	double silt = ero_params[2];
	double sand = ero_params[3];
	double gravel = ero_params[4];
	double cf = ero_params[5];
	double cr = ero_params[6];
	double p = ero_params[7];//solid density
    double y = ero_params[8];//i ??density of water constan
 
    double o= ero_params[9];  //rill erodibility(kg s/m-4)
    double c= ero_params[10];  //interrill erodibility (s m-1)chien_coefficient
    double pr= ero_params[11];//pa
   
	//double ws = ero_params[8]; //i ??
  
	//find neighbour cells
	vector<Cell*> neighbours;
	int y_c = cell->GetY();
	int x_c = cell->GetX();
	for (int y = y_c - 1; y <= y_c + 1; y++) {
		for (int x = x_c - 1; x <= x_c + 1; x++) {
			Cell* neighbour = cells_[y][x];
			neighbours.push_back(neighbour);
		}
	}

    
    vector<double> Stotb(9, 0.); //total sediment (from all sources and all fractions) transported to neighbours - i
    vector<double>SKR_i(4,0.);
	//1. compute total volume of outgoing sediment SKR
	//for each outgoing cell n :
	double SKRb = 0; //total volume of outgoing sediment
    
    int count1(0); 
    double sum_s(0.);
    double sum_d(0.);
    double mean_s(0.); //mean slope
    double mean_d(0.); //mean distance
    double CE= cell->GetLandElevation();// elevation of central cell
    vector<double> slope = cell->GetSlope();
   // double SN = cell-> GetSlope();
     for (int n=0;n<9;n++) {
        if(neighbours[n]->GetValid()) {
             //vector<double> flow = cell->GetFlows();
            //double U = flow [n];
            double EN = neighbours[n]->GetLandElevation(); // elevation of neighbour cell
            double SN =slope[n];
            
            if(EN < CE) {
                
                //sum_s+=EN; Qi
                count1++;
                double distance = cell_size_; 
                if(n==0 || n==2 || n==6 || n==8) distance *= 1.414; ///digonal neighbours
                sum_d += distance;
                //vector<double> SN(9,0.); 
                // SN[n] = fabs(CE-EN)/(distance); //slope  Shahla added
                sum_s+=SN;                           
            }
        } 
    } 
    if(count1>0) {
        mean_s = sum_s / count1;
        mean_d = sum_d / count1;
    } else {
        return 0.;
    }
	
   // cout<<"passed 4.1"<<endl;

	double water_depth = cell->GetWaterDepth();
    double h = cell->GetWaterDepth()/10.;//Q: mm->cm
    
    //particlesize
    vector<double> d(4, 0.);
    d[0] = 0.001/10.; //clay mm->cm 0.0001        Ycr=0  
    d[1] = 0.016/10.; //silt mm->cm  0.0016       Ycr=10^-1 or 0.1    -> 0.001   (gr/cm2)
    d[2] = 0.35/10.; //sand  mm->cm  0.035        Ycr=0.15 (N/M^2)   -> 0.0015  (gr/cm2)
    d[3] = 3/10.; //graval mm->cm  0.3            Ycr=2 (N/M^2)      -> 0.0203 (gr/cm2)
    // critical shear strees
    vector<double> YCR(4, 0.);
    /*
    YCR[0]=0.;
    YCR[1]=0.001;
    YCR[2]= 0.000000000015 ;        //0.0015;
    YCR[3]=0.000000000203; 
    */
   
   

    double Shear_stres = y *h *mean_s; //shear stress
    double Shear_velocity = pow(( g * mean_s * h),0.5); //shear velocity
    

    
    //vector<double> beta (4, 0.);   

    vector<double> Sbm_nb(9, 0.); //transported bed material (in all fractions) to neighbour cells 
    //bedload sediment flow to each neighbour in each fraction
    vector<vector<double>> Sbm_i_flowb(4, vector<double>(9,0.));  
    vector<double> BmVolb = cell->GetBedMaterialVolumesb(); //bed material volumes in different fractions    
   // cout<<"passed 4.2"<<endl;
    vector<double> SKR_left (4, 0.); //transport capacity for each of size fractions
    //boys equation
    double Y=0.54/(p -y)*g;  // or Y=0.173/d^3/4
	

    for (int i = 0; i < 4; i++) {  
        //if(i=4) if gravel
        if(i>1 ){
           // YCR[i]=((p/y)-1)*d[i]/13; // computing critical shear stress with kreys formula
           //m YCR[i]=55.7 * ((p/y)-1.)*o*d[i];  //computing critical shear_ stress  by Akiand sato,s formul , o is coifficent depend on size
             //YCR[i]= ((p/y)-1.)*(d[i]/13);  //computing critical shear_ stress  by Krey’s Formula 
             
               //YCR[i]= ((100*((p/y)-1.)*pow(d[i],(6/5)))/3)*((2+o)/ (1+(2*o)));     // vector<double> YCR(4, 0.);
           
           // YCR[i]= y* U * c^2/(y-p)*d[i];
            //YCR[i]= 16.8*(pow(d[i],1.262));
            YCR[i]=3.54*(pow(10,(-28.1*o)));  //Sg = specific gravity of soil=2.65-2.80,d50=>o

            double Dimensionles_ShearStress =( (Shear_stres/y)/((G-1) * g * d[i])); //Dimensionless shear stress
           // vector<double> Y (4, 0.);
         
           
           // cout<< " num2= " << num1;
         
          // SKR_i[i]=Y*Shear_stres*(Shear_stres-YCR[i]); // Du Boys formul

           //SKR_i[i]=10 *((Shear_stres-YCR[i])/(p-y)*d[i])*discharge*SLOPE  //shields formula

           SKR_i[i]= c * pow((Shear_stres-YCR[i]),(3/2));  //Chein (1956) formula
           // SKR_i[i] =1.28*c*(sqrt(Shear_stres)-0.7*sqrt(YCR[i])) *  (Shear_stres-YCR[i]) ; //Hunziker [9],geosciences-10-00368 no good

           // SKR_i[i]=c* pow(Shear_stres,1.5)*(1-(sqrt(YCR[i]/Shear_stres)))*((1-YCR[i])/Shear_stres);

           //################# Schoklitch (1949)
           //qc[i]= (0.00001944 *d[i])/slope;
            // SKR_i[i]=7000*((pow(slope,(3/2)))/pow(d[i],0.5))*(discharge - qc[i]);   or SKR_i[i]=2500*pow(slope,(3/2))*(discharge - qc[i])


  
          //######################



           

          
           //if (i>=2) cout<< " SKR_i[i]= " << SKR_i[i]<<endl;
            SKR_i[i] *= ( (mean_d*100.)* (time_increment_/60.) / p );//from gr/min/cm ->cm3
            SKR_i[i]/= 1.0e6; //m3
            SKRb += SKR_i[i];
            
           
            
           // cout<< " SKRb= " << SKRb<<endl;
          // cout<<"passed 4.1"<<endl;
           // if (SKRb == 0) cout<<"passed SKRb > 0"<<endl;

            //compute bedload transport from deposition part
            
            //vector<double> VOL (4, 0.);
            //VOL[i] = std::min(SKRb, BmVolb[i]);
             double VOL = std::min(SKRb, BmVolb[i]);
            SKR_left[i] = SKRb - VOL;
            

            for (int n = 0; n < 9; n++) {
                 if (cell->GetValid() ) {
                    
                    if(n!=4){
                        double EN = neighbours[n]->GetLandElevation(); // elevation of neighbour cell
                        if(EN < CE) {
                            
                                 double distance = cell_size_; 
                                if(n==0 || n==2 || n==6 || n==8) {distance *= 1.414;} ///digonal neighbours
                                    sum_d += distance;
                                    //vector<double> SN(9,0.); 
                                    //SN[n] = fabs(CE-EN)/(distance);
                            
                                     double SN = slope[n];
                                     Sbm_i_flowb[i][n] = VOL * SN / sum_s; 
                                   // Sbm_i_flowb[i][n] = VOL[i] * cell-> GetSlope()[n] / sum_s; 

                                    //Sbm_i_flowb[i][n] = VOL[i] * neighbours[n]-> GetSlope() / sum_s; Qi
                                    Sbm_nb[n] += Sbm_i_flowb[i][n];
                                    Stotb[n] += Sbm_i_flowb[i][n];
                                        Sbmb += Sbm_i_flowb[i][n];   
                        }
                    }

                }
            }
        }
    }
    //cout<<"passed 4.2"<<endl;
	//4. compute transport of parent material (erosion)
	//update total excess capacity totXSScap
    double Remaincap(0.);
    for (int i = 0; i < 4; i++) { 
        Remaincap += SKR_left[i];
    }

    vector<double> Sero_nb(9, 0.); //transported eroded material (in all fractions) to neighbour cells 
    vector<vector<double>> Sero_i_flowb(4, vector<double>(9,0.)); 
    if(Remaincap>0) { 		    	
        for (int i = 2; i < 4; i++) { 
            for (int n = 0; n < 9; n++) {
                if(n!=4){
                    double EN = neighbours[n]->GetLandElevation(); // elevation of neighbour cell
                    if(EN < CE) {
                        //
                                    double distance = cell_size_; 
                                    if(n==0 || n==2 || n==6 || n==8) distance *= 1.414; ///digonal neighbours
                                    sum_d += distance;
                                   // vector<double> SN(9,0.); 
                                    //SN[n] = fabs(CE-EN)/(distance);
                                    double SN = slope[n];
                        //
                                    Sero_i_flowb[i][n] = SKR_left[i] * SN / sum_s;
                                    //Sero_i_flowb[i][n] = SKR_left[i] * cell->GetSlope()[n] / sum_s;
                                   //Sero_i_flowb[i][n] = SKR_left[i] * neighbours[n]->GetSlope() / sum_s;  Qi
                                    Sero_nb[n] += Sero_i_flowb[i][n];
                                    Stotb[n] += Sero_i_flowb[i][n];
                                    Serob += Sero_i_flowb[i][n];   
                    }
                }

            }
        }
    }
    //cout<<"passed 4.3"<<endl;
	//computing total sediment outflows for central cell (n=4)/ to each neighbour in each fraction
    for (int n = 0; n < 9 ; n++) {
      if(n!=4){
        Stotb[4] -= Stotb[n]; 
        for (int i = 0; i < 4; i++) {     
			Sbm_i_flowb[i][4] -= Sbm_i_flowb[i][n];
		    Sero_i_flowb[i][4] -= Sero_i_flowb[i][n];
            //if(Sero_i_flow[i][n]!=0.) cout<<"Sero_i_flow[i][n] = "<<Sero_i_flow[i][n]<<" Sero_i_flow[i][4] = "<<Sero_i_flow[i][4]<<" i = "<<i<<" n = "<<n<<endl;
        }
	  }
	}
      
	cell->SettotalSedimentTransportb(Stotb);   //totaltransportedmaterial
    //addd
    cell->SetSbm_i_flowb(Sbm_i_flowb);
	cell->SetSero_i_flowb(Sero_i_flowb);
	//cell->SetSbm_i_bedload(Sbm_i_flowb);
	//cell->SetSero_i_bedload(Sero_i_flowb);
    
    
	double transportedsedimentb(0.);
    if(y_c==ny_) transportedsedimentb+=(Stotb[6]+Stotb[7]+Stotb[8]); //m3   
   // cout<<"passed 4.4"<<endl;
    return transportedsedimentb;
}

 




///....compute Accumulate bedload transport
void RunCA::AccumulateBedloadTransport(Cell* cell) {

	int y = cell->GetY();
	int x = cell->GetX();

	//neighbourhood
	//(y-1, x-1)  (y-1, x)  (y-1, x+1)
	//( y , x-1)  ( y , x)  ( y , x+1)
	//(y+1, x-1)  (y+1, x)  (y+1, x+1)

	//flows
	//0 1 2
	//3 4 5
	//6 7 8

   
    //total accumualted sediment amount
    double com_accum_totsedb = cells_[y + 1][x + 1]->GettotalSedimentTransportb()[0] //down right neighbour, up left flow 
			                + cells_[y + 1][x]->GettotalSedimentTransportb()[1] //down neighbour, up flow
			                + cells_[y + 1][x - 1]->GettotalSedimentTransportb()[2] //down left neighbour, up right flow
			                + cells_[y][x + 1]->GettotalSedimentTransportb()[3] //right neighbour, left flow
			                + cells_[y][x]->GettotalSedimentTransportb()[4] //central cell, outflow
			                + cells_[y][x - 1]->GettotalSedimentTransportb()[5] //left neighbour, right flow
			                + cells_[y - 1][x + 1]->GettotalSedimentTransportb()[6] //up right neighbour, down left flow
			                + cells_[y - 1][x]->GettotalSedimentTransportb()[7] //up neighbour, down flow
			                + cells_[y - 1][x - 1]->GettotalSedimentTransportb()[8]; //up left neighbour, down right flow
    

	//accumulate sediment flow for each fraction size
    //vector<double> com_accum_sus_sedi(4, 0.), com_accum_bm_sedi(4, 0.), com_accum_ero_sedi(4, 0.); 

      vector<double>  bm_sedi_inflowb(4, 0.), ero_sedi_inflowb(4, 0.); //Qi added               sus_sedi_inflow(4, 0.),
	  vector<double>  bm_sedi_outflowb(4, 0.), ero_sedi_outflowb(4, 0.); //Qi added             sus_sedi_outflow(4, 0.),
	
    for (int i = 0; i < 4; i++) { //i: use a single loop over fractions   
                 
		//////// calculate accumulated bed load sedimenttransport for each fraction size ()
		
		//com_accum_bm_sedi[i] = cells_[y + 1][x + 1]->GetSbm_i_flow()[i][0] //down right neighbour, up left flow 
        bm_sedi_inflowb[i]    = cells_[y + 1][x + 1]->GetSbm_i_flowb()[i][0] //down right neighbour, up left flow //Qi added
				             + cells_[y + 1][x]->GetSbm_i_flowb()[i][1] //down neighbour, up flow
				             + cells_[y + 1][x - 1]->GetSbm_i_flowb()[i][2] //down left neighbour, up right flow
				             + cells_[y][x + 1]->GetSbm_i_flowb()[i][3] //right neighbour, left flow
				             //+ cells_[y][x]->GetSbm_i_flow()[i][4] //central cell, outflow
				             + cells_[y][x - 1]->GetSbm_i_flowb()[i][5] //left neighbour, right flow
				             + cells_[y - 1][x + 1]->GetSbm_i_flowb()[i][6] //up right neighbour, down left flow
				             + cells_[y - 1][x]->GetSbm_i_flowb()[i][7] //up neighbour, down flow
				             + cells_[y - 1][x - 1]->GetSbm_i_flowb()[i][8];
        
        bm_sedi_outflowb[i]   = cells_[y][x]->GetSbm_i_flowb()[i][4]; //central cell, outflow //Qi added
			
        
        //////// calculate accumulated eroded sedimenttransport for each fraction size ()
		
		//com_accum_ero_sedi[i] = cells_[y + 1][x + 1]->GetSero_i_flow()[i][0] //down right neighbour, up left flow 
        ero_sedi_inflowb[i]      = cells_[y + 1][x + 1]->GetSero_i_flowb()[i][0] //down right neighbour, up left flow //Qi added
		                 		+ cells_[y + 1][x]->GetSero_i_flowb()[i][1] //down neighbour, up flow
				                + cells_[y + 1][x - 1]->GetSero_i_flowb()[i][2] //down left neighbour, up right flow
				                + cells_[y][x + 1]->GetSero_i_flowb()[i][3] //right neighbour, left flow
				                //+ cells_[y][x]->GetSero_i_flow()[i][4] //central cell, outflow
				                + cells_[y][x - 1]->GetSero_i_flowb()[i][5] //left neighbour, right flow
				                + cells_[y - 1][x + 1]->GetSero_i_flowb()[i][6] //up right neighbour, down left flow
				                + cells_[y - 1][x]->GetSero_i_flowb()[i][7] //up neighbour, down flow
				                + cells_[y - 1][x - 1]->GetSero_i_flowb()[i][8];

        ero_sedi_outflowb[i]     = cells_[y][x]->GetSero_i_flowb()[i][4]; //central cell, outflow //Qi added
    
          //Qi added
          
        //if(x==5) {
            //if(cells_[y + 1][x + 1]->GetSero_i_flowb()[i][0] != 0.)
                //cout<<"i = "<<i<<", y = "<<y<<", down right neighbour, up left suspended flow = "<<cells_[y + 1][x + 1]->GetSero_i_flowb()[i][0]<<endl;
           // if(cells_[y + 1][x]->GetSero_i_flowb()[i][1] != 0.)
               // cout<<"i = "<<i<<", y = "<<y<<", down neighbour, up suspended flow = "<<cells_[y + 1][x]->GetSero_i_flowb()[i][1]<<endl;
            //if(cells_[y + 1][x - 1]->GetSero_i_flowb()[i][2] != 0.)
               // cout<<"i = "<<i<<", y = "<<y<<", down left neighbour, up right suspended flow = "<<cells_[y + 1][x - 1]->GetSero_i_flowb()[i][2]<<endl;
           // if(cells_[y][x + 1]->GetSero_i_flowb()[i][3] != 0.)
               // cout<<"i = "<<i<<", y = "<<y<<", right neighbour, left suspended flow = "<<cells_[y][x + 1]->GetSero_i_flowb()[i][3]<<endl;
           // if(cells_[y][x]->GetSero_i_flowb()[i][4] != 0.)
               // cout<<"i = "<<i<<", y = "<<y<<", central cell, suspended outflow = "<<cells_[y][x]->GetSero_i_flowb()[i][4]<<endl;
           // if(cells_[y][x - 1]->GetSero_i_flowb()[i][5] != 0.)
               // cout<<"i = "<<i<<", y = "<<y<<", left neighbour, right suspended flow = "<<cells_[y][x - 1]->GetSero_i_flowb()[i][5]<<endl;
            //if(cells_[y - 1][x + 1]->GetSero_i_flowb()[i][6] != 0.)
               // cout<<"i = "<<i<<", y = "<<y<<", up right neighbour, down left suspended flow = "<<cells_[y - 1][x + 1]->GetSero_i_flowb()[i][6]<<endl;
            //if(cells_[y - 1][x]->GetSero_i_flowb()[i][7] != 0.)
               // cout<<"i = "<<i<<", y = "<<y<<", up neighbour, down suspended flow = "<<cells_[y - 1][x]->GetSero_i_flowb()[i][7]<<endl;
           // if(cells_[y - 1][x - 1]->GetSero_i_flowb()[i][8] != 0.)
              //  cout<<"i = "<<i<<", y = "<<y<<", up left neighbour, down right suspended flow = "<<cells_[y - 1][x - 1]->GetSero_i_flowb()[i][8]<<endl;
       // } 
            
        
     
	    
    }

  
    //Qi added
    cell->Setaccum_totsedb(com_accum_totsedb);//updated sediment transport(m)
   
   
    cell->SetBmSediInflowb(bm_sedi_inflowb);
    cell->SetBmSediOutflowb(bm_sedi_outflowb);
    cell->SetEroSediInflowb(ero_sedi_inflowb);
    cell->SetEroSediOutflowb(ero_sedi_outflowb);
       	
}


//
double RunCA::ComputeGlobalTimeStepSize(bool initial_guess)
{
    vector<double> timesteps;     
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];              
            if(cell->GetValid()) {
                double ori_depth = cell->GetWaterDepth();
                if(initial_guess) cell->SetWaterDepth(1.);
                if(cell->GetWaterDepth() != 0.) { //there is water to flow
                    //find neighbour cells
                    vector<Cell*> neighbours;
                    int y_c = cell->GetY();
                    int x_c = cell->GetX();
                    for(int y=y_c-1;y<=y_c+1;y++){
                        for(int x=x_c-1;x<=x_c+1;x++){
                            Cell* neighbour = cells_[y][x];
                            neighbours.push_back(neighbour);
                        }
                    }    

                    int soiltype = cell->GetSoilType();
                    double manning_n = iparams_[soiltype-1][8];
                    double water_depth = cell->GetWaterDepth(); //mm
                    double wh_c = cell->GetWaterSurfaceElevation(); //surface water elevation of central cell 
                    double h_c = cell->GetLandElevation(); //added 
                    for (int n=0;n<9;n++) {
                        if(n!=4 && neighbours[n]->GetValid()) {
                            double wh_n = neighbours[n]->GetWaterSurfaceElevation(); //surface water elevation of neighbour cell
                            double h_n = neighbours[n]->GetLandElevation(); //added
                            if(wh_c > wh_n) {  
                                double distance = cell_size_;
                                if(n==0 || n==2 || n==6 || n==8) distance *= 1.414;
                                //double depth = max(water_depth, 1.0); //added
                                double depth = water_depth; //Qi added - 5/11
                                double time_step=distance*manning_n/pow(depth/1000.,0.667)/pow((wh_c-wh_n)/1000./distance,0.5);
                                timesteps.push_back(time_step);  
                            }         
                        }
                    }
                    if(initial_guess) cell->SetWaterDepth(ori_depth);
                }
            }
        }
    }	

    sort(timesteps.begin(), timesteps.end());
   int five_percent = timesteps.size() / 20.;
   // int five_percent = timesteps.size()/30. ;

    if(initial_guess) {
        string file_name; 
        file_name = "initial_time_steps.txt";
        ofstream ifs(file_name.c_str());
        if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}    
        for(auto n : timesteps) ifs << n << endl;
        ifs.close();
    }

    return timesteps[five_percent];
}


//***********************************************************************
double RunCA::ComputeTimeStepSize(Cell* cell)
{
    double min_time_step (numeric_limits<double>::max());

    if(cell->GetWaterDepth() == 0.) { //there is no water to flow
        return min_time_step;
    }

    //find neighbour cells
    vector<Cell*> neighbours;
    int y_c = cell->GetY();
    int x_c = cell->GetX();
    for(int y=y_c-1;y<=y_c+1;y++){
        for(int x=x_c-1;x<=x_c+1;x++){
            Cell* neighbour = cells_[y][x];
            neighbours.push_back(neighbour);
        }
    }    

    int soiltype = cell->GetSoilType();
    double manning_n = iparams_[soiltype-1][8];
    double water_depth = cell->GetWaterDepth(); //mm
    double wh_c = cell->GetWaterSurfaceElevation(); //surface water elevation of central cell 
    double h_c = cell->GetLandElevation(); //added 
    for (int n=0;n<9;n++) {
        if(n!=4 && neighbours[n]->GetValid()) {
        //if(n!=4 && neighbours[n]->GetValid() && neighbours[n]->GetY()!=ny_+1) { //Qi added
            double wh_n = neighbours[n]->GetWaterSurfaceElevation(); //surface water elevation of neighbour cell
            double h_n = neighbours[n]->GetLandElevation(); //added
            if(wh_c > wh_n) {
            //if(h_c > h_n) {    
                double distance = cell_size_;
                if(n==0 || n==2 || n==6 || n==8) distance *= 1.414;
                double depth = max(water_depth, 1.0); //added
                //double depth = water_depth; //Qi added - 5/11
                double time_step=distance*manning_n/pow(depth/1000.,0.667)/pow((wh_c-wh_n)/1000./distance,0.5);
                //double time_step=distance*manning_n/pow(depth/1000.,0.667)/pow((h_c-h_n)/distance,0.5); //added
                //double adjusted_n = manning_n * pow(water_depth/10., -1./3.); 
                //double time_step=distance*adjusted_n/pow(depth/1000.,0.667)/pow((h_c-h_n)/distance,0.5); //added 
                min_time_step = min(min_time_step, time_step);


              
                //shahla added sort timestep (31/05/2021).
               // vector<double> sort_time_step ;
               // sort_time_step.push_back(time_step);
                
                //if(sort_time_step[0]>0.) cout<<"sort_time_step = " <<sort_time_step[0]<<endl;
                //sort(sort_time_step.begin(),sort_time_step.end());
                //double si= sort_time_step.size();
                //if(si>0.) cout<<"si = " <<si<<endl;
                /*
                if(sort_time_step[0]>0.) cout<<"sort_time_step[0] = " <<sort_time_step[0]<<endl;
                if(sort_time_step[1]>0.) cout<<"sort_time_step[1] = " <<sort_time_step[1]<<endl;
                if(sort_time_step[2]>0.) cout<<"sort_time_step[2] = " <<sort_time_step[2]<<endl;
             */
               
                
           

            }
        }
    }

    return min_time_step;    
}




void RunCA::SaveInputDataToVti()
{
    string file_name(model_name_);
    int minute = round(model_time_/60.);
    string time = to_string(minute);
    file_name += ("-input_data.vti");
    ofstream out_file(file_name.c_str());

    out_file<<"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"<< endl;
    out_file<<"<ImageData WholeExtent=\"0 "<<nx_<<" 0 "<<ny_<<" 0 " <<0<<"\" Origin=\""<<0<<" "<<0<<" "<<0<<"\" Spacing=\""\
            <<cell_size_<<" "<<cell_size_<<" "<<cell_size_<<"\">" << endl;
    out_file << "<Piece Extent=\"0 "<<nx_<<" 0 "<<ny_<<" 0 "<<0<<"\">" << endl;
    out_file << "<CellData>" << endl;

    //land evelations
    out_file << "<DataArray type=\"Float32\" Name=\"land elevation(m)\" format=\"ascii\">" << endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            if(cell->GetValid()) out_file<<std::setprecision(10)<<cell->GetLandElevation()<<" ";
            else out_file<<"-9999"<<" ";
        }
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;
    
    //soil types
    out_file << "<DataArray type=\"Float32\" Name=\"soil types\" format=\"ascii\">" << endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            if(cell->GetValid()) out_file<<cell->GetSoilType()<<" ";
            else out_file<<"-9999"<<" ";
        }
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;
    /*
    //out_elevaion_change
    out_file<<"<DataArray type=\"Float32\" Name=\"output-Elevation_change(y)\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            if(cell->GetValid()) out_file << cell->GetLandElavationChange() << " ";
            else out_file << "nan" << " ";
        }
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;    
    //
    */
    out_file << "</CellData>" << endl;
    out_file << "</Piece>" << endl;
    out_file << "</ImageData>" << endl;
    out_file << "</VTKFile>" << endl;

    // close file
    out_file.close();    
}


void RunCA::SaveResultsToVti()
{
    string file_name(model_name_);
    int minute = round(model_time_/60.);
    string time = to_string(minute);
    file_name += ("-output-results-" + time + "min.vti");
    ofstream out_file(file_name.c_str());    

    out_file<<"<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">"<< endl;
    out_file<<"<ImageData WholeExtent=\"0 "<<nx_<<" 0 "<<ny_<<" 0 " <<0<<"\" Origin=\""<<0<<" "<<0<<" "<<0<<"\" Spacing=\""\
            <<cell_size_<<" "<<cell_size_<<" "<<cell_size_<<"\">" << endl;
    out_file << "<Piece Extent=\"0 "<<nx_<<" 0 "<<ny_<<" 0 "<<0<<"\">" << endl;
    out_file << "<CellData>" << endl;

    //water depth(mm)
    out_file << "<DataArray type=\"Float32\" Name=\"water depth(mm)\" format=\"ascii\">" << endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10) << cell->GetWaterDepth() << " ";
            if(cell->GetValid()) {
                double value = cell->GetWaterDepth();
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }            
            else out_file << "-9999" << " ";
        }
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;

    //water surface elevation(m)
    out_file << "<DataArray type=\"Float32\" Name=\"water surface elevation(m)\" format=\"ascii\">" << endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10) << cell->GetWaterSurfaceElevation()/1000. << " ";
            if(cell->GetValid()) {
                double value = cell->GetWaterSurfaceElevation()/1000.;
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        }
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;

    //infiltration rate(mm/h)
    out_file << "<DataArray type=\"Float32\" Name=\"infiltration rate(mm/h)\" format=\"ascii\">" << endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetInfiltrationRate() << " ";
            if(cell->GetValid()) {
                double value = cell->GetInfiltrationRate();
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        }
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;    

    //infiltration capacity(mm/h)
    out_file << "<DataArray type=\"Float32\" Name=\"infiltration capacity(mm/h)\" format=\"ascii\">" << endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetInfiltrationCapacity() << " ";
            if(cell->GetValid()) {
                double value = cell->GetInfiltrationCapacity();
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        }
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;      

    //drainage rate(mm/h)
    out_file << "<DataArray type=\"Float32\" Name=\"drainage rate(mm/h)\" format=\"ascii\">" << endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetDrainageRate() << " ";
            if(cell->GetValid()) {
                double value = cell->GetDrainageRate();
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        }
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;      

    //soil water storage potential(mm)
    out_file << "<DataArray type=\"Float32\" Name=\"soil water storage potential(mm)\" format=\"ascii\">" << endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetStoragePotential() << " ";
            if(cell->GetValid()) {
                double value = cell->GetStoragePotential();
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        }
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl; 
    //total accumulation of sediment volume
    out_file<<"<DataArray type=\"Float32\" Name=\"totalaccum_sediment\" format=\" ascii\">"<<endl;//(m3)volume
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->Getaccum_totsed() << " ";
            if(cell->GetValid()) {
                double value = cell->Getaccum_totsed();
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";

        }
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl; 
   
    //i added
    //vector<double> SusVol = cell->GetSuspendedSedimentVolumes(); //^
    //vector<double> BmVol = cell->GetBedMaterialVolumes(); //^
    //vector<double> EroVol = cell->GetErodedMaterialVolumes(); //^
    //suspended sediment volumes
    //clay
    out_file<<"<DataArray type=\"Float32\" Name=\"suspended_volume_clay\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetSuspendedSedimentVolumes()[0] << " ";
            if(cell->GetValid()) {
                double value = cell->GetSuspendedSedimentVolumes()[0];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    } 
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;     
    //silt
    out_file<<"<DataArray type=\"Float32\" Name=\"suspended_volume_silt\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetSuspendedSedimentVolumes()[1] << " ";
            if(cell->GetValid()) {
                double value = cell->GetSuspendedSedimentVolumes()[1];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl; 
    //sand
    out_file<<"<DataArray type=\"Float32\" Name=\"suspended_volume_sand\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetSuspendedSedimentVolumes()[2] << " ";
            if(cell->GetValid()) {
                double value = cell->GetSuspendedSedimentVolumes()[2];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    } 
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl; 
    //gravel
    out_file<<"<DataArray type=\"Float32\" Name=\"suspended_volume_gravel\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetSuspendedSedimentVolumes()[3] << " ";
            if(cell->GetValid()) {
                double value = cell->GetSuspendedSedimentVolumes()[3];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;     

    //bed material volumes
    //clay
    out_file<<"<DataArray type=\"Float32\" Name=\"bed_volume_clay\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetBedMaterialVolumes()[0] << " ";
            if(cell->GetValid()) {
                double value = cell->GetBedMaterialVolumes()[0];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    } 
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;     
    //silt
    out_file<<"<DataArray type=\"Float32\" Name=\"bed_volume_silt\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetBedMaterialVolumes()[1] << " ";
            if(cell->GetValid()) {
                double value = cell->GetBedMaterialVolumes()[1];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl; 
    //sand
    out_file<<"<DataArray type=\"Float32\" Name=\"bed_volume_sand\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetBedMaterialVolumes()[2] << " ";
            if(cell->GetValid()) {
                double value = cell->GetBedMaterialVolumes()[2];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    } 
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl; 
    //gravel
    out_file<<"<DataArray type=\"Float32\" Name=\"bed_volume_gravel\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetBedMaterialVolumes()[3] << " ";
            if(cell->GetValid()) {
                double value = cell->GetBedMaterialVolumes()[3];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;   

    //eroded material volumes
    //clay
    out_file<<"<DataArray type=\"Float32\" Name=\"eroded_volume_clay\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetErodedMaterialVolumes()[0] << " ";
            if(cell->GetValid()) {
                double value = cell->GetErodedMaterialVolumes()[0];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    } 
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;     
    //silt
    out_file<<"<DataArray type=\"Float32\" Name=\"eroded_volume_silt\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetErodedMaterialVolumes()[1] << " ";
            if(cell->GetValid()) {
                double value = cell->GetErodedMaterialVolumes()[1];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl; 
    //sand
    out_file<<"<DataArray type=\"Float32\" Name=\"eroded_volume_sand\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetErodedMaterialVolumes()[2] << " ";
            if(cell->GetValid()) {
                double value = cell->GetErodedMaterialVolumes()[2];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    } 
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl; 
    //gravel
    out_file<<"<DataArray type=\"Float32\" Name=\"eroded_volume_gravel\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetErodedMaterialVolumes()[3] << " ";
            if(cell->GetValid()) {
                double value = cell->GetErodedMaterialVolumes()[3];
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        } 
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;        

    //change of land elevation
       out_file<<"<DataArray type=\"Float32\" Name=\"elevation_change(y)\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetLandElavationChange() << " ";
            if(cell->GetValid()) {
                double value = cell->GetLandElavationChange();
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        }
     out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl; 

    //cumulative_elevaion_change
    out_file<<"<DataArray type=\"Float32\" Name=\"cumulative_elevation_change(m)\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetCumulativLandElavationChange() << " ";
            if(cell->GetValid()) {
                double value = cell->GetCumulativLandElavationChange();
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        }
        out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl;  

    ////Save land elevation update to vti
      out_file<<"<DataArray type=\"Float32\" Name=\"land_elevation(m)\" format=\"ascii\">"<<endl;
    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            //if(cell->GetValid()) out_file <<std::setprecision(10)<< cell->GetLandElevation() << " ";
            if(cell->GetValid()) {
                double value = cell->GetLandElevation();
                if(fabs(value) < 1.0e-12 && value != 0.) value = value/fabs(value) * 1.0e-12;
                out_file <<std::setprecision(10)<< value << " ";
            }               
            else out_file << "-9999" << " ";
        }
     out_file << endl;
    }
    out_file << endl;
    out_file << "</DataArray>" << endl <<endl; 
    //
    out_file << "</CellData>" << endl;
    out_file << "</Piece>" << endl;
    out_file << "</ImageData>" << endl;
    out_file << "</VTKFile>" << endl;

    // close file
    out_file.close();    
}


void RunCA::SaveWaterDepthToText()
{
    string file_name(model_name_);
    int minute = round(model_time_/60.);
    string time = to_string(minute);
    file_name += ("-output-water_depth-" + time + "min.txt");
    ofstream fflowdepth(file_name.c_str());
    fflowdepth << "ncols " << nx_ << endl << "nrows " << ny_ << endl;
    fflowdepth << "xllcorner 0" << endl << "yllcorner 0" << endl;
    fflowdepth << "cellsize " << cell_size_ << endl;
    fflowdepth << "NODATA_value -9999" << endl;
    for (int y=1;y<(ny_+1);y++) {
        for (int x=1;x<(nx_+1);x++) {
            if (!cells_[y][x]->GetValid()) {
                fflowdepth << -9999 << " ";
            }
            else {
                 fflowdepth << setprecision(10) << cells_[y][x]->GetWaterDepth()<< " ";
            }
        }
        fflowdepth << endl;
    }
      
    fflowdepth.close();   
}

// save cumulative elevation change to text to test effect of land elevation change on water depth
void RunCA::SaveCumulativeLandElavationChangeToText()
{
    string file_name(model_name_);
    int minute = round(model_time_/60.);
    string time = to_string(minute);
    file_name += ("-output-cumulativeelavtion_change-" + time + "min.txt");
    ofstream fcelevation(file_name.c_str());
    fcelevation << "ncols " << nx_ << endl << "nrows " << ny_ << endl;
    fcelevation << "xllcorner 0" << endl << "yllcorner 0" << endl;
    fcelevation << "cellsize " << cell_size_ << endl;
    fcelevation << "NODATA_value -9999" << endl;
    for (int y=1;y<(ny_+1);y++) {
        for (int x=1;x<(nx_+1);x++) {
            if (!cells_[y][x]->GetValid()) {
                fcelevation << -9999 << " ";
            }
            else {
                 fcelevation << setprecision(10) << cells_[y][x]->GetCumulativLandElavationChange()<< " ";
            }
        }
        fcelevation << endl;
    }
      
   fcelevation.close();   
}

//*****
void RunCA::SaveElevationUpdateToText(){
     string file_name(model_name_);
    int minute = round(model_time_/60.);
    string time = to_string(minute);
    file_name += ("output-Elevation_update" + time + "min.txt");
    ofstream felevationupdate(file_name.c_str());
    felevationupdate << "ncols " << nx_ << endl << "nrows " << ny_ << endl;
    felevationupdate << "xllcorner 0" << endl << "yllcorner 0" << endl;
    felevationupdate << "cellsize " << cell_size_ << endl;
    felevationupdate << "NODATA_value -9999" << endl;
    for (int y=1;y<(ny_+1);y++) {
        for (int x=1;x<(nx_+1);x++) {
            if (!cells_[y][x]->GetValid()) {
                felevationupdate << -9999 << " ";
            }
            else {
                 felevationupdate << setprecision(10) << cells_[y][x]-> GetLandElevation()<< " ";
            }
        }
        felevationupdate << endl;
    }
      
    felevationupdate.close(); 

}

void RunCA::SaveTotalAccumulationSedimentVolumeToText(){
     string file_name(model_name_);
    int minute = round(model_time_/60.);
    string time = to_string(minute);
    file_name += ("-output-Total_Accumulation_Sediment_volume-" + time + "min.txt");
    ofstream ftotaccumulationsediment(file_name.c_str());
    ftotaccumulationsediment << "ncols " << nx_ << endl << "nrows " << ny_ << endl;
    ftotaccumulationsediment << "xllcorner 0" << endl << "yllcorner 0" << endl;
    ftotaccumulationsediment << "cellsize " << cell_size_ << endl;
    ftotaccumulationsediment << "NODATA_value -9999" << endl;
    for (int y=1;y<(ny_+1);y++) {
        for (int x=1;x<(nx_+1);x++) {
            if (!cells_[y][x]->GetValid()) {
                ftotaccumulationsediment << -9999 << " ";
            }
            else {
                 ftotaccumulationsediment << setprecision(10) << cells_[y][x]->Getaccum_totsed()<< " ";
            }
        }
        ftotaccumulationsediment << endl;
    }
      
    ftotaccumulationsediment.close(); 


}

//*****


void RunCA::WriteOutputFileHeaders(){
    //create output files and write headers
    string file_name; 
    //infiltration rate
    file_name = model_name_;
    file_name += "-output-infiltration_rate.txt";
    ofstream ifs(file_name.c_str());
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    ifs << "Time(min)" << "\t";
    for(int n=0;n<iparams_.size();n++) ifs<<"Soiltype"<<n+1<<"\t";
    ifs<<"Mean"<<endl;
    ifs.close();
    ifs.clear();
    //infiltration capacity
    file_name = model_name_;
    file_name += "-output-infiltration_capacity.txt";
    ifs.open(file_name.c_str());
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    ifs << "Time(min)" << "\t";
    for(int n=0;n<iparams_.size();n++) ifs<<"Soiltype"<<n+1<<"\t";
    ifs<<"Mean"<<endl;
    ifs.close();
    ifs.clear();
    //drainage rate
    file_name = model_name_;
    file_name += "-output-drainage_rate.txt";
    ifs.open(file_name.c_str());
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    ifs << "Time(min)" << "\t";
    for(int n=0;n<iparams_.size();n++) ifs<<"Soiltype"<<n+1<<"\t";
    ifs<<"Mean"<<endl;
    ifs.close();
    ifs.clear();
    //runoff
    file_name = model_name_;
    file_name += "-output-runoff.txt";
    ifs.open(file_name.c_str());
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    /*ifs<<"Time(min)\tRainfall(mm/h)\tRunoff(cm2/min)\tCumulative_rain(mm)\tCumulative_runoff(mm)\t"//*/
    ifs<<"Time(min)\tRainfall(mm/h)\tRunoff(cm2/min)\trunoff_litre\trunoff_mm_hour_cellc\trunoff_timestep_mm_min\ttrunoff_timestep_L_s\trunoff_timestep_mm_hour\tCumulative_rain(mm)\tCumulative_runoff(cm2)\t"// added
       <<"Cumulative_infiltration(mm)\tCumulative_drainage(mm)\tPonded_water(mm)\tAccumulated_runoff(L)\trunoff_ratem3"<<endl;
   
    ifs<<"0\t0\t0\t0\t0\t0\t0\t0"<<endl; //Qi added - 6/11   
    ifs.close();
    ifs.clear(); 
//output-transportedsediment_rate
    file_name = model_name_;
    file_name += "-output-transportedsediment_rate.txt";
    ifs.open(file_name.c_str());
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    //ifs<<"Time(min)\tsum_transportedsediment()\ttransportedsediment_rate()"<<endl;
    //ifs<<"Time(min)\tsediment_yield_rate(g/min/cm2)"<<endl; //Qi: we only need sediment_yield_rate
    ifs<<"Time(min)\tsediment_yield_rate(g/min/m)\ttransportedsediment_rateg_m_m2\ttransportedsediment_rateg_min\t transportedsediment_kg\ttransportedsediment_timestep_gmin\taccumulated_sediment(g/m2)\tsediment_concetration(g/L)"<<endl; //added 
    //ifs<<"0\t0"<<endl; //Qi added - 6/11
    ifs<<"0\t0\t0"<<endl; //added
    ifs.close();
    ifs.clear(); 
    //slope
    file_name = model_name_;
    file_name += "-output-slope.txt";
    ifs.open(file_name.c_str());
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    //ifs<<"Time(min)\tsum_transportedsediment()\ttransportedsediment_rate()"<<endl;
    ifs<<"Time(min)\tslope"<<endl; //Qi: we only need sediment_yield_rate
    ifs.close();
    ifs.clear();  

       //#####################  //output-bedloadtransportedsediment_rate
        file_name = model_name_;
    file_name += "-output-bedloadtransportedsediment_rate.txt";
    ifs.open(file_name.c_str());
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
   
    ifs<<"Time(min)\tsediment_yield_rate(g/min/m)b\ttransportedsediment_rateg_m_m2b\ttransportedsediment_rateg_minb\t transportedsediment_kgb"<<endl; //added 
    //ifs<<"0\t0"<<endl; //Qi added - 6/11
    ifs<<"0\t0\t0"<<endl; //added
    ifs.close();
    ifs.clear(); 
   
    //#####################    
    
    
}


void RunCA::WriteResults()
{
    //infiltration and drainage recorders
    int soiltypes=iparams_.size();
    vector<double> Ie_rec(soiltypes,0.), Ic_rec(soiltypes,0.), Dr_rec(soiltypes,0.);
    vector<int> Cells_rec(soiltypes,0);
    tot_ponding_ = 0.;

    for(int y=1;y<ny_+1;y++){
        for(int x=1;x<nx_+1;x++){
            Cell* cell = cells_[y][x];
            if(cell->GetValid()){
                int soiltype = cell->GetSoilType();
                Cells_rec.at(soiltype-1) = Cells_rec[soiltype-1] + 1;
                Ie_rec.at(soiltype-1) = Ie_rec[soiltype-1] + cell->GetInfiltrationRate();
                Ic_rec.at(soiltype-1) = Ic_rec[soiltype-1] + cell->GetInfiltrationCapacity();
                Dr_rec.at(soiltype-1) = Dr_rec[soiltype-1] + cell->GetDrainageRate();
                tot_ponding_ += cell->GetWaterDepth();
            }
        }
    }
    
    //write infiltration and drainage rates to files
    string file_name; 
    //infiltration rate
    file_name = model_name_;
    file_name += "-output-infiltration_rate.txt";
    ofstream ifs;
    ifs.open(file_name.c_str(), ofstream::app);
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    ifs << model_time_/60. << "\t";
    double sum_ie(0.);
    for(int n=0;n<soiltypes;n++) {ifs<<Ie_rec[n]/Cells_rec[n]<<"\t"; sum_ie+=Ie_rec[n]/Cells_rec[n];}
    ifs<<sum_ie/soiltypes<<endl;
    ifs.close();
    ifs.clear();   
    //infiltration capacity
    file_name = model_name_;
    file_name += "-output-infiltration_capacity.txt";
    ifs.open(file_name.c_str(), ofstream::app);
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    ifs << model_time_/60. << "\t";
    double sum_ic(0.);
    for(int n=0;n<soiltypes;n++) {ifs<<Ic_rec[n]/Cells_rec[n]<<"\t"; sum_ic+=Ic_rec[n]/Cells_rec[n];}
    ifs<<sum_ic/soiltypes<<endl;
    ifs.close();
    ifs.clear();   
    //drainage rate
    file_name = model_name_;
    file_name += "-output-drainage_rate.txt";
    ifs.open(file_name.c_str(), ofstream::app);
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    ifs << model_time_/60. << "\t";
    double sum_dr(0.);
    for(int n=0;n<soiltypes;n++) {ifs<<Dr_rec[n]/Cells_rec[n]<<"\t"; sum_dr+=Dr_rec[n]/Cells_rec[n];}
    ifs<<sum_dr/soiltypes<<endl;
    ifs.close();
    ifs.clear();       
    //cout<<"  passed write infiltration"<<endl;
    //write rain and runoff
    int minute = round(model_time_/60.);
    //double rain_rate = rainfall_[minute-1]; //mm/h   
    double rain_rate = rainfall_[minute]; //mm/h  
    double sum_runoff(0.);
    double sum_runoffL(0.);
     double sum_runoffm3(0.);
    double sum_runoffm(0.);
    double sum_runoffmm(0.);
    //double runoff_mm(0.);
    int cell_count(0); //Qi added
     double runoff_litre(0);
    for (int y=ny_,x=1,runoff=0;x<(nx_+1);x++) {
    //for (int y=ny_-100,x=1,runoff=0;x<(nx_+1);x++) { //added
      
        Cell* cell=cells_[y][x];
        if(cell->GetValid()){
            double water_depth = cell->GetWaterDepth(); //mm
            vector<double> flows = cell->GetFlows();
            double outflow = flows[6]/(cell_size_*1.414) + flows[7]/cell_size_ + flows[8]/(cell_size_*1.414); //m2
            sum_runoff += outflow * 100. * 100.; //(m2 -> cm2)
            //sum_runoff = sum_runoff/10. * cell_size_*100.; //cm2 //added 27/10
            //sum_runoff = sum_runoff/10.; //cm //added


            //shahla added for getting runoff as a liter and liter/s
           double outflowm3 = flows[6] + flows[7]+ flows[8]; //m3-
           double outflowL=(flows[6]+flows[7]+flows[8])*1000. ; //LITER
         
           double outflow_mm=(flows[6]+flows[7]+flows[8]) / cell_size_ / cell_size_/effective_cells_ * 1000; //runoff amount mm
          
           sum_runoffL+=outflowL;//L
          sum_runoffm3+=outflowm3;//m3
          sum_runoffmm+=outflow_mm; //mm
          
           
          
          

            
        }
        cell_count++; //Qi added
    }
     
    
    //double runoff_rate = sum_runoff/(time_increment_/60.); //cm2/min //added 27/10
    double runoff_rate = sum_runoff/cell_count/(time_increment_/60.); //cm2/min //added
    
   //test21 double runoff_ratem3=sum_runoff/100./100.*cell_size_*cell_size_*cell_count;// cm2/min -> m3
    double runoff_ratem3=sum_runoffm3; //m3
   
    //shahla added
    //double runoff_mm=sum_runoffm3/ cell_size_ / cell_size_/effective_cells_ * 1000.;//runoffrate mm
     double runoff_mm=sum_runoffm3/ cell_size_ / cell_size_/cell_count * 1000.;//runoffrate mm
    runoff_litre=sum_runoffL ; //litre
    
     
     double runoff_mm_hour_cellc=sum_runoffm3/ cell_size_ / cell_size_/effective_cells_ * 1000/cell_count/(time_increment_/3600.); //mm/h
     double runoff_timestep_mm_min=sum_runoffm3/ (cell_size_ * cell_size_*effective_cells_) * 1000. / (time_increment_/60.); //mm/min
     //double runoff_timestep_mm_min=sum_runoffm3/ (cell_size_ * cell_size_*effective_cells_) * 1000./cell_count / (time_increment_/60.); //mm/min
     
     double runoff_timestep_L_s=sum_runoffm3 * 1000. / (time_increment_); //L/s
     //double runoff_timestep_L_s=sum_runoffm3 * 1000. /cell_count/ (time_increment_); //L/s
     double runoff_timestep_mm_hour=sum_runoffm3/ (cell_size_ * cell_size_*effective_cells_) * 1000. / (time_increment_/3600.);; //mm/hour
 
    file_name = model_name_;
    file_name += "-output-runoff.txt";
    ifs.open(file_name.c_str(), ofstream::app);
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    ifs<<model_time_/60.<<"\t"<<rain_rate<<"\t"<<runoff_rate<<"\t"<<runoff_litre<<"\t"<<runoff_mm_hour_cellc<<"\t"<<runoff_timestep_mm_min<<"\t"<< runoff_timestep_L_s<<"\t"<<runoff_timestep_mm_hour<<"\t"<<tot_rain_/effective_cells_<<"\t"//
      
       <<tot_runoff_/10./cell_count*cell_size_*100.<<"\t"<<tot_infiltration_/effective_cells_<<"\t"// added
       <<tot_drainage_/effective_cells_<<"\t"<<tot_ponding_/effective_cells_<<"\t"<<tot_runoff_<<"\t"<<runoff_ratem3<<endl;
    ifs.close();
    ifs.clear(); 
    //cout<<"  passed write runoff"<<endl;
     double sum_concentration(0.);  
     // write sedimentoutput
    double sum_transportedsediment(0.);
    double totalvolume(0.);
    //int cell_count(0); //Qi added
    for (int y=ny_,x=1,transportedsediment=0;x<(nx_+1);x++) {
    //for (int y=ny_-100,x=1,transportedsediment=0;x<(nx_+1);x++) {
        //for (int y=ny_-3,x=1,transportedsediment=0;x<(nx_+1);x++) {
        Cell* cell=cells_[y][x];
        if(cell->GetValid()){
            vector<double> Stot = cell->GettotalSedimentTransport() ;
            ////*****conagust
            vector<double>sus_sedi_outflow = cell->GetSusSediOutflow();
            double concen = sus_sedi_outflow[6]+sus_sedi_outflow[7]+sus_sedi_outflow[8];//m3
            sum_concentration+=concen;
            ////**********
            
            //sum_transportedsediment += (Stot[6]+Stot[7]+Stot[8]); //sedimentdischargeamount amount (m3)
            double sediment = (Stot[6]+Stot[7]+Stot[8]); //m3
            totalvolume+=(Stot[6]+Stot[7]+Stot[8]);////test
            //m3 -> g/min/m2
	        int soiltype = cell->GetSoilType();
	        double rho = erosion_para_[soiltype - 1][7];//solid density g/cm3
        
            rho *= 1.0e6; //g/cm3 -> g/m3           
            //sediment /= (time_increment_/60.); //m3/min
            sediment *= rho; //m3 -> g
            //sediment /= (cell_size_*cell_size_); //g/min -> g/min/m2
            sum_transportedsediment += sediment; //g
            cell_count++; //Qi added
            
           ///*****************
             
            
           ///*****************



        }
    }
    
    //double transportedsediment_rate = sum_transportedsediment/time_increment_*3600./effective_cells_; //m3/h??
    //double transportedsediment_rate = sum_transportedsediment / (time_increment_/60.) / (effective_cells_ * cell_size_*cell_size_); //g -> g/min/m2
    //Qi added: use sediment discharge cells instead of all cells
    double transportedsediment_rate (0.);
    double transportedsediment_rateg_m_m2 (0.);
    double sediment_concentration (0.);
    double transportedsediment_kg (0.);
    double transportedsediment_timestep_gmin (0.);
   double transportedsediment_rateg_min (0.);
    if(cell_count != 0) {
        double tsm3=sum_transportedsediment / (time_increment_/60.) ;//gr ->gr/min
        double tsm2 = sum_transportedsediment / (time_increment_/60.) / (cell_count* cell_size_*cell_size_); //g -> g/min/m2
         double ts = sum_transportedsediment / (time_increment_/60.) / (cell_count* cell_size_); //g -> g/min/m
        //double ts = sum_transportedsediment/cell_count*100.; //g/min/cm -> g/min/m //added 
        //transportedsediment_rate =  ts/10000.; //g/min/m2 -> g/min/cm2shahla adeed 9/11/2020
        //transportedsediment_rate =  ts/1000.; //g/min/m2 -> g/min/cm2
        transportedsediment_rate =  ts; //added
       transportedsediment_kg=sum_transportedsediment /1000.; //sedimen gr->kg
       transportedsediment_timestep_gmin=sum_transportedsediment/ (time_increment_/60.) ; //sediment gr/min
       transportedsediment_rateg_m_m2 =tsm2;//g/min/m2
        transportedsediment_rateg_min =tsm3;//g/min
       // if(transportedsediment_rateg_min_m2> 0.)cout<<"i = "<<transportedsediment_rateg_min_m2<<endl;
                
        sediment_concentration = sum_transportedsediment/(sum_runoffm3 * 1000.); //g/L 
          //sediment_concentration =totalvolume/ runoff_ratem3;
          //sediment_concentration=sum_concentration;
      
    }

    //double accumulated_sediment = tot_transportedsediment_; //g //added
    //accumulated_sediment = accumulated_sediment/cell_count/(cell_size_*cell_size_); //g/m2 //added  

    double accumulated_sediment = tot_transportedsediment_/cell_count*100.; //g/cm2 -> g/m2 //added  

     

    file_name = model_name_;
    file_name += "-output-transportedsediment_rate.txt";
    ifs.open(file_name.c_str(), ofstream::app);
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    //
    //ifs<<model_time_/60.<<"\t"<<transportedsediment_rate<<"\t"<<tot_transportedsediment_/effective_cells_<<endl;
    //ifs<<model_time_/60.<<"\t"<<transportedsediment_rate<<endl; //Qi: we only need sediment_yield_rate
    ifs<<model_time_/60.<<"\t"<<transportedsediment_rate<<"\t"<<transportedsediment_rateg_m_m2<<"\t"<<transportedsediment_rateg_min<<"\t"<<transportedsediment_kg<<"\t"<<transportedsediment_timestep_gmin<<"\t"<<accumulated_sediment<<"\t"<<sediment_concentration<<endl;
    //
    
    ifs.close();
    ifs.clear(); 
    //cout<<"  passed write erosion"<<endl;
    /////*******************************************************************************

    
    //####################################### write Bedload result
// write  bed load sedimentoutput
    double sum_transportedsedimentb(0.);
    double totalvolumeb(0.);
    //int cell_count(0); //Qi added
    for (int y=ny_,x=1,transportedsedimentb=0;x<(nx_+1);x++) {
   
        Cell* cell=cells_[y][x];
        if(cell->GetValid()){
            vector<double> Stotb = cell->GettotalSedimentTransportb() ;
            
            vector<double>sus_sedi_outflow = cell->GetSusSediOutflow();
            double concen = sus_sedi_outflow[6]+sus_sedi_outflow[7]+sus_sedi_outflow[8];//m3
            sum_concentration+=concen;
            ////**********
            
            //sum_transportedsediment += (Stot[6]+Stot[7]+Stot[8]); //sedimentdischargeamount amount (m3)
            double sedimentb = (Stotb[6]+Stotb[7]+Stotb[8]); //m3
            totalvolumeb+=(Stotb[6]+Stotb[7]+Stotb[8]);////test
            //m3 -> g/min/m2
	        int soiltype = cell->GetSoilType();
	        double rhob = erosion_para_[soiltype - 1][7];//solid density g/cm3
        
            rhob *= 1.0e6; //g/cm3 -> g/m3           
            //sediment /= (time_increment_/60.); //m3/min
            sedimentb *= rhob; //m3 -> g
            //sediment /= (cell_size_*cell_size_); //g/min -> g/min/m2
            sum_transportedsedimentb += sedimentb; //g
            cell_count++; //Qi added
            
         


        }
    }
    
 
    double transportedsediment_rateb (0.);
    double transportedsediment_rateg_m_m2b (0.);
   
    double transportedsediment_kgb (0.);
    double transportedsediment_timestep_gminb (0.);
   double transportedsediment_rateg_minb (0.);
    if(cell_count != 0) {
        double tsm3b=sum_transportedsedimentb / (time_increment_/60.) ;//gr ->gr/min
        double tsm2b = sum_transportedsedimentb / (time_increment_/60.) / (cell_count* cell_size_*cell_size_); //g -> g/min/m2
         double tsb = sum_transportedsedimentb / (time_increment_/60.) / (cell_count* cell_size_); //g -> g/min/m
      
        transportedsediment_rateb =  tsb; //added
       transportedsediment_kgb=sum_transportedsedimentb /1000.; //sedimen gr->kg
       transportedsediment_timestep_gminb=sum_transportedsedimentb/ (time_increment_/60.) ; //sediment gr/min
       transportedsediment_rateg_m_m2b =tsm2b;//g/min/m2
        transportedsediment_rateg_minb =tsm3b;//g/min
      
    }
        file_name = model_name_;
    file_name += "-output-Bedloadtransportedsediment_rate.txt";
    ifs.open(file_name.c_str(), ofstream::app);
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    //
   
    ifs<<model_time_/60.<<"\t"<<transportedsediment_rateb<<"\t"<<transportedsediment_rateg_m_m2b<<"\t"<<transportedsediment_rateg_minb<<"\t"<<transportedsediment_kgb<<endl;
    //
    
    ifs.close();
    ifs.clear();
  //cout<<"  passed write bedload"<<endl;
    ////
    /*
    
    ///write elevation change
    //slope out put
    double sum_slope(0.);
    cell_count=0; //Qi added
    for (int y=ny_,x=1,transportedsediment=0;x<(nx_+1);x++) {
        Cell* cell=cells_[y][x];
        if(cell->GetValid()){
            vector<double> Slope = cell->GetSlope() ;
            sum_slope += Slope[7]; //downslope
            cell_count++; //Qi added
        }
    }

    if(cell_count != 0) {
        sum_slope /= cell_count;
    }

    file_name = model_name_;
    file_name += "-output-slope.txt";
    ifs.open(file_name.c_str(), ofstream::app);
    if (!ifs.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}
    //
    //ifs<<model_time_/60.<<"\t"<<transportedsediment_rate<<"\t"<<tot_transportedsediment_/effective_cells_<<endl;
    ifs<<model_time_/60.<<"\t"<<sum_slope<<endl; //Qi: we only need sediment_yield_rate
    //
    
    ifs.close();
    ifs.clear(); 
    cout<<"  passed write slope"<<endl;
    */
    /////  
}


void RunCA::WriteSummary()
{
    string file_name; 
    file_name = model_name_;
    file_name += "-output-summary.txt";
    ofstream fsummary;
    fsummary.open(file_name.c_str());
    if (!fsummary.is_open()) {cerr<<"ERROR: unable to open "<<file_name<<endl; terminate();}

    fsummary << "\n******** Configurations ********" << endl;
    fsummary <<"Plot_size = "<<cell_size_*nx_<<" m(X) x "<<cell_size_*ny_<<" m(Y) = "<<cell_size_*nx_*cell_size_*ny_<<" m2"<<endl;
    fsummary <<"Cell_side_length = "<<cell_size_<<" m"<<endl;
    fsummary <<"Cell_numbers = "<<nx_<<" (X) x "<<ny_<<" (Y) = "<<nx_*ny_<<endl;
    fsummary <<"Effective_cell_numbers = "<<effective_cells_<<endl;
    fsummary <<"Soil_types = "<<iparams_.size()<<endl;
    fsummary <<"Simulation_duration = "<<duration_<<" min"<<endl<<endl;

    fsummary << "\n******** Results ********" << endl;
    fsummary <<"Iteration_steps = "<<iteration_<<endl;
    fsummary <<"Mean_time_step = "<<duration_*60./iteration_<<" sec"<<endl; 
    fsummary <<"Adjusting_coefficient = "<<adjust_coefficient_<<endl;   
    fsummary <<"Execution_time = "<<run_time_/60.<<" min"<<endl; 
    fsummary <<"Number of threads in use = "<<threads_<<endl<<endl;
    fsummary.setf(ios::showpoint);
    fsummary.setf(ios::fixed);
    fsummary << "Total_rainfall_amount = "<<setprecision(3)<<tot_rain_/effective_cells_<<" mm"<<endl;
    fsummary << "Total_runoff_amount = "<<setprecision(3)<<tot_runoff_/effective_cells_<<" mm"<<endl;
    fsummary << "Total_infiltration_amount = "<<setprecision(3)<<tot_infiltration_/effective_cells_<<" mm"<<endl;
    fsummary << "Total_drainage_amount = "<<setprecision(3)<<tot_drainage_/effective_cells_<<" mm"<<endl;
    fsummary << "Total_ponding_amount = "<<setprecision(3)<<tot_ponding_/effective_cells_<<" mm"<<endl;
    fsummary << "Runoff_coefficient = "<<setprecision(3)<<tot_runoff_/tot_rain_*100.<<"%"<<endl;
    fsummary << "Ponding_coefficient = "<<setprecision(3)<<tot_ponding_/tot_rain_*100.<< "%"<<endl;
    fsummary << "Infiltration_coefficient = "<<setprecision(3)<<tot_infiltration_/tot_rain_*100.<< "%"<<endl;
    fsummary << "Drainage_coefficient = "<<setprecision(3)<<tot_drainage_/tot_rain_*100.<< "%"<<endl;
    //*****summary of sediment output (total)
    fsummary << "Total_sediment_amount = "<<setprecision(3)<<tot_transportedsediment_/effective_cells_<<" "<<endl;
    //how much sediment left the plot
 
    fsummary.close();
}




