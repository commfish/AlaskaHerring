// -------------------------------------------------------------------------- //
//               Age-structured model for Alaska herring stocks               //
//                                                                            //
//                               VERSION 0.1                                  //
//                                Jan  2015                                   //
//                                                                            //
//                                 AUTHORS                                    //
//                              Sherri Dressel                                //                                 
//                        sherri.dressel@alaska.gov                           //
//                               Sara Miller                                  //
//                          sara.miller@alaska.gov                            //
//                               Kray Van Kirk                                //
//                          kray.vankirk@alaska.gov                           //
//                                                                            //
//                   Built on code developed by Peter Hulson                  //
//                            pete.hulson@noaa.gov                            //
//                                                                            //
//                           Layout and references                            //
//                              Steven Martell                                //
//                          martell.steve@gmail.com                           //
//                                                                            //
// CONVENTIONS: Formatting conventions are based on                           //
//              The Elements of C++ Style (Misfeldt et al. 2004)              //
//                                                                            //
//                                                                            //
// NAMING CONVENTIONS:                                                        //
//             Macros       -> UPPERCASE                                      //
//             Constants    -> UpperCamelCase                                 //
//             Functions    -> lowerCamelCase                                 //
//             Variables    -> lowercase                                      //
//                                                                            //
//                                                                            //
//                                                                            //
// -------------------------------------------------------------------------- //
//-- CHANGE LOG:                                                             -//
//--  Jan 2015 - revision of legacy code::                                   -//
//--               :variable naming conventions                              -//
//--               :intra-annual calendar                                    -//
//--               :standardization of units across stocks                   -//
//--               :modification for potential code distribution             -//
//--  May 2015 - added code and scripts to a git-hub repo.                   -//
//--                                                                         -//
// -------------------------------------------------------------------------- //
// small tweak

DATA_SECTION

// |---------------------------------------------------------------------------|
// | CHECK FOR OPTIONAL COMMAND LINE ARGUMENTS & SET FLAGS
// |---------------------------------------------------------------------------|
// | b_simulation_flag  -> flag for running in simulation mode
// | rseed      -> random number seed for simulation 
  int b_simulation_flag;
  int rseed;
  int retro_yrs;
  LOCAL_CALCS
    int on = 0;
    retro_yrs = 0;
    rseed  = 0;

    b_simulation_flag = 0;
    if (ad_comm::argc > 1)
    {
      int on = 0;
      if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-sim")) > -1 )
      {
        b_simulation_flag = 1;
        rseed = atoi(ad_comm::argv[on+1]);
      }

      if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-retro")) > -1 )
      {
        retro_yrs = atoi(ad_comm::argv[on+1]);
        cout<<"|—————————————————————————————————————————————————|\n";
        cout<<"| Implementing Retrospective analysis             |\n";
        cout<<"|—————————————————————————————————————————————————|\n";
        cout<<"| Number of retrospective years = "<<retro_yrs<<endl;
      }

    }

  END_CALCS
// |---------------------------------------------------------------------------|


// |---------------------------------------------------------------------------|
// | STRINGS FOR INPUT FILES                                                   |
// |---------------------------------------------------------------------------|
// |- The files below are listed in the 'model.dat' file;
// |   nothing else should be in the 'model.dat' file
// |- These files should be named by the stock they are modeling;
// |   example: "sitka.dat", "seymour.dat"
// |- DO NOT use two word names such as "seymour cove.dat"
// |
// | DEBUG_FLAG             : Boolean Flag used for manual debugging
// | DataFile               : data to condition the assessment model    
// | ControlFile            : controls for years, phases, and block options 
  init_int DEBUG_FLAG;
  init_adstring DataFile;      
  init_adstring ControlFile;


// |---------------------------------------------------------------------------|
// | READ CONTENTS OF DATA FILE
// |---------------------------------------------------------------------------|
  !! ad_comm::change_datafile_name(DataFile);

// |---------------------------------------------------------------------------|
// | Model dimensions.
// |---------------------------------------------------------------------------|
// | dat_syr      -> first year of data
// | day_nyr      -> last year of data
// | mod_syr      -> first year of model
// | mod_nyr      -> last year of model
// | sage         -> first age class
// | nage         -> plus group age class
  init_int dat_syr;
  init_int dat_nyr;
  init_int mod_syr;
  init_int mod_nyr;
  init_int sage;
  init_int nage;
  int rec_syr;
  !!  rec_syr = mod_syr + sage;

  vector age(sage,nage);
  !! age.fill_seqadd(sage,1);

// |---------------------------------------------------------------------------|
// | Fecundity regression coefficients
// |---------------------------------------------------------------------------|
// | - 
  init_int nFecBlocks;
  init_ivector nFecBlockYears(1,nFecBlocks);
  init_vector fec_slope(1,nFecBlocks);
  init_vector fec_inter(1,nFecBlocks);
  // !!COUT(fec_slope)

// |---------------------------------------------------------------------------|
// | Time series data.  Catch in short tons. Comps in proportions.
// |---------------------------------------------------------------------------|
// | data_catch   -> Catch: colnames(Year,Catch,log.se) units are short tons
// | data_sp_waa  -> Spawn Weight-at-age (Year,weight-at-age)
// | data_cm_waa  -> Commercial catch weight-at-age (Year, weight-at-age)
// | data_cm_comp -> Commercial catch composition (Year, age proporitons)
// | data_sp_comp -> Spawn Sample age composition (Year, age proportions)
// | data_egg_dep -> Egg deposition survey (Year, Index, log.se)
// | avg_sp_waa   -> Average weight-at-age for spawning biomass.
  init_matrix  data_ct_raw(dat_syr,dat_nyr,1,3);
  init_matrix  data_sp_waa(dat_syr,dat_nyr,sage-1,nage);
  init_matrix  data_cm_waa(dat_syr,dat_nyr,sage-1,nage);
  init_matrix data_cm_comp(dat_syr,dat_nyr,sage-1,nage);
  init_matrix data_sp_comp(dat_syr,dat_nyr,sage-1,nage);
  init_matrix data_egg_dep(dat_syr,dat_nyr,1,3);
  init_matrix data_mileday(dat_syr,dat_nyr,1,3);
  
  // Calculate average spawner weight-at-age for SR parameters
  vector avg_sp_waa(sage,nage);
  LOCAL_CALCS
    int n = data_sp_waa.rowmax() - data_sp_waa.rowmin() + 1;
    avg_sp_waa = colsum(data_sp_waa)(sage,nage) / n;
  END_CALCS

  // Calculate Fecundity-at-age based on regression coefficients.
  matrix Eij(mod_syr,mod_nyr,sage,nage);
  LOCAL_CALCS
    int iyr = mod_syr;
    
    for(int h = 1; h <= nFecBlocks; h++){
      do{
        Eij(iyr) = 1.e-6 *
                (data_sp_waa(iyr)(sage,nage) * fec_slope(h) - fec_inter(h));
        iyr ++;
      }while(iyr <= nFecBlockYears(h));
    }
  END_CALCS

// |---------------------------------------------------------------------------|
// | END OF DATA FILE
// |---------------------------------------------------------------------------|
  init_int dat_eof; 
  !! if(dat_eof != 999){cout<<"Error reading data file, aborting."<<endl; exit(1);}
  

// BUG TEEST


// |---------------------------------------------------------------------------|
// | READ CONTENTS OF CONTROL FILE
// |---------------------------------------------------------------------------|
  !! ad_comm::change_datafile_name(ControlFile);


// |-------------------------------------------------------------------------|
// | DESIGN MATRIX FOR PARAMETER CONTROLS                                    |
// |-------------------------------------------------------------------------|
// | - theta_DM -> theta is a vector of estimated parameters.
  int n_theta;
  !! n_theta = 6;
  init_matrix theta_DM(1,n_theta,1,7);
  vector    theta_ival(1,n_theta);
  vector      theta_lb(1,n_theta);
  vector      theta_ub(1,n_theta);
  ivector    theta_phz(1,n_theta);
  ivector theta_iprior(1,n_theta);
  vector      theta_p1(1,n_theta);
  vector      theta_p2(1,n_theta);
  !! theta_ival = column(theta_DM,1);
  !! theta_lb  = column(theta_DM,2);
  !! theta_ub  = column(theta_DM,3);
  !! theta_phz = ivector(column(theta_DM,4));
  !! theta_iprior = ivector(column(theta_DM,5));
  !! theta_p1 = column(theta_DM,6);
  !! theta_p2 = column(theta_DM,7);
  

// |---------------------------------------------------------------------------|
// | Controls for time-varying maturity
// |---------------------------------------------------------------------------|
// | Oct 12, 2016:
// |  - SJDM Added option to specify fixed maturity for each block.
// | 
  init_int nMatBlocks;
  init_matrix maturity_cont(1,nMatBlocks,1,4);
  vector mat_a50;
  vector mat_a95;
  ivector mat_phz(1,nMatBlocks);
  ivector nMatBlockYear(1,nMatBlocks);
  LOCAL_CALCS
    mat_a50 = column(maturity_cont,1);
    mat_a95 = column(maturity_cont,2);
    mat_phz = ivector(column(maturity_cont,3));
    nMatBlockYear = ivector(column(maturity_cont,4));
  END_CALCS


// |---------------------------------------------------------------------------|
// | Controls for natural mortality rate deviations in each block.
// |---------------------------------------------------------------------------|
  init_int mort_type;
  init_int mort_dev_phz;
  init_int nMortBlocks;
  init_ivector nMortBlockYear(1,nMortBlocks);
  

// |---------------------------------------------------------------------------|
// | Controls for selectivity parameters
// |---------------------------------------------------------------------------|
// | - nSlxCols     » number of columns in selectivity design matrix
// | - nSlxBlks       » number of selectivity blocks/patterns 
// | - selex_cont     » matrix of controls to be read in from control file.
  int nSlxCols;
  !! nSlxCols = 9;
  init_int nSlxBlks;
  init_matrix selex_cont(1,nSlxBlks,1,nSlxCols);
  ivector       nSelType(1,nSlxBlks);
  ivector       nslx_phz(1,nSlxBlks);
  ivector      nslx_rows(1,nSlxBlks);
  ivector      nslx_cols(1,nSlxBlks);
  ivector       nslx_syr(1,nSlxBlks);
  ivector       nslx_nyr(1,nSlxBlks);


  LOCAL_CALCS
    nSelType = ivector(column(selex_cont,2));
    nslx_phz = ivector(column(selex_cont,7));
    nslx_syr = ivector(column(selex_cont,8));
    nslx_nyr = ivector(column(selex_cont,9));

    // determine dimensions for log_slx_pars ragged object.
    for(int h = 1; h <= nSlxBlks; h++){
      nslx_rows(h) = 1;
      switch(nSelType(h)){
        case 1: // logistic 2-parameters
          nslx_cols = int(2);
        break;
      }
    }
  END_CALCS

// |---------------------------------------------------------------------------|
// | Miscellaneous Controls
// |---------------------------------------------------------------------------|
// | nMiscCont  » Number of controls to read in.
// | dMiscCont  » Vector of miscelaneous controls,
// | pos 1 » Catch scaler.
// | pos 2 » Condition on catch, condition on effort.

  init_int nMiscCont;
  init_vector dMiscCont(1,nMiscCont);

  // Rescale the catch data as needed.
  matrix  data_catch(dat_syr,dat_nyr,1,3);
  LOCAL_CALCS
    data_catch = data_ct_raw;
    for( int i = dat_syr; i <= dat_nyr; i++ ) {
      data_catch(i,2) = dMiscCont(1) * data_ct_raw(i,2);
    }
    if(dMiscCont(2)) cout<<"Condition model on Ft"<<endl;
  END_CALCS


// |---------------------------------------------------------------------------|
// | END OF Control FILE
// |---------------------------------------------------------------------------|
  init_int ctl_eof;
  LOCAL_CALCS
    if(ctl_eof != 999){
      cout<<"Error reading control file, aborting."<<ctl_eof<<endl; 
      exit(1);
    }
  END_CALCS

// |---------------------------------------------------------------------------|
// | RETROSPECTIVE ADJUSTMENTS
// |---------------------------------------------------------------------------|
  !! mod_nyr = mod_nyr - retro_yrs;
  int nf;

  !! nf = 0;

// |---------------------------------------------------------------------------|
// | VARIABLES FOR STORING SIMULATED VALUES
// |---------------------------------------------------------------------------|
  vector sim_spawners(rec_syr,mod_nyr+1);
  vector sim_recruits(rec_syr,mod_nyr+1);

INITIALIZATION_SECTION
  theta theta_ival;


PARAMETER_SECTION

// |---------------------------------------------------------------------------|
// | POPULATION PARAMETERS
// |---------------------------------------------------------------------------|
// | - theta(1) -> log natural mortality
// | - theta(2) -> log initial average age-3 recruitment for ages 4-9+ in dat_styr
// | - theta(3) -> log average age-3 recruitment from dat_styr to dat_endyr
// | - theta(4) -> log of unfished recruitment.
// | - theta(5) -> log of recruitment compensation (reck > 1.0)
// | - theta(6) -> log of simga R
  init_bounded_number_vector theta(1,n_theta,theta_lb,theta_ub,theta_phz);
  number log_natural_mortality;
  number log_rinit;
  number log_rbar;
  number log_ro;
  number log_reck;
  number log_sigma_r;
  init_bounded_dev_vector log_rinit_devs(sage+1,nage,-15.0,15.0,2);
  init_bounded_dev_vector log_rbar_devs(mod_syr,mod_nyr+1,-15.0,15.0,2);
  LOCAL_CALCS
    if( global_parfile ) {
      theta_ival = value(theta);
    }
  END_CALCS
  

// |---------------------------------------------------------------------------|
// | MATURITY PARAMETERS
// |---------------------------------------------------------------------------|
// | TO BE DEPRECATED
// | - mat_params[1] -> Age at 50% maturity
// | - mat_params[2] -> Slope at 50% maturity
  init_bounded_vector_vector mat_params(1,nMatBlocks,1,2,0,100,mat_phz);
  //init_bounded_matrix mat_params(1,nMatBlocks,1,2,0,10,mat_phz);
  matrix mat(mod_syr,mod_nyr,sage,nage);
  LOCAL_CALCS
    //cout<<"Good to here"<<endl;
    //cout<<mat_params(1)<<endl;
    if( !global_parfile ) {
      for(int h = 1; h <= nMatBlocks; h++){
        mat_params(h,1) = mat_a50(h);
        mat_params(h,2) = mat_a95(h);       
      }
    }
    //cout<<mat_params(1)<<endl;
    
  END_CALCS

// |---------------------------------------------------------------------------|
// | NATURAL MORTALITY PARAMETERS
// |---------------------------------------------------------------------------|
// | - log_m_devs   -> deviations in natural mortality for each block.
// | - Mij          -> Array for natural mortality rate by year and age.
  init_bounded_dev_vector log_m_devs(1,nMortBlocks,-15.0,15.0,mort_dev_phz);
  matrix Mij(mod_syr,mod_nyr,sage,nage);

// |---------------------------------------------------------------------------|
// | SELECTIVITY PARAMETERS
// |---------------------------------------------------------------------------|
// | - log_slx_pars » parameters for selectivity models (ragged object).
  init_bounded_matrix_vector log_slx_pars(1,nSlxBlks,1,nslx_rows,1,nslx_cols,-25,25,nslx_phz);
  matrix log_slx(mod_syr,mod_nyr,sage,nage);
  LOCAL_CALCS
   if( ! global_parfile ){
     for(int h = 1; h <= nSlxBlks; h++){
       switch(nSelType(h)){
         case 1: //logistic
           log_slx_pars(h,1,1) = log(selex_cont(h,3));
           log_slx_pars(h,1,2) = log(selex_cont(h,4));
         break; 
       }
     } 
   }
  END_CALCS

// |---------------------------------------------------------------------------|
// | FISHING MORTALITY RATE PARAMETERS
// |---------------------------------------------------------------------------|
// |
  !! int phz; phz = dMiscCont(2)==0?-1:1;
  init_bounded_vector log_ft_pars(mod_syr,mod_nyr,-30.,3.0,phz);
  LOCAL_CALCS
    if(b_simulation_flag && ! global_parfile) log_ft_pars = log(0.2);
  END_CALCS

// |---------------------------------------------------------------------------|
// | VARIABLES
// |---------------------------------------------------------------------------|
// |- fore_sb spawning biomass
// |- fore_sb vulnerable biomass
// |- ghl guidline harvest level
  number ro;  
  number reck;
  number so;  
  number beta;
  number fore_sb;   
  number fore_vb;   
  number ghl;       


// |---------------------------------------------------------------------------|
// | VECTORS
// |---------------------------------------------------------------------------|
// | - ssb      » spawning stock biomass at the time of spawning.
// | - recruits » vector of sage recruits predicted by S-R curve.
// | - spawners » vector of ssb indexed by brood year.
// | - resd_rec » vector of residual process error (log-normal).
  vector ssb(mod_syr,mod_nyr);
  vector recruits(rec_syr,mod_nyr+1);
  vector spawners(rec_syr,mod_nyr+1);
  vector resd_rec(rec_syr,mod_nyr+1);

  vector pred_egg_dep(mod_syr,mod_nyr);
  vector resd_egg_dep(mod_syr,mod_nyr);

  vector pred_mileday(mod_syr,mod_nyr);
  vector resd_mileday(mod_syr,mod_nyr);

  vector pred_catch(mod_syr,mod_nyr);
  vector resd_catch(mod_syr,mod_nyr);

// |---------------------------------------------------------------------------|
// | MATRIXES
// |---------------------------------------------------------------------------|
// | - Nij      » numbers-at-age N(syr,nyr,sage,nage)
// | - Oij      » mature numbers-at-age O(syr,nyr,sage,nage)
// | - Pij      » numbers-at-age P(syr,nyr,sage,nage) post harvest.
// | - Sij      » selectivity-at-age 
// | - Qij      » vulnerable proportions-at-age
// | - Cij      » predicted catch-at-age in numbers.
  matrix Nij(mod_syr,mod_nyr+1,sage,nage);
  matrix Oij(mod_syr,mod_nyr+1,sage,nage);
  matrix Pij(mod_syr,mod_nyr+1,sage,nage);
  matrix Sij(mod_syr,mod_nyr+1,sage,nage);
  matrix Qij(mod_syr,mod_nyr+1,sage,nage);
  matrix Cij(mod_syr,mod_nyr+1,sage,nage);
  matrix Fij(mod_syr,mod_nyr+1,sage,nage);

  matrix pred_cm_comp(mod_syr,mod_nyr,sage,nage);
  matrix resd_cm_comp(mod_syr,mod_nyr,sage,nage);
  matrix pred_sp_comp(mod_syr,mod_nyr,sage,nage);
  matrix resd_sp_comp(mod_syr,mod_nyr,sage,nage);

// |---------------------------------------------------------------------------|
// | OBJECTIVE FUNCTION VALUE
// |---------------------------------------------------------------------------|
  vector nll(1,7);
  objective_function_value f;

  number fpen;
  sdreport_number sd_terminal_ssb;
  sdreport_number sd_forecast_ssb;
  sdreport_number sd_projected_ssb;
  sdreport_vector sd_ssb(mod_syr,mod_nyr);


PRELIMINARY_CALCS_SECTION

  /* 
   * SIMULATION MODEL SWITCH
   */
  if( b_simulation_flag && rseed >= 0) {
    cout<<"|--------------------------|"<<endl;
    cout<<"| RUNNING SIMULATION MODEL |"<<endl;
    cout<<"|--------------------------|"<<endl;
    
    char type;
    do
    {
        cout<<"Theta \n"<<theta<<endl;
        cout<<"| Continue? [y]es or [n]o "<<endl;
        cin >> type;
    }
    while( !cin.fail() && type!='y' && type!='n' );
    if( type =='y' ){
      runSimulationModel(rseed);
    } else {
      exit(1);
    }

  } else if ( b_simulation_flag && rseed < 0 ){
    theta_ival = value(theta);
    runSimulationModel(rseed);
  }


PROCEDURE_SECTION
 
// |---------------------------------------------------------------------------|
// | RUN STOCK ASSEAAMENT MODEL ROUTINES
// |---------------------------------------------------------------------------|
// | PSUEDOCODE:
// | - initialize model parameters.
// | - initialize Maturity Schedule information.
// | - get natural mortality schedules.
// | - get fisheries selectivity schedules.
// | - initialize State variables
// | - update State variables:
// |    - calculate spawning stock biomass
// |    - calculate age-composition residuals
// | - stock-recruitment relationship
// | - observation models:
// |    - calculate age-composition residuals
// |    - calculate egg-deposition survey residuals
// |    - calculate mile-days of milt residuals
// |    - calculate catch residuals if model is conditioned on effort
// | - calculate objective function value
// | - TERMINAL PHASE:
// |    - run forecast model
// |    - save posterior samples from -mceval cmd line option
// |---------------------------------------------------------------------------|
  
  initializeModelParameters();
  if(DEBUG_FLAG) cout<<"--> Ok after initializeModelParameters      <--"<<endl;

  initializeMaturitySchedules();
  if(DEBUG_FLAG) cout<<"--> Ok after initializeMaturitySchedules    <--"<<endl;

  calcNaturalMortality();
  if(DEBUG_FLAG) cout<<"--> Ok after calcNaturalMortality           <--"<<endl;
  
  calcSelectivity();
  if(DEBUG_FLAG) cout<<"--> Ok after calcSelectivity                <--"<<endl;
  
  if( dMiscCont(2) ) {
    calcFishingMortalitiy();
    if(DEBUG_FLAG) cout<<"--> Ok after calcFishingMortalitiy          <--"<<endl;  
  }

  initializeStateVariables();
  if(DEBUG_FLAG) cout<<"--> Ok after initializeStateVariables       <--"<<endl;
  
  updateStateVariables();
  if(DEBUG_FLAG) cout<<"--> Ok after updateStateVariables           <--"<<endl;
  
  calcSpawningStockRecruitment();
  if(DEBUG_FLAG) cout<<"--> Ok after calcSpawningStockRecruitment   <--"<<endl;

  calcAgeCompResiduals();
  if(DEBUG_FLAG) cout<<"--> Ok after calcAgeCompResiduals           <--"<<endl;

  calcEggSurveyResiduals();
  if(DEBUG_FLAG) cout<<"--> Ok after calcEggSurveyResiduals         <--"<<endl;

  calcMiledaySurveyResiduals();
  if(DEBUG_FLAG) cout<<"--> Ok after calcMiledaySurveyResiduals     <--"<<endl;

  //if( dMiscCont(2) ) {
    calcCatchResiduals();
    if(DEBUG_FLAG) cout<<"--> Ok after calcCatchResiduals           <--"<<endl; 
  // }

  calcObjectiveFunction(); nf++;
  if(DEBUG_FLAG) cout<<"--> Ok after calcObjectiveFunction          <--"<<endl;


  if( last_phase() ) {
    runForecast();
  }

  if(mceval_phase()) {
    writePosteriorSamples();
  }
  if(sd_phase()) {
    sd_terminal_ssb = ssb(mod_nyr);
    sd_ssb  = ssb(mod_syr,mod_nyr);
  }

// |---------------------------------------------------------------------------|

FUNCTION void runForecast()
  /** 
  Conduct a 1-year ahead forcasts based on SR
  PSUEDOCODE:
    -1 declare variables for forcasting.
      * recruitment, spawning biomass, numbers-at-age
      * selectivity curve
    -2 Calculate F-at-age conditional on harvest rule.
      * harvest_rate = 20% || set TAC option
    -3 Update state variables from pyr=mod_nyr+1 to pyr+2
      * recruitment based on ssb(pyr-sage)
    -4 Compute GHLs given threshold and target harvest rate
      * user specifies threshold and harvest rate in control file.

  **/
  int nyr = mod_nyr;
  int pyr = nyr+1;
  dvariable fore_rt;  // sage recruits

  dvar_vector fore_nj(sage,nage); //numbers-at-age
  dvar_vector fore_cj(sage,nage); //catch-at-age

  fore_rt = so * ssb(pyr-sage) * exp(-beta*ssb(pyr-sage));
  fore_nj = Nij(pyr); fore_nj(sage) = fore_rt;
  fore_vb = fore_nj * elem_prod(Sij(nyr),data_cm_waa(nyr)(sage,nage));
  fore_sb = fore_nj * elem_prod(mat(nyr),data_sp_waa(nyr)(sage,nage));
  sd_forecast_ssb = fore_sb;

  // GHL for pyr
  double ssb_threshold = dMiscCont(3);
  double target_hr = dMiscCont(4);
  dvariable hr = (2.0 + 8.0*fore_sb / dMiscCont(5))/100.0;
  if( hr > target_hr) {
    hr = target_hr;
  } else if( fore_sb < ssb_threshold ) {
    hr = 0.0;
  }

  ghl = hr * fore_sb;
  //cout<<"harvest rate = "<<hr<<" GHL = "<<ghl<<endl;

  // update state variables to pyr+1 so you can predict
  // the effect of the 2016 fishery on the 2017 spawning stock.
  // predicted catch-at-age
  dvar_vector pa = elem_prod(fore_nj,Sij(nyr));
  pa /= sum(pa);
  dvariable wbar = pa * data_cm_waa(nyr)(sage,nage);
  fore_cj = ghl/wbar * pa;
  fore_nj = elem_prod(fore_nj - fore_cj,mfexp(-Mij(nyr)));
  fore_nj(sage) = so * ssb(pyr+1-sage) * exp(-beta*ssb(pyr-sage));
  fore_sb = fore_nj * elem_prod(mat(nyr),data_sp_waa(nyr)(sage,nage));
  sd_projected_ssb = fore_sb;

FUNCTION void writePosteriorSamples()
  /**
  - This function is only envoked when the -mceval
    command line option is implemented.
  */
  if(nf==1){
    ofstream ofs("ssb.ps");
  }
  ofstream ofs("ssb.ps",ios::app);
  ofs<<ssb<<endl;
  
FUNCTION void runSimulationModel(const int& rseed)
  /*
    PSUEDOCODE:
    1) initialize model parameters based on pin file, or control file.
    2) initialize maturity schedules for all block years.
    3) initialize natural mortality schedules for all block years.
    4) calculate selectivity parameters.
    5) generate random normal deviates for process/observation errors.
    6) initialize state variables.
    7) update state variables conditioned on the observed catch data.
    8) calculate age-composition residuals and over-write input data in memory.
    9) calculate egg survey residuals and simulate fake survey.
  */
  if(global_parfile) {
    cout<<"\nUsing pin file for simulation parameter values.\n"<<endl;
  }

  // 1) initialize model parameters based on pin file, or control file.
  initializeModelParameters();

  // 2) initialize maturity schedules for all block years.
  initializeMaturitySchedules();

  // 3) initialize natural mortality schedules for all block years.
  calcNaturalMortality();

  // 4) calculate selectivity parameters.
  calcSelectivity();
  
  // 4b) calculate fishing mortality.
  calcFishingMortalitiy();

  // 5) generate random normal deviates for process errors:
  //    - log_m_devs
  //    - log_rinit_devs
  //    - log_rbar_devs
  random_number_generator rng(rseed);
  dvector epsilon_m_devs(1,nMortBlocks);
  dvector epsilon_rbar_devs(mod_syr,mod_nyr+1);
  dvector epsilon_rinit_devs(sage+1,nage);

  double sigma_m_devs = dMiscCont(6);
  double sigma_rbar_devs = value(exp(log_sigma_r));
  double sigma_rinit_devs = value(exp(log_sigma_r));

  if ( rseed == 0 ) {
    sigma_m_devs = 0.0;
    sigma_rbar_devs = 0.0;
    sigma_rinit_devs = 0.0;

  }

  epsilon_m_devs.fill_randn(rng);
  epsilon_rbar_devs.fill_randn(rng);
  epsilon_rinit_devs.fill_randn(rng);

  log_m_devs = dvar_vector(epsilon_m_devs * sigma_m_devs - 0.5 * square(sigma_m_devs));
  log_rbar_devs = dvar_vector(epsilon_rbar_devs * sigma_rbar_devs - 0.5 * square(sigma_rbar_devs));
  log_rinit_devs = dvar_vector(epsilon_rinit_devs * sigma_rinit_devs - 0.5 * square(sigma_rinit_devs));

  
  // 6) initialize state variables.
  initializeStateVariables();

  // 7) update state variables conditioned on the observed catch data.
  updateStateVariables();

  // 7a) calcSpawningRecruitment
  calcSpawningStockRecruitment();
  sim_spawners = value(spawners);
  sim_recruits = value(column(Nij,sage)(rec_syr,mod_nyr+1));


  // 8) calculate age-composition residuals and over-write input data in memory.
  //    - Note the age-comps are sampled from a multivariate logistic dist.
  calcAgeCompResiduals();
  double sigma_a = 0.07;
  if( rseed == 0 ) sigma_a = 0.0;
  for(int i = mod_syr; i <= mod_nyr; i++) {
    dvector t1 = value(pred_cm_comp(i));
    dvector t2 = value(pred_sp_comp(i));
    
    data_cm_comp(i)(sage,nage) = rmvlogistic(t1,sigma_a,rseed + i);
    data_sp_comp(i)(sage,nage) = rmvlogistic(t2,sigma_a,rseed - i);
  }


  // 9) calculate egg survey residuals and simulate fake survey.
  calcEggSurveyResiduals();
  
  dvector epsilon_egg_dep(mod_syr,mod_nyr);
  epsilon_egg_dep.fill_randn(rng);
  for(int i = mod_syr; i <= mod_nyr; i++) {
    if( data_egg_dep(i,2) > 0 ) {
      double sd = data_egg_dep(i,3);
      //double log_it = value(log(pred_egg_dep(i)));
      data_egg_dep(i,2) = value(pred_egg_dep(i)) 
                          * exp(epsilon_egg_dep(i) * sd - 0.5 * square(sd));  
      //data_egg_dep(i,2) = exp(log_it + epsilon_egg_dep(i) * sd) - 05 * square(sd));
    }
  }

  // 10) calculate mile_days
  calcMiledaySurveyResiduals();
  
  dvector epsilon_mileday(mod_syr,mod_nyr);
  epsilon_mileday.fill_randn(rng);
  for(int i = mod_syr; i <= mod_nyr; i++) {
    if( data_mileday(i,2) > 0 ) {
      double sd = data_mileday(i,3);
      data_mileday(i,2) = value(pred_mileday(i))
                          * exp(epsilon_mileday(i) * sd - 0.5 * square(sd));
    }
  }

  // 11) calculate predicted catch
  calcCatchResiduals();

  dvector epsilon_catch(mod_syr,mod_nyr);
  epsilon_catch.fill_randn(rng);
  for(int i = mod_syr; i <= mod_nyr; i++) {
    if( data_ct_raw(i,2) > 0 ) {
      double sd = data_ct_raw(i,3);
      data_ct_raw(i,2) = value(pred_catch(i))
                          * exp(epsilon_catch(i) * sd - 0.5 * square(sd));
    }
  }


FUNCTION void initializeModelParameters()
  fpen = 0;
  log_natural_mortality = theta(1);
  log_rinit             = theta(2);
  log_rbar              = theta(3);
  log_ro                = theta(4);
  log_reck              = theta(5);
  log_sigma_r           = sqrt(1.0/theta(6));
  // COUT(theta);



FUNCTION void initializeMaturitySchedules() 
  int iyr = mod_syr;
  int jyr = 0;
  mat.initialize();
  
  for(int h = 1; h <= nMatBlocks; h++) {
    dvariable mat_a50 = mat_params(h,1);
    dvariable mat_a95 = mat_params(h,2);

    jyr = (h != nMatBlocks) ? nMatBlockYear(h) : nMatBlockYear(h)-retro_yrs;
    
    // fill maturity array using logistic function
    do{
      mat(iyr++) = plogis95(age,mat_a50,mat_a95);
    } while(iyr <= jyr);  
  }
  

FUNCTION void calcNaturalMortality()
  
  int iyr = mod_syr;
  int jyr;
  Mij.initialize();
  switch(mort_type) {
    case 1: // constant M within block.
      for(int h = 1; h <= nMortBlocks; h++){
        dvariable mi = mfexp(log_natural_mortality + log_m_devs(h));
        
        jyr = h != nMortBlocks?nMortBlockYear(h):nMortBlockYear(h)-retro_yrs;
        // fill mortality array by block
        do{
          Mij(iyr++) = mi;
        } while(iyr <= jyr);
      }
    break;
    case 2: // cubic spline  
      dvector iiyr = dvector((nMortBlockYear - mod_syr)/(mod_nyr-mod_syr));
      dvector jjyr(mod_syr,mod_nyr);
      jjyr.fill_seqadd(0,double(1.0/(mod_nyr-mod_syr)));
      dvar_vector mi = log_natural_mortality + log_m_devs;
      vcubic_spline_function cubic_spline_m(iiyr,mi);
      dvar_vector mtmp = cubic_spline_m(jjyr);
      for(int i=mod_syr; i <= mod_nyr; i++)
      {
        Mij(i) = mfexp(mtmp(i));
      }
    break;
  }
  
  

FUNCTION void calcSelectivity()
  /**
    - Loop over each of the selectivity block/pattern
    - Determine which selectivity type is being used.
    - get parameters from log_slx_pars
    - calculate the age-dependent selectivity pattern
    - fill selectivty array for that block.
    - selectivity is scaled to have a mean = 1 across all ages.
  */
  dvariable p1,p2;
  dvar_vector slx(sage,nage);
  log_slx.initialize();
  
  

  for(int h = 1; h <= nSlxBlks; h++){

    switch(nSelType(h)){
      case 1: //logistic
        p1  = mfexp(log_slx_pars(h,1,1));
        p2  = mfexp(log_slx_pars(h,1,2));
        slx = plogis(age,p1,p2) + TINY;
      break;
    }
    
    int jyr = h != nSlxBlks ? nslx_nyr(h):nslx_nyr(h)-retro_yrs;
    for(int i = nslx_syr(h); i <= jyr; i++){
      log_slx(i) = log(slx) - log(mean(slx));
    }
  }
  Sij.sub(mod_syr,mod_nyr) = mfexp(log_slx);


FUNCTION void calcFishingMortalitiy()
  /**
    - Calculate Fishing mortality, and then Zij = Mij + Fij
    */

  for(int i = mod_syr; i <= mod_nyr; i++) {
    Fij(i) = exp(log_ft_pars(i)) * Sij(i);
  }



FUNCTION void initializeStateVariables()
  /**
    - Set initial values for numbers-at-age matrix in first year
      and sage recruits for all years.
    */
  Nij.initialize();

  // initialize first row of numbers-at-age matrix
  // lx is a vector of survivorship (probability of surviving to age j)
  dvar_vector lx(sage,nage);

  for(int j = sage; j <= nage; j++){
    
    lx(j) = exp(-Mij(mod_syr,j)*(j-sage));
    if( j==nage ) lx(j) /= (1.0-exp(-Mij(mod_syr,j)));


    if( j > sage ){
      Nij(mod_syr)(j) = mfexp(log_rinit + log_rinit_devs(j)) * lx(j);     
    }
  } 
    


  // iniitialize first column of numbers-at-age matrix
  for(int i = mod_syr; i <= mod_nyr + 1; i++){
    Nij(i,sage) = mfexp(log_rbar + log_rbar_devs(i));
  }
  //COUT(lx);
  //COUT(Nij);


FUNCTION void updateStateVariables()
  /**
  - Update the numbers-at-age conditional on the catch-at-age.
  - Assume a pulse fishery.
  - step 1 calculate a vector of vulnerable-numbers-at-age
  - step 2 calculate vulnerable proportions-at-age.
  - step 3 calc average weight of catch (wbar) conditional on Qij.
  - step 4 calc catch-at-age | catch in biomass Cij = Ct/wbar * Qij.
  - step 5 condition on Ft or else condition on observed catch.
  - step 6 update numbers-at-age (using a very dangerous difference eqn.)

  Nov 30, 2016
  - added options for simulation model:
    - b_simulation_flag is true, then condition model on catch & S-R curve.

  */

  Qij.initialize();
  Cij.initialize();
  Pij.initialize();
  dvariable wbar;   // average weight of the catch.
  dvar_vector vj(sage,nage);
  //dvar_vector pj(sage,nage);
  dvar_vector sj(sage,nage);
  

  for(int i = mod_syr; i <= mod_nyr; i++){


    // step 1.
    vj = elem_prod(Nij(i),Sij(i));

    // step 2.
    Qij(i) = vj / sum(vj);

    // ADF&G's approach.
    if( !dMiscCont(2) ) {
      // step 3.
      dvector wa = data_cm_waa(i)(sage,nage);
      wbar = wa * Qij(i);
    
      // step 4.
      Cij(i) = data_catch(i,2) / wbar * Qij(i);
      // should use posfun here
      Pij(i) = posfun(Nij(i) - Cij(i),0.01,fpen); 
      
      // step 6. update numbers at age    
      sj = mfexp(-Mij(i));
      Nij(i+1)(sage+1,nage) =++ elem_prod(Pij(i)(sage,nage-1),
                                          sj(sage,nage-1));
      Nij(i+1)(nage) += Pij(i,nage) * sj(nage);       
    } 
    // step 5.
    // Condition on Ft
    else { 
      // add flexibility here for Popes approx, or different seasons               
      Pij(i) = elem_prod( Nij(i), exp(-Fij(i)) );
      Cij(i) = elem_prod( Nij(i), 1.-exp(-Fij(i)));
      
      // the following assumes fishery is at the start of the year.
      dvar_vector zi = Mij(i) + Fij(i);
      Nij(i+1)(sage+1,nage) = ++ elem_prod(Nij(i)(sage,nage-1),
                                           mfexp(-zi(sage,nage-1)));
      Nij(i+1)(nage) += Nij(i,nage) * mfexp(-zi(nage));
    }
  }

  // cross check... Looks good.
  // COUT(Cij(mod_syr) * data_cm_waa(mod_syr)(sage,nage));


FUNCTION void calcSpawningStockRecruitment()
  /**
  The functional form of the stock recruitment model follows that of a 
  Ricker model, where R = so * SSB * exp(-beta * SSB).  The two parameters
  so and beta where previously estimated as free parameters in the old
  herring model.  Herein this fucntion I derive so and beta from the 
  leading parameters Ro and reck; Ro is the unfished sage recruits, and reck
  is the recruitment compensation parameter, or the relative improvement in
  juvenile survival rates as the spawning stock SSB tends to 0. Simply a 
  multiple of the replacement line Ro/Bo.

  At issue here is time varying maturity and time-varying natural mortality.
  When either of these two variables are assumed to change over time, then
  the underlying stock recruitment relationship will also change. This 
  results in a non-stationary distribution.  For the purposes of this 
  assessment model, I use the average mortality and maturity schedules to
  derive the spawning boimass per recruit, which is ultimately used in 
  deriving the parameters for the stock recruitment relationship.
    */

  /*
    Spoke to Sherri about this. Agreed to change the equation to prevent
    substracting immature fish from the numbers before calculating SSB.
  */
  for(int i = mod_syr; i <= mod_nyr; i++){
    //Oij(i) = elem_prod(mat(i),Nij(i));
    //ssb(i) = (Oij(i) - Cij(i)) * data_sp_waa(i)(sage,nage);

    Oij(i) = elem_prod(mat(i),Nij(i)-Cij(i));
    ssb(i) = Oij(i) * data_sp_waa(i)(sage,nage);
  }
  

  // average natural mortality
  dvar_vector mbar(sage,nage);
  int n = Mij.rowmax() - Mij.rowmin() + 1;
  mbar  = colsum(Mij)/n;

  // average maturity
  dvar_vector mat_bar(sage,nage);
  mat_bar = colsum(mat)/n;


  // unfished spawning biomass per recruit
  dvar_vector lx(sage,nage);
  lx(sage) = 1.0;
  for(int j = sage + 1; j <= nage; j++){
    lx(j) = lx(j-1) * mfexp(-mbar(j-1));
    if(j == nage){
      lx(j) /= 1.0 - mfexp(-mbar(j));
    }
  }
  dvariable phie = lx * elem_prod(avg_sp_waa,mat_bar);

  // Ricker stock-recruitment function 
  // so = reck/phiE; where reck > 1.0
  // beta = log(reck)/(ro * phiE)

  // Beverton Holt use:
  // beta = (reck - 1.0)/(ro *phiE)
  ro   = mfexp(log_ro);
  reck = mfexp(log_reck) + 1.0;
  so   = reck/phie;
  beta = log(reck) / (ro * phie);

  spawners = ssb(mod_syr,mod_nyr-sage+1).shift(rec_syr);
  recruits = elem_prod( so*spawners , mfexp(-beta*spawners) );
  resd_rec = log(column(Nij,sage)(rec_syr,mod_nyr+1)+TINY) 
             - log(recruits+TINY);



FUNCTION void calcAgeCompResiduals()
  /**
  - Commercial catch-age comp residuals
  - Spawning survey catch-age comp residuals.
  */

  resd_cm_comp.initialize();
  resd_sp_comp.initialize();
  for(int i = mod_syr; i <= mod_nyr; i++){
    
    // commercial age-comp prediction 
    pred_cm_comp(i) = Qij(i);
    if( data_cm_comp(i,sage) >= 0 ){
      resd_cm_comp(i) = data_cm_comp(i)(sage,nage) - pred_cm_comp(i);
    }

    // spawning age-comp prediction
    pred_sp_comp(i) = (Oij(i)+TINY) / sum(Oij(i)+TINY);
    if( data_sp_comp(i,sage) >= 0 ){
      resd_sp_comp(i) = data_sp_comp(i)(sage,nage) - pred_sp_comp(i);
    }
  }

  //COUT(resd_sp_comp);

FUNCTION void calcEggSurveyResiduals()
  /**
  - Observed egg data is in trillions of eggs
  - Predicted eggs is the mature female numbers-at-age multiplied 
    by the fecundity-at-age, which comes from a regession of 
    fecundity = slope * obs_sp_waa - intercept
  - Note Eij is the Fecundity-at-age j in year i.
  - 
  */
  resd_egg_dep.initialize();
  for(int i = mod_syr; i <= mod_nyr; i++){
    pred_egg_dep(i) = (0.5 * Oij(i)) * Eij(i);
    if(data_egg_dep(i,2) > 0){
      resd_egg_dep(i) = log(data_egg_dep(i,2)) - log(pred_egg_dep(i));
    }
  }

FUNCTION void calcMiledaySurveyResiduals()
  /**
  - Assumed index from aerial survey is a relative abundance
  index. The slope of the regression ln(SSB) = q * ln(MileMilt) + 0
  is computed on the fly in case of missing data. 
  - See Walters and Ludwig 1994 for more details.
  */
  resd_mileday.initialize();
  pred_mileday.initialize();
  int n = 1;
  dvar_vector zt(mod_syr,mod_nyr); zt.initialize();
  dvariable zbar = 0;
  for(int i = mod_syr; i <= mod_nyr; i++){
    if(data_mileday(i,2) > 0){
      zt(i) = log(data_mileday(i,2)) - log(ssb(i));
      zbar  = zt(i)/n + zbar *(n-1)/n;
      n++ ;
    }   
  }
  pred_mileday = ssb * exp(zbar);
  resd_mileday = zt - zbar;
  // COUT(resd_mileday);

FUNCTION void calcCatchResiduals()
  /**
  - Catch residuals assuming a lognormal error structure.
  */
  pred_catch.initialize();
  resd_catch.initialize();
  for(int i = mod_syr; i <= mod_nyr; i++) {
    if(data_catch(i,2) > 0) {
      pred_catch(i) = Cij(i) * data_cm_waa(i)(sage,nage);
      resd_catch(i) = log(data_catch(i,2)) - log(pred_catch(i));
    }
  }

FUNCTION void calcObjectiveFunction()
  /**
  - THIS FUNCTION IS ORGANIZED AS FOLLOWS:
    1. Negative logliklihoods (nll)
    2. Penalized loglikelihoods (penll)

  */



  // 1. Negative loglikelihoods
  nll.initialize();

  // Mulitvariate logistic likelihood for composition data.
  double sp_tau2;
  double minp = 0.00;
  int t1 = mod_syr;
  int t2 = mod_nyr;
  dmatrix d_sp_comp = trans(trans(data_sp_comp).sub(sage,nage)).sub(t1,t2);
  nll(1) = dmvlogistic(d_sp_comp,pred_sp_comp,resd_sp_comp,sp_tau2,minp);
  
  // Mulitvariate logistic likelihood for composition data.
  double cm_tau2;
  dmatrix d_cm_comp  = trans(trans(data_cm_comp).sub(sage,nage)).sub(t1,t2);
  nll(2) = dmvlogistic(d_cm_comp,pred_cm_comp,resd_cm_comp,cm_tau2,minp);

  // Negative loglikelihood for egg deposition data
  dvector std_egg_dep = TINY + column(data_egg_dep,3)(t1,t2);
  nll(3) = dnorm(resd_egg_dep,std_egg_dep);

  // Negative loglikelihood for milt mile day
  dvector std_mileday = TINY + column(data_mileday,3)(t1,t2);
  nll(4) = dnorm(resd_mileday,std_mileday);

  // Negative loglikelihood for stock-recruitment data
  dvariable std_rec = log_sigma_r;
  nll(5) = dnorm(resd_rec,std_rec);
  
  // Negative loglikelihood for catch data
  dvector std_catch = column(data_catch,3)(t1,t2);
  nll(6) = dnorm(resd_catch,std_catch);

  



  // 2. Penalized logliklelihoods

  
  dvar_vector penll(1,4);
  penll.initialize();

  // | Terminal phase of estimation.
  dvariable log_fbar = mean(log_ft_pars);
  if( last_phase() ){
    penll(1) = dnorm(log_rinit_devs,0.0,5.0);
    penll(2) = dnorm(log_rbar_devs,0.0,5.0);
    penll(3) = dnorm(log_fbar,0.2,2.0);

  } else {
    penll(1) = dnorm(log_rinit_devs,0.0,1.0);
    penll(2) = dnorm(log_rbar_devs,0.0,1.0);
    penll(3) = dnorm(log_fbar,0.2,0.07); 
  }

  if( active(log_m_devs) ) {
    penll(4) = dnorm(log_m_devs,0.0,dMiscCont(6));     
  }

  // | OBJECTIVE FUNCTION THAT ADMB WILL MINIMIZE
  f = sum(nll) + sum(penll) + sum(calcPriors()) + fpen;

  if(DEBUG_FLAG){
    COUT(f);
    COUT(nll);
    COUT(penll);
    COUT(calcPriors());
    if(fpen > 0 ){COUT(fpen);}
  }  


  //FUNCTION dvariable dnorm(const dvar_vector residual, const dvector std)
  //RETURN_ARRAYS_INCREMENT();
  //double pi=3.141593;
  //int n=size_count(residual);
  //dvector var=elem_prod(std,std);
  //dvar_vector SS=elem_prod(residual,residual);
  //RETURN_ARRAYS_DECREMENT();
  //return(0.5*n*log(2.*pi)+sum(log(std))+sum(elem_div(SS,2.*var)));

FUNCTION dvar_vector calcPriors()
  // |---------------------------------------------------------------------------------|
  // | PRIORS FOR LEADING PARAMETERS p(theta)
  // |---------------------------------------------------------------------------------|
  // | - theta_prior is a switch to determine which statistical distribution to use.
  // |
  dvariable ptmp; 
  dvar_vector priors(1,n_theta);
  priors.initialize();

  for(int i=1;i<=n_theta;i++)
  {
    ptmp = 0;
    // for(j=1;j<=ipar_vector(i);j++)
    // {
      if( active(theta(i)) )
      {
        switch(theta_iprior(i))
        {
        case 1:   //normal
          ptmp += dnorm(theta(i),theta_p1(i),theta_p2(i));
          break;
          
        case 2:   //lognormal CHANGED RF found an error in dlnorm prior. rev 116
          ptmp += dlnorm(theta(i),theta_p1(i),theta_p2(i));
          break;
          
        case 3:   //beta distribution (0-1 scale)
          double lb,ub;
          lb=theta_lb(i);
          ub=theta_ub(i);
          ptmp += dbeta((theta(i)-lb)/(ub-lb),theta_p1(i),theta_p2(i));
          break;
          
        case 4:   //gamma distribution
          ptmp += dgamma(theta(i),theta_p1(i),theta_p2(i));
          break;
          
        default:  //uniform density
          if(theta_lb(i) > theta_ub(i)){
            ptmp += log(1./(theta_lb(i)-theta_ub(i)));
          } else {
            ptmp += log(1./(theta_ub(i)-theta_lb(i)));
          }
          break;
        }
      }
    // }
    priors(i) = ptmp; 
  }
  return(priors);


GLOBALS_SECTION
  #include <admodel.h>
  #include <string.h>
  #include <time.h>

  #undef EOF
  #define EOF 999

  #undef TINY
  #define TINY  1.0e-10

  #undef REPORT
  #define REPORT(object) report << #object "\n" << setw(8) \
  << setprecision(4) << setfixed() << object << endl;

  #undef COUT
  #define COUT(object) cout << #object "\n" << setw(6) \
  << setprecision(3) << setfixed() << object << endl;

  template<typename T>
  dvar_vector plogis(const dvector x, T location, T scale)
  {
    return(1.0 / (1.0 + mfexp(-(x-location)/scale)));
  }

  dvar_vector plogis(const dvector x, dvariable location, dvariable scale)
  {
    return(1.0 / (1.0 + mfexp(-(x-location)/scale)));
  }

  dvar_vector plogis(const dvector x, prevariable location, prevariable scale)
  {
    return(1.0 / (1.0 + mfexp(-(x-location)/scale)));
  }

  template<typename T>
  dvar_vector plogis95(const dvector x, T a50, T a95)
  {
    return(1.0 / (1.0 + mfexp(-log(19)*(x-a50)/(a95-a50))));
  }


  dvariable dmvlogistic(const dmatrix o, const dvar_matrix& p,dvar_matrix& nu, double& tau2,const double minp)
  { //returns the negative loglikelihood using the MLE for the variance
  /*
    This is a modified version of the dmvlogistic negative log likelihood
    where proportions at age less than minp are pooled into the consecutive 
    age-classes.  See last paragraph in Appendix A of Richards, Schnute and
    Olsen 1997. 
    
    NB minp must be greater than 0, otherwise this algorithm returns an
    error if one of the observed proportions is zero.
    
    -1) first count the number of observations for each year > minp
    -2) normalized observed and predicted age-proportions
    -3) loop over ages, and check if observed proportion is < minp
        -if yes, then add observed proprtion to bin k
        -if no then add observed proportion to bin k and increment
         bin k if k is currently less than the number of bins.
    -4) do the same grouping for the predicted proportions.
    -5) use ivector iiage to correctly assign residuals into nu
    
    FEB 8, 2011.  Fixed a bug in the variance calculation & 
    likelihood scaling that was discovered at the 2011 Hake
    assessment STAR panel in Seattle.
  */
  RETURN_ARRAYS_INCREMENT();
  int i,j,k,n;
  int age_counts=0;
  int a = o.colmin();
  int A=o.colmax();
  double tiny=0.001/(A-a+1);
  int t=o.rowmin();
  int T=o.rowmax();
  dvariable tau_hat2;   //mle of the variance
  //dvar_matrix nu(t,T,a,A);
  nu.initialize();
  

  for(i=t; i<=T; i++)
  { 
    n=0;
    dvector oo = o(i)/sum(o(i));
    dvar_vector pp = p(i)/sum(p(i));
    
    //count # of observations greater than minp (2% is a reasonable number)
    for(j=a;j<=A;j++)
      if(oo(j) > minp)n++;
    
    ivector iiage(1,n);
    dvector o1(1,n); o1.initialize();
    dvar_vector p1(1,n); p1.initialize();
    k=1;
    for(j=a;j<=A;j++)
    {
      if(oo(j)<=minp)
      {
        o1(k)+=oo(j);
        p1(k)+=pp(j);
      }
      else
      {
        o1(k)+=oo(j);
        p1(k)+=pp(j);
        if(k<=n)iiage(k)=j;   //ivector for the grouped residuals
        if(k<n) k++;
      }
    }
    
    //assign residuals to nu based on iiage index
    dvar_vector t1 = log(o1)-log(p1) - mean(log(o1)-log(p1));
    

    for(j=1;j<=n;j++)
      nu(i)(iiage(j))=t1(j);
    
    age_counts += n-1;
  }
  //Depricated  Wrong Variance & likelihood calculation.
  //tau_hat2 = 1./(age_counts*(T-t+1))*norm2(nu);
  //dvariable nloglike =(age_counts*(T-t+1))*log(tau_hat2);
  
  //Feb 8, 2011  Fixed variance & likelihood
  tau_hat2 = 1./(age_counts)*norm2(nu);
  dvariable nloglike =(age_counts)*log(tau_hat2);
  tau2=value(tau_hat2); //mle of the variance 
  RETURN_ARRAYS_DECREMENT();
  return(nloglike);
  }
  double dicValue;
  double dicNoPar;
  void function_minimizer::mcmc_eval(void)
  {
    // |---------------------------------------------------------------------------|
    // | Added DIC calculation.  Martell, Jan 29, 2013                             |
    // |---------------------------------------------------------------------------|
    // | DIC = pd + dbar
    // | pd  = dbar - dtheta  (Effective number of parameters)
    // | dbar   = expectation of the likelihood function (average f)
    // | dtheta = expectation of the parameter sample (average y) 

    gradient_structure::set_NO_DERIVATIVES();
    initial_params::current_phase=initial_params::max_number_phases;
    uistream * pifs_psave = NULL;

  #if defined(USE_LAPLACE)
  #endif

  #if defined(USE_LAPLACE)
      initial_params::set_active_random_effects();
      int nvar1=initial_params::nvarcalc(); 
  #else
    int nvar1=initial_params::nvarcalc(); // get the number of active parameters
  #endif
    int nvar;
    
    pifs_psave= new
      uistream((char*)(ad_comm::adprogram_name + adstring(".psv")));
    if (!pifs_psave || !(*pifs_psave))
    {
      cerr << "Error opening file "
              << (char*)(ad_comm::adprogram_name + adstring(".psv"))
         << endl;
      if (pifs_psave)
      {
        delete pifs_psave;
        pifs_psave=NULL;
        return;
      }
    }
    else
    {     
      (*pifs_psave) >> nvar;
      if (nvar!=nvar1)
      {
        cout << "Incorrect value for nvar in file "
             << "should be " << nvar1 << " but read " << nvar << endl;
        if (pifs_psave)
        {
          delete pifs_psave;
          pifs_psave=NULL;
        }
        return;
      }
    }
  
    int nsamp = 0;
    double sumll = 0;
    independent_variables y(1,nvar);
    independent_variables sumy(1,nvar);

    do
    {
      if (pifs_psave->eof())
      {
        break;
      }
      else
      {
        (*pifs_psave) >> y;
        sumy = sumy + y;
        if (pifs_psave->eof())
        {
          double dbar = sumll/nsamp;
          int ii=1;
          y = sumy/nsamp;
          initial_params::restore_all_values(y,ii);
          initial_params::xinit(y);   
          double dtheta = 2.0 * get_monte_carlo_value(nvar,y);
          double pd     = dbar - dtheta;
          double dic    = pd + dbar;
          dicValue      = dic;
          dicNoPar      = pd;

          cout<<"Number of posterior samples    = "<<nsamp    <<endl;
          cout<<"Expectation of log-likelihood  = "<<dbar     <<endl;
          cout<<"Expectation of theta           = "<<dtheta   <<endl;
          cout<<"Number of estimated parameters = "<<nvar1    <<endl;
          cout<<"Effective number of parameters = "<<dicNoPar <<endl;
          cout<<"DIC                            = "<<dicValue <<endl;
          break;
        }
        int ii=1;
        initial_params::restore_all_values(y,ii);
        initial_params::xinit(y);   
        double ll = 2.0 * get_monte_carlo_value(nvar,y);
        sumll    += ll;
        nsamp++;
        // cout<<sumy(1,3)/nsamp<<" "<<get_monte_carlo_value(nvar,y)<<endl;
      }
    }
    while(1);
    if (pifs_psave)
    {
      delete pifs_psave;
      pifs_psave=NULL;
    }
    return;
  }

REPORT_SECTION
// Write out Raw Data (Useful for simulation studies)
// Note that I use a Macro called REPORT here to ensure a standard format.
  REPORT(dat_syr);
  REPORT(dat_nyr);
  REPORT(mod_syr);
  REPORT(mod_nyr);
  REPORT(sage);
  REPORT(nage);
  REPORT(nFecBlocks);
  REPORT(nFecBlockYears);
  REPORT(fec_slope);
  REPORT(fec_inter);
  REPORT(data_ct_raw);
  REPORT(data_sp_waa);
  REPORT(data_cm_waa);
  REPORT(data_cm_comp);
  REPORT(data_sp_comp);
  REPORT(data_egg_dep);
  REPORT(data_mileday);
  REPORT(EOF)
// END of Replicated Data File. (run model with -noest to check data)
  
// Negative loglikelihoods
  REPORT(nll);

// Vectors of years.
  ivector year(mod_syr,mod_nyr);
  year.fill_seqadd(mod_syr,1);
  
  ivector years(mod_syr,mod_nyr+1);
  years.fill_seqadd(mod_syr,1);

  ivector rec_years(rec_syr,mod_nyr+1);
  rec_years.fill_seqadd(rec_syr,1);

  ivector iage = ivector(age);

  REPORT(iage);
  REPORT(year);
  REPORT(years);
  REPORT(rec_years);
  REPORT(data_catch);
  REPORT(pred_catch);

// SSB, recruits, spawners,
  REPORT(ssb);
  REPORT(spawners);
  REPORT(recruits);
  REPORT(pred_egg_dep);

// Numbers-at-age of various flavors.
  REPORT(Nij);
  REPORT(Oij);
  REPORT(Pij);
  REPORT(Cij);

// Selectivity and vulnerable proportion-at-age.
  REPORT(Sij);
  REPORT(Qij);

// Natural mortality
  REPORT(Mij);

// Fishing mortality
  dvector ft = value(exp(log_ft_pars));
  REPORT(ft);
  REPORT(Fij);


// Residuals
  REPORT(resd_egg_dep);
  REPORT(resd_mileday);
  REPORT(resd_rec);
  REPORT(resd_catch);
  REPORT(resd_sp_comp);
  REPORT(resd_cm_comp);
  
// Forecast output
  REPORT(fore_sb);
  REPORT(fore_vb);
  REPORT(ghl);

// Initial parameter values from simulation studies.
  REPORT(theta_ival);
  REPORT(theta);

// Simulation spawning stock biomass
  if(b_simulation_flag) {
    REPORT(sim_spawners);
    REPORT(sim_recruits);
  }


FINAL_SECTION
  //system("cp ham.rep run1.rep");


























