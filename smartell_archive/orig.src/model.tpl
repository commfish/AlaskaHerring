// ----------------------------------------------------------------------------- //
//               Age-structured model for Alaska herring stocks                  //
//                                                                               //
//                               VERSION 0.1                                     //
//                                Jan  2015                                      //
//                                                                               //
//                                 AUTHORS                                       //
//                              Sherri Dressel                                   //                                 
//                        sherri.dressel@alaska.gov                              //
//                               Sara Miller                                     //
//                          sara.miller@alaska.gov                               //
//                               Kray Van Kirk                                   //
//                          kray.vankirk@alaska.gov                              //
//                                                                               //
//                   Built on code developed by Peter Hulson                     //
//                            pete.hulson@noaa.gov                               //
//                                                                               //
//                           Layout and references                               //
//                              Steven Martell                                   //
//                          martell.steve@gmail.com                              //
//                                                                               //
// CONVENTIONS: Formatting conventions are based on                              //
//              The Elements of C++ Style (Misfeldt et al. 2004)                 //
//                                                                               //
//                                                                               //
// NAMING CONVENTIONS:                                                           //
//             Macros       -> UPPERCASE                                         //
//             Constants    -> UpperCamelCase                                    //
//             Functions    -> lowerCamelCase                                    //
//             Variables    -> lowercase                                         //
//                                                                               //
//                                                                               //
//                                                                               //
// ----------------------------------------------------------------------------- //
//-- CHANGE LOG:                                                               --//
//--  Jan 2015 - revision of legacy code::                                     --//
//--               :variable naming conventions                                --//
//--               :intra-annual calendar                                      --//
//--               :standardization of units across stocks                     --//
//--               :modification for potential code distribution               --//
//--  May 2015 - added code and scripts to a git-hub repo.                     --//
//--                                                                           --//
// ----------------------------------------------------------------------------- //


DATA_SECTION

  // |--------------------------------------------------------------------------|
  // |MCMC OUTPUT FILE
  // |--------------------------------------------------------------------------|

     !!CLASS ofstream evalout("evalout.prj");

  // |--------------------------------------------------------------------------|
  // | STRINGS FOR INPUT FILES                                                  |
  // |--------------------------------------------------------------------------|
  // |- The files below are listed in the 'model.dat' file;
  // |   nothing else should be in the 'model.dat' file
  // |- These files should be named by the stock they are modeling;
  // |   example: "sitka.dat", "seymour.dat"
  // |- DO NOT use two word names such as "seymour cove.dat"
  // |
  // | DataFile               : data to condition the assessment model    
  // | ControlFile            : controls for years, phases, and block options 
  // | Graphics               : vectors for some stock assessment graphics

  init_adstring DataFile;      
  init_adstring ControlFile;
  init_adstring Graphics;    

  // | BaseFileName           : file prefix used for all  model output
  // | ReportFileName         : file name to which report file is printed

  !! BaseFileName = stripExtension(DataFile);  
  !! ReportFileName = BaseFileName + adstring(".rep");

  !! cout<<"You are modeling the "<<BaseFileName<<" stock of herring"<<endl;
  !! cout<<""<<endl;
  !! time(&start);

	
  // |--------------------------------------------------------------------------|
  // | MODEL DATA FROM CONTROL FILE                                             |
  // |--------------------------------------------------------------------------|
  // | This calls the control file
  // | This file controls model years, estimation phases, blocks for
  // | selectivity, maturity, and mortality (and fecundity if needed)
  // | and anything else that is NOT observed data (catch, weight, eggs, etc.)

  !! ad_comm::change_datafile_name(ControlFile);

  // |--------------------------------------------------------------------------|
  // | MODEL STRUCTURAL DIMENSIONS                                              |
  // |--------------------------------------------------------------------------|
  // | 
  // |-nages      -> number of ages
  // |-dat_styr   -> data start year
  // |-dat_endyr  -> data end year
  // |-mod_styr   -> model start year
  // |-mod_endyr  -> model end year
  // |-dyrs       -> data year index
  // |-myrs       -> model year index
  // |-md_offset  -> offset data to model
  // |-Year       -> year sequence for graphics (model + 1)
  // | 
  // | year         i
  // | age          j
  // | 

  init_int nages
  init_int dat_styr
  init_int dat_endyr
  init_int mod_styr
  init_int mod_endyr
 

  int dyrs
  int myrs
  int md_offset
  vector Year(mod_styr,mod_endyr+1);       


  int i;
  int j;
  
 LOCAL_CALCS

  dyrs=dat_endyr-dat_styr+1;
  myrs=mod_endyr-mod_styr+1;
  md_offset=mod_styr-dat_styr;
  Year.fill_seqadd(mod_styr,1);   

 END_CALCS 


  // |--------------------------------------------------------------------------|
  // | FORECASTING                                                              |
  // |--------------------------------------------------------------------------|
  // | 
  // | These are used to forecast net year's spawning biomass and whether
  // | that biomass exceeds the minimum threshold below which no commercial
  // | fishery is implemented 
  // |  
  // | -Thresh     ->  Threshold (in short tons) below which no fishery occurs
           
     init_number Thresh;
     init_number Thresh_denom;
     init_vector Threshold(1,myrs);
     init_vector fw_a_a(1,nages);
   !!cout <<myrs<<endl; 
  // |--------------------------------------------------------------------------|
  // | ARRAY INDEXING                                                           |
  // |--------------------------------------------------------------------------|
  // | 
  // | These are the points at which the model allows natural mortality (M),
  // | maturity-at-age (mat), gear selectivity-at-age (gs) to change based on
  // | climate changes indicated by shifts in the Pacific Decadal Oscillation. 
  // | Fecundity, at this point, remains stable with one estimate 
  // |  
  // | -F_Bk       ->  number of fecundity blocks (1 split = 2 blocks)
  // | -F_Bk_Yrs   ->  specific years in which fecundity changes
  // | -F_slope    ->  Slope of fecundity-at-age regression
  // | -F_inter    ->  Intercept of fecundity
  // | -M_Bk       ->  number of mortality blocks (1 split = 2 blocks)
  // | -M_Bk_Yrs   ->  specific years in which mortality M changes
  // | -mat_Bk     ->  number of maturity blocks (1 split = 2 blocks)
  // | -mat_Bk_Yrs ->  specific years in which maturity-at-age changes
  // | -gs_Bk      ->  number of gear selectivity blocks (1 split = 2 blocks)
  // | -gs_Bk_Yrs  ->  specific years in which gear-selectivity-at-age changes
 
  init_number F_Bk
  init_vector F_Bk_Yrs(1,F_Bk+1)
  init_vector F_slope(1,F_Bk)
  init_vector F_inter(1,F_Bk)
  vector F_Bk_Idx(1,F_Bk+1)

  init_number S_Bk
  init_vector S_Bk_Yrs(1,S_Bk+1)
  vector S_Bk_Idx(1,S_Bk+1)

  init_number mat_Bk
  init_vector mat_Bk_Yrs(1,mat_Bk+1)
  vector mat_Bk_Idx(1,mat_Bk+1)

  init_number gs_Bk
  init_vector gs_Bk_Yrs(1,gs_Bk+1)
  vector gs_Bk_Idx(1,gs_Bk+1)


  // |--------------------------------------------------------------------------|
  // | ESTIMATION PHASES                                                        |
  // |--------------------------------------------------------------------------|
  // |
  // | These govern the phases in which given parameters are estimated.
  // | While these should likely not need adjusting, if you encounter problems,
  // | you might explore some changes. A parameter whose estimation phase
  // | is > 1 remains at its starting value for all phases < estimation phase.
  // | A parameter with a negative phase is not estimated.
  // |
  // | Phase 1: -ph_Int    -> initial population
  // |          -ph_S      -> natural mortality
  // |          -ph_mat_a  -> maturity inflection
  // |          -ph_gs_a   -> gear selectivity inflection
  // |          -ph_gs_b   -> gear selectivity slope
  // | Phase 2: -ph_mat_b  -> maturity slope
  // | Phase 3: -ph_Rec    -> recruitment (age 3)
  // |          -ph_Ric    -> Ricker function

  init_number ph_Int
  init_number ph_S
  init_number ph_mat_a
  init_number ph_gs_a
  init_number ph_gs_b
  init_number ph_mat_b
  init_number ph_Rec
  init_number ph_Ric
  init_number ph_md		


  // |--------------------------------------------------------------------------|
  // | OBJECTIVE FUNCTION WEIGHTS                                               |
  // |--------------------------------------------------------------------------|
  // |

  init_number lC                          //Catch age composition
  init_number lS                          //Spawning age composition
  init_number lR                          //Ricker spawner-recruit//check every few years to make sure weight has no influence on model (see Ricker comparison 2011 for sitka)
                                          // run model with 0 weight and then compare to different weights, want the smallest weight above zero with no influence on recruits
  init_vector lE(1,dyrs)                  //Egg deposition
  init_vector lM(1,dyrs)         	  //Mile-days

  // |--------------------------------------------------------------------------|
  // | END OF FILE MARKER                                                       |
  // |--------------------------------------------------------------------------|
  init_number eof1

 LOCAL_CALCS

    if(eof1==42) cout << BaseFileName<<".ctl has been read correctly!"<<endl;
    else 
    {    
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|      Red alert! Captain to bridge! The .ctl file is compromised!     |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof1<<", but the file *should* end with 42      |"<<endl;
         cout <<"| Please check the .ctl file for errors and make sure the above calls  |"<<endl;
         cout <<"|              are matched exactly by the file's contents              |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
    exit(1); 
    }



   // | The lines below populate the indices for fecundity, mortality,
   // | selectivity and maturity

   for (int i=1;i<=F_Bk+1;i++)
     {
       F_Bk_Idx(i)=F_Bk_Yrs(i)-mod_styr+1;
     } 

   for (int i=1;i<=S_Bk+1;i++)
     {
       S_Bk_Idx(i)=S_Bk_Yrs(i)-mod_styr+1;
     }

   for (int i=1;i<=mat_Bk+1;i++)
     {
       mat_Bk_Idx(i)=mat_Bk_Yrs(i)-mod_styr+1;
     }

   for (int i=1;i<=gs_Bk+1;i++)
     {
       gs_Bk_Idx(i)=gs_Bk_Yrs(i)-mod_styr+1;
     }

 END_CALCS

  // |--------------------------------------------------------------------------|
  // | MODEL DATA FROM DATA FILE                                                |
  // |--------------------------------------------------------------------------|
  // | This calls the data file

  !! ad_comm::change_datafile_name(DataFile);

  // |--------------------------------------------------------------------------|
  // | MODEL DATA                                                               |
  // |--------------------------------------------------------------------------|
  // | 
  // |-tcb            -> total annual commercial catch biomass (short tons)
  // |-obs_sp_waa     -> spawning (cast net) weight-at-age
  // |-obs_c_waa      -> commercial catch weight-at-age
  // |-obs_c_comp     -> commercial catch age composition
  // |-obs_sp_comp    -> spawner (cast net) age composition
  // |-obs_egg        -> egg deposition
  // |-mile_days      -> mile days of milt (not always used)
  // | 

  init_vector tcb(1,dyrs)                
  init_matrix obs_sp_waa(1,dyrs,1,nages) 
  init_matrix obs_c_waa(1,dyrs,1,nages)
  init_matrix obs_c_comp(1,dyrs,1,nages)
  init_matrix obs_sp_comp(1,dyrs,1,nages)
  init_vector tot_obs_egg(1,dyrs)
  init_vector mile_days(1,dyrs)
  vector tcbm(1,dyrs)

  // |--------------------------------------------------------------------------|
  // | END OF FILE MARKER                                                       |
  // |--------------------------------------------------------------------------|

     init_number eof2

 LOCAL_CALCS

    if(eof2==42) cout << BaseFileName<<".dat has been read correctly!"<<endl;
    else 
    {       
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|   ** Red alert! Captain to bridge! The .dat file is compromised! **  |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof2<<", but the file *should* end with 42      |"<<endl;
         cout <<"| Please check the .dat file for errors and make sure the above calls  |"<<endl;
         cout <<"|              are matched exactly by the file's contents              |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
    exit(1); 
    }

    // |This line converts annual catch biomass from short tons to metric tons
    for (int i=1;i<=dyrs;i++)
    {
       tcbm(i)=0.90718*tcb(i);
    } 

    // |This line extracts the last line of spawning weight-at-age for forecasting
    //fw_a_a = obs_sp_waa(1,dyrs);
   
 END_CALCS


  // |--------------------------------------------------------------------------|
  // | HISTORICAL ESTIMATES FROM GRAPHICS FILE                                  |
  // |--------------------------------------------------------------------------|
  // | This calls the graphics file that holds historical model estimates
  // | for the R graphics template
  // | THis file will be updated each year (by you!)

  !! ad_comm::change_datafile_name(Graphics);

  // |--------------------------------------------------------------------------|
  // | GRAPHICS CONTENTS                                                        |
  // |--------------------------------------------------------------------------|
  // |-model years
  // |-yminusthree       -> ASA model (pre-fishery biomass) three years ago
  // |-yminustwo         -> ASA model (pre-fishery biomass) two years ago
  // |-yminusone         -> ASA model (pre-fishery biomass) one year ago
  // |-yminusthreeFOR    -> pre-biomass forecast (short tons) three years ago
  // |-yminustwoFOR      -> pre-biomass forecast (short tons) two years ago
  // |-yminusoneFOR      -> pre-biomass forecast (short tons) one year ago
  // | 
  init_vector mod_yrs(1,myrs)
  init_vector yminusthree(1,myrs)
  init_vector yminustwo(1,myrs)
  init_vector yminusone(1,myrs)    
  init_number yminusthreeFOR
  init_number yminustwoFOR
  init_number yminusoneFOR
  init_number eof3

 LOCAL_CALCS

    if(eof3==42) cout << BaseFileName<<"_graphics.ctl has been read correctly!"<<endl;
    else 
    {       
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|   Red alert! Captain to bridge! The graphics file is compromised!    |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof3<<", but the file *should* end with 42      |"<<endl;
         cout <<"|  Please check the graphics file for errors and make sure the above   |"<<endl;
         cout <<"|          calls are matched exactly by the file's contents            |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
    exit(1); 
    }

   
 END_CALCS

PARAMETER_SECTION

  // |---------------------------------------------------------------------------------|
  // | INITIAL POPULATION PARAMETERS
  // |---------------------------------------------------------------------------------|
  // | 
  // |- initial age 3 abundance
  // |- initial population abundance (ages 4 - 8+)
	
  init_bounded_vector init_age_3(1,myrs,0,2500,ph_Rec)
   
  init_bounded_vector init_pop(1,nages-1,0,500,ph_Int) 


  // |---------------------------------------------------------------------------------|
  // | SPAWNER-RECRUIT FUNCTION PARAMETERS
  // |---------------------------------------------------------------------------------|
  // | 
  // |- log_alpha: this is a natural log-scale parameter 
  // |- log_beta:  this is a natural log-scale parameter 

  init_bounded_number log_alpha(-10,0,ph_Ric);
  init_number log_beta(ph_Ric);                                  
  number alpha;
  number beta;
 // init_bounded_number md_c(200,600,ph_md)	    //Mile-days coefficient
 

  // |---------------------------------------------------------------------------------|
  // | MORTALITY / SURVIVAL PARAMETERS
  // |---------------------------------------------------------------------------------|
  // | 
  // |- natural mortality
  // |- survival matrix
  // |- survival vector
  // |- mean survival (for forecasting)

  init_bounded_vector M(1,S_Bk,.0001,1,ph_S);                   
  matrix Sur(1,myrs,1,nages);    
  vector S(1,S_Bk);               
  number S_for;   
              

  // |---------------------------------------------------------------------------------|
  // | MATURITY & SELECTIVITY  PARAMETERS
  // |---------------------------------------------------------------------------------|
  // | 
  // |- age at 50% maturity
  // |- maturity-at-age slope
  // |- age at 50% selectivity
  // |- selectivity-at-age slope
  // |- mile-days coefficient

  init_bounded_vector mat_a(1,mat_Bk,1,10,ph_mat_a);    
  init_bounded_vector mat_b(1,mat_Bk,0,5,ph_mat_b);                    
  init_bounded_vector gs_a(1,gs_Bk,2,10,ph_gs_a);
  init_bounded_vector gs_b(1,gs_Bk,0,5,ph_gs_b);
  init_bounded_number md_c(200,600,ph_Int);
  matrix GS(1,myrs,1,nages);
  matrix GS_Sc(1,myrs,1,nages);
  matrix Mat(1,myrs,1,nages);
  vector int1(1,myrs);              //Maximum of gear selectivity, used to scale
  matrix int2(1,myrs,1,nages);      //Proportion of pop rec to gear*Weight_at_age
  vector int3(1,myrs);              //sum over age of int2
  vector mat_for(1,nages);


  // |---------------------------------------------------------------------------------|
  // | ESTIMATED AND DERIVED POPULATION MATRICES
  // |---------------------------------------------------------------------------------|
  // | 
  // |- naa            total population [mature + immature]               [millions]       
  // |- sel_naa        N selected by gear                           
  // |- sel_naa_prop   proportion of N selected by gear            
  // |- est_c_naa      catch composition-at-age                                  
  // |- est_sp_naa     spawning numbers-at-age                            [millions]            
  // |- est_sp_comp    spawners numbers-at-age                            [proportion]            
  // |- est_sp_baa     spawning biomass-at-age                            [metric tons]
  // |- est_mat_naa    mature abundance-at-age                            [millions]                       
  // |- est_mat_baa    mature biomass-at-age                              [metric tons]          
  // |- post_naa       total numbers-at-age [mature+ immature] - catch    [millions]     
  // |- est_egg_naa    egg production-at-age                              [trillions]     


  matrix naa(1,myrs,1,nages);           
  matrix sel_naa(1,myrs,1,nages);
  matrix sel_naa_prop(1,myrs,1,nages);
  matrix est_c_naa(1,myrs,1,nages);
  matrix est_sp_naa(1,myrs,1,nages);
  matrix est_sp_comp(1,myrs,1,nages);
  matrix est_sp_baa(1,myrs,1,nages);
  matrix est_mat_naa(1,myrs,1,nages);
  matrix est_mat_baa(1,myrs,1,nages);
  matrix post_naa(1,myrs,1,nages);
  matrix est_egg_naa(1,myrs,1,nages);
     
  
  // |---------------------------------------------------------------------------------|
  // | ESTIMATED AND DERIVED POPULATION VECTORS
  // |---------------------------------------------------------------------------------|
  // | 
  // |- tot_sel_N    total N selected by gear                    [millions]                      
  // |- tot_sp_B     total spawning biomass                      [metric tons]
  // |- tot_mat_N    total mature abundance                      [millions]
  // |- tot_sp_N     total spawning abundance                    [millions]                    
  // |- tot_mat_B    total mature biomass                        [metric tons]   
  // |- tot_post_N   total population [mature+ immature] - catch [millions]     
  // |- tot_est_egg  total egg deposition                        [trillions]                          
  // |- N            total abundance (mature + immature)         [millions]
  // |
  // |- Ricker spawner-recruit
  // |- Mile days of milt
 
  vector          tot_sel_N(1,myrs);
  vector          tot_sp_B(1,myrs);               
  vector          tot_mat_N(1,myrs);             
  sdreport_vector tot_sp_N(1,myrs);     
  sdreport_vector tot_mat_B(1,myrs); 
  sdreport_vector tot_post_N(1,myrs); 
  sdreport_vector tot_est_egg(1,myrs);
  sdreport_vector N(1,myrs);
  vector          SR(1,myrs-3);                 //Ricker Spawner-Recruit
  vector M_D(1,myrs);	 	      //Estimated Mile-days of milt
 

  // |---------------------------------------------------------------------------------|
  // | FORECAST QUANTITIES
  // |---------------------------------------------------------------------------------|
  // | 
  // |- for_naa         forecast numbers-at-age  [mature + immature]       [millions]
  // |- for_mat_naa     forecast mature numbers-at-age                     [millions]
  // |- for_mat_baa     forecast mature biomass-at-age                     [metric tons]  
  // |- for_mat_prop    forecast mature proportion-at-age [by number]      [proportion]
  // |- for_mat_b_prop  forecast mature proportion-at-age [by biomass]      [proportion]
  // |- for_mat_B       total forecast mature biomass                      [metric tons]
  // |- for_mat_B_st    total forecast mature biomass            [short dweeby US tons]  
  // |- for_tot_mat_N   total mature numbers-at-age                   
  // |- HR              harvest rate                                 
  // |- HR_p            harvest rate sliding proportion               
  // |- GHL             general harvest limit                         

  vector for_naa(1,nages)
  vector for_mat_naa(1,nages)        
  vector for_mat_baa(1,nages)      
  vector for_mat_prop(1,nages)
  vector for_mat_b_prop(1,nages)
  number for_mat_B                  
  number for_mat_B_st		  
  number for_tot_mat_N
  number HR		         
  number HR_p                  
  number GHL ;                 

  // |---------------------------------------------------------------------------------|
  // | GRAPHICAL CONSTRUCTS
  // |---------------------------------------------------------------------------------|
  // | 
  // |- These matrices are read in R for standardized graphical analyses and figures
  // |- FIGDATA
  // |- FIGDATAAGE
  // |- FIGDATA2

   matrix FIGDATA(1,myrs,1,42);
   matrix FIGDATAAGE(1,nages,1,3);


  // |---------------------------------------------------------------------------------|
  // | OBJECTIVE FUNCTION COMPONENTS
  // |---------------------------------------------------------------------------------|
  // | 
  // |* Residuals *
  // |- catch age composition      
  // |- spawner age composition
  // |- egg deposition
  // |- spawner-recruit function
  // |
  // |* Sums of squares *
  // |- catch age composition
  // |- spawner age composition
  // |- weighted egg deposition vector
  // |- egg deposition
  // |- spawner-recruit
  // |- Mile-days of milt vector
  // |- MD SSM
  // |
  // |- ADMB minimization target f

  matrix res_c_comp(1,myrs,1,nages);
  matrix res_sp_comp(1,myrs,1,nages);
  vector res_tot_egg(1,myrs);
  vector res_SR(1,myrs-3);

  number SSQC;
  number SSQSp;
  vector wSSQE(1,myrs);
  number WSSQE;
  number SSQR;
  vector M_DR(1,myrs);
  vector wSSQM(1,myrs);        
  number WSSQM;	          

  objective_function_value f;




  // |---------------------------------------------------------------------------------|
  // | AICc CONSTRUCTS
  // |---------------------------------------------------------------------------------|
  // | 
  // |* Residuals *
  // |- number of observations: catch age composition      
  // |- number of observations: spawner age composition
  // |- number of observations: egg deposition
  // |
  // |- vector of sample size
  // |- weighting vector
  // |- log-likelihood vector
  // |
  // |* Sums of squares *
  // |- sum catch age observations
  // |- sum spawner age observations
  // |- sum egg observations
  // |- sum all observations
  // |
  // |- AIC weighting terms
  // |- Total L(like)
  // |- AIC
  // |- AICc
  // |- number of parameters
  // |

  matrix n_obs_c_comp(1,myrs,1,nages);
  matrix n_obs_sp_comp(1,myrs,1,nages); 
  vector n_tot_obs_egg(1,myrs)

  vector n_d(1,myrs+3)
  vector w_d(1,myrs+3)
  vector lnL_d(1,myrs+3)

  number n_C
  number n_S
  number n_R
  number n

  number sig_1
  number lnL
  number AIC
  number AICc
  number p

  vector n_My(1,myrs)			//Index vector for number of observations
  number n_M				//Number of mile-days observations
  vector sig_d(1,myrs+3)		//Sigma vector

  // |---------------------------------------------------------------------------------|
  // | BOOTSTRAP ELEMENTS - not currently implemented 1/13/2015 -kvk
  // |---------------------------------------------------------------------------------|
  // |
  // |- Define bootstrap data matrix
  // |- Index
  // |- Index
  // |- Index
  // |- Index
  // |- Rejection criteria

  //matrix BSDATA(1,myrs,1,26)
  //sdreport_number Data_start_year
  //sdreport_number Data_end_year   
  //sdreport_number Model_start_year
  //sdreport_number Model_end_year
  //sdreport_number obj


  // |---------------------------------------------------------------------------------|
  // | A few SD report calls
  // |---------------------------------------------------------------------------------|
  // |
  // |- Spawning biomass
  // |- Forecase spawning biomass
  // |- GHL
  // |- AIC

  sdreport_vector SpawnBio(1,myrs)
  sdreport_number SpawnBioFor
  sdreport_number GHLsd
  sdreport_number AICcsd

 LOCAL_CALCS

  //Convert observed spawner-at-age to proportions
  //for (int i=1; i<=dyrs; i++)
  //{
   //obs_sp_comp(i)/=sum(obs_sp_comp(i)+0.00001);
  //}

 END_CALCS

PRELIMINARY_CALCS_SECTION

  // |---------------------------------------------------------------------------------|
  // | Needed? .PIN file better option?
  // |---------------------------------------------------------------------------------|

  log_beta=-10;                       
  mat_b=1;


PROCEDURE_SECTION

  get_parameters();
  Time_Loop();
  get_residuals();
  evaluate_the_objective_function();
  if(last_phase())
    {
      get_forecast();
    } 


  if(sd_phase())
    {
      compute_AICc();
     get_FIGDATA();
     output_FIGDATA();
     get_FIGDATAAGE();
     output_FIGDATAAGE();
      get_sdreports();
      get_report();
    }


  if(mceval_phase())
    {
     evalout<<for_mat_B_st<<" "
          <<init_age_3<<" "
          <<init_pop<<" "<<endl;
    }


FUNCTION get_parameters

  /**
    * Defines matrices for survival, maturity, and selectivity
    * based on the values entered in the Control File for:
    * S_Bk, mat_BK, and gs_BK
  **/

  Sur.initialize();
  Mat.initialize();
  GS.initialize();
  GS_Sc.initialize();
  int1.initialize();
  S_for.initialize();
  mat_for.initialize();

  //Survival
  for (int t=1;t<=S_Bk;t++)
  {
    S(t)=exp(-M(t));
      for (int i=S_Bk_Idx(t);i<=S_Bk_Idx(t+1);i++)
      {
         for (int j=1;j<=nages;j++)
         {
           Sur(i,j)=exp(-M(t));
         }
      }
  }

  S_for=Sur(myrs,1);

  //Maturity
  for (int t=1;t<=mat_Bk;t++)
  {
    for (int i=mat_Bk_Idx(t);i<=mat_Bk_Idx(t+1);i++)
       {
         for (int j=1;j<=nages;j++)
          {
            Mat(i,j)=1/(1+exp(-1.0*mat_b(t)*((j+2)-mat_a(t))));
          }
       }
   }

   mat_for = Mat(myrs);

  //Gear selectivity-at-age (scaled)
  for (int t=1;t<=gs_Bk;t++)
    {
      for (int i=gs_Bk_Idx(t);i<=gs_Bk_Idx(t+1);i++)
        {
          for (int j=1;j<=nages;j++)
            {
              GS(i,j)=1/(1+exp(-1.0*gs_b(t)*((j+2)-gs_a(t))));
            }
        }
    }

  for (int i=1;i<=myrs;i++)
    {
      int1(i)=max(GS(i));   // This is not differentiab
        for (int j=1;j<=nages;j++)
        {
          GS_Sc(i,j)=GS(i,j)/int1(i);
        }
     }


FUNCTION Time_Loop

  naa.initialize();
  N.initialize();
  sel_naa.initialize();
  tot_sel_N.initialize();
  sel_naa_prop.initialize();
  est_c_naa.initialize();
  est_sp_naa.initialize();
  tot_sp_N.initialize();
  est_sp_comp.initialize();
  est_sp_baa.initialize();
  tot_sp_B.initialize();
  est_mat_baa.initialize();
  tot_mat_B.initialize();
  post_naa.initialize(); 
  tot_post_N.initialize();
  est_egg_naa.initialize();
  tot_est_egg.initialize();
  SR.initialize();
  int2.initialize();
  int3.initialize();
  M_D.initialize();


  /**
    *Pre-Fishery & Catch: This loop is the basis for the model. To perform
    *calculations we first set up the first year, get catch, then proceed to
    *calculate for remaining years
  **/

  //----------------------------------------------------------------------------
  // YEAR ONE
  //----------------------------------------------------------------------------
  
  for(int i=1;i<=1;i++)
   {
        naa(i,1)=init_age_3(i);
     for(int j=2;j<=nages;j++)
       {  //recruitment vector - year 1
         naa(i,j)=init_pop(j-1);              //initial population - year 1
       }
     for(int j=1;j<=nages;j++)
       {
         sel_naa(i,j)=naa(i,j)*GS_Sc(i,j);    //numbers-at-age vulnerable to gear (make sure this statement is separate from above)
       }
   
     tot_sel_N=rowsum(sel_naa);               //total numbers vulnerable to gear

   
    for(int j=1;j<=nages;j++)
     {
       sel_naa_prop(i,j)=sel_naa(i,j)/tot_sel_N(i);          //proportion-at-age: gear
       int2(i,j)=sel_naa_prop(i,j)*obs_c_waa(i+md_offset,j); //proporton * wt_g
     }

     int3=rowsum(int2);


    for(int j=1;j<=nages;j++)
      {
        est_c_naa(i,j)=tcbm(i+md_offset)*sel_naa_prop(i,j)/int3(i); //est. caa as % from tcbm
      }


    for(int j=1;j<=nages;j++)
      {
        post_naa(i,j)=naa(i,j)-est_c_naa(i,j); //numbers - catch
      }
  } 

  //----------------------------------------------------------------------------
  // END FIRST YEAR LOOP
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  // ALL OTHER YEARS
  //----------------------------------------------------------------------------


  for(int i=2;i<=myrs;i++)
    {
    
    for(int j=2;j<=nages;j++)
      {   
        naa(i,1)=init_age_3(i);                 //recruitment vector, years > 1     
        naa(i,j)=post_naa(i-1,j-1)*Sur(i,j-1);  //naa: (numbers-catch)*survival
      }


    for(int j=nages;j<=nages;j++)
      {
        naa(i,j)=post_naa(i-1,j-1)*Sur(i,j-1)+post_naa(i-1,j)*Sur(i,j); //+ class, naa
      }


    for(int j=1;j<=nages;j++)
      {
        sel_naa(i,j)=naa(i,j)*GS_Sc(i,j); //naa vulnerable to gear
      }
    
     tot_sel_N=rowsum(sel_naa);


    for(int j=1;j<=nages;j++)
      {
        sel_naa_prop(i,j)=sel_naa(i,j)/tot_sel_N(i);           //proportion-at-age: gear
        int2(i,j)=sel_naa_prop(i,j)*obs_c_waa(i+md_offset,j);  //proporton * wt_kg
      }

      int3=rowsum(int2);


    for(int j=1;j<=nages;j++)
      {
        est_c_naa(i,j)=tcbm(i+md_offset)*sel_naa_prop(i,j)/int3(i); //est. caa as % from tcb
      }


    for(int j=1;j<=nages;j++)
      {
        post_naa(i,j)=naa(i,j)-est_c_naa(i,j); //numbers-catch
      }
      
  }  //End other years Loop

    N=rowsum(naa);                  // total abundance (millions)
    tot_post_N = rowsum(post_naa);  // total abundance [numbers - catch] (millions)

  

  /**
    *spawning population quantities
    *Recall that "spawning population" is distinct from "mature population"
    *because spawning = mature - fishery
    *If there is NO fishery, spawning population = mature population
  **/

  
  //----------------------------------------------------------------------------
  // SPAWNERS & MATURE
  //----------------------------------------------------------------------------

  for (int i=1;i<=myrs;i++)
    {
      
      for(int j=1;j<=nages;j++)
        {
          est_sp_naa(i,j)=(Mat(i,j)*naa(i,j))-est_c_naa(i,j);          //spawning numbers-at-age (differs for stocks) For Craig should be Mat *(naa-est_c_na)a
        }}

    tot_sp_N=rowsum(est_sp_naa);                                       //total spawning numbers
  
  for (int i=1;i<=myrs;i++)
    {
      
      for(int j=1;j<=nages;j++)
        {
          est_sp_baa(i,j)=obs_sp_waa(i+md_offset,j)*est_sp_naa(i,j);   //spawning biomass-at-age
        }}

    tot_sp_B=rowsum(est_sp_baa);                                       //total spawning biomass


  for (int i=1;i<=myrs;i++)
    {
      
      for(int j=1;j<=nages;j++)
        {
          est_mat_naa(i,j)=(Mat(i,j)*naa(i,j));                        //mature-at-age

        }}

    tot_mat_N=rowsum(est_mat_naa);                                     //total mature numbers
      

  for (int i=1;i<=myrs;i++)
    {
      
      for(int j=1;j<=nages;j++)
        {
          est_mat_baa(i,j)=est_c_naa(i,j)*obs_c_waa(i+md_offset,j)+est_sp_naa(i,j)*obs_sp_waa(i+md_offset,j); //mature biomass at age
        }}

    tot_mat_B=rowsum(est_mat_baa);                                     //total mature biomass


  for (int i=1;i<=myrs;i++)
    {
   for(int j=1;j<=nages;j++)
        {
          est_sp_comp(i,j)=est_sp_naa(i,j)/tot_sp_N(i);                //spawning % at age
        }}


  //----------------------------------------------------------------------------
  // EGG DEPOSITION
  //----------------------------------------------------------------------------

  for (int t=1;t<=F_Bk;t++)
    {
      for (int i=F_Bk_Idx(t);i<=F_Bk_Idx(t+1);i++)
        {
          for (int j=1;j<=nages;j++)
            {
              est_egg_naa(i,j)=0.5*est_sp_naa(i,j)*(F_slope(t)*obs_sp_waa(i+md_offset,j)-F_inter(t))*0.000001;
            }
        //cout<<est_egg_naa(i)<<endl;
        }
    }
    //exit(1);
  tot_est_egg=rowsum(est_egg_naa);

  //----------------------------------------------------------------------------
  // SR FUNCTION
  //----------------------------------------------------------------------------

  alpha=exp(log_alpha);
  beta=exp(log_beta);

  for(int i=1;i<=myrs-3;i++)
    {
      SR(i)=alpha*tot_sp_B(i)*exp(-1.0*beta*tot_sp_B(i));
    }
 
  //Mile-days of milt
    for (int i=1;i<=myrs;i++){
    M_D(i)=0.5*tot_sp_B(i)/md_c;}


FUNCTION get_residuals

  res_c_comp.initialize();
  res_sp_comp.initialize();
  res_tot_egg.initialize();
  res_SR.initialize();
  M_DR.initialize();


  //----------------------------------------------------------------------------
  // CATCH AGE COMPOSITION
  //----------------------------------------------------------------------------

  for (int i=1;i<=myrs;i++) 
    {
      for (int j=1;j<=nages;j++)
        {
          if (obs_c_comp(i+md_offset,1)<0)
            {
              res_c_comp(i,j)=0;
            }
          else
            {
              res_c_comp(i,j)=obs_c_comp(i+md_offset,j)-sel_naa_prop(i,j);
            }
        }
    }

  //----------------------------------------------------------------------------
  // SPAWNING AGE COMPOSITION
  //----------------------------------------------------------------------------

  for (int i=1;i<=myrs;i++)
    {
     for (int j=1;j<=nages;j++)
       {
         if (obs_sp_comp(i+md_offset,j)<0)
           {
             res_sp_comp(i,j)=0;
           }
        else
          {
            res_sp_comp(i,j)=obs_sp_comp(i+md_offset,j)-est_sp_comp(i,j);
          }
       }
    }


  //----------------------------------------------------------------------------
  // EGG DEPOSITION
  //----------------------------------------------------------------------------

  for (int i=1;i<=myrs;i++)
    {
      if (tot_obs_egg(i+md_offset)<0){res_tot_egg(i)=0;}
      else{res_tot_egg(i)=log(tot_obs_egg(i+md_offset))-log(tot_est_egg(i));
      wSSQE(i)=lE(i+md_offset)*square(res_tot_egg(i));

    }}

  //----------------------------------------------------------------------------
  // SPAWNER-RECRUIT FUNCTION
  //----------------------------------------------------------------------------

  for (int i=1;i<=myrs-3;i++)
   {
     res_SR(i)=log(naa(i+3,1))-log(SR(i));
   }


  //----------------------------------------------------------------------------
  // MILE DAYS OF MILT
  //----------------------------------------------------------------------------

   for (int i=1;i<=myrs;i++)
    {
      if (mile_days(i+md_offset)<0)
        {
          M_DR(i)=0;
        }
      else
        {
          M_DR(i)=log(mile_days(i+md_offset))-log(M_D(i));
        }
      wSSQM(i)=lM(i+md_offset)*square(M_DR(i));
    }


  WSSQM=sum(wSSQM);
  SSQC=norm2(res_c_comp);
  SSQSp=norm2(res_sp_comp);
  WSSQE=sum(wSSQE);
  SSQR=norm2(res_SR);


FUNCTION evaluate_the_objective_function

  f=lC*SSQC+lS*SSQSp+WSSQE+lR*SSQR;


FUNCTION get_forecast

  for_naa(1)=alpha*tot_sp_B(myrs-2)*exp(-1.0*beta*tot_sp_B(myrs-2)); //forecast age 3 numbers
  
  for (int j=2;j<=nages-1;j++)
    {
      for_naa(j)=post_naa(myrs,j-1)*S_for;                           //forecast naa, ages 4 - 7
    }
    // Stupid!  No need for a loop for the plus group.
  for (int j=nages;j<=nages;j++)
    {
      for_naa(j)=post_naa(myrs,j-1)*S_for+post_naa(myrs,j)*S_for;    //forecast naa, age 8
    }


  for (int j=1;j<=nages;j++)
    {
      for_mat_naa(j)=for_naa(j)*Mat(myrs,j);    //forecast mature numbers at age
    }

  for_tot_mat_N=sum(for_mat_naa);

  for (int j=1;j<=nages;j++)
    {
      for_mat_baa(j)=for_mat_naa(j)*fw_a_a(j);  //forecast mature biomass at age
    }

  for_mat_B=sum(for_mat_baa);

  for (int j=1;j<=nages;j++)
    {
      for_mat_prop(j)=for_mat_naa(j)/for_tot_mat_N;  //forecast % mature at age
    }



  for (int j=1;j<=nages;j++)
    {
      for_mat_b_prop(j) = for_mat_baa(j)/for_mat_B;  //forecast % mature at age biomass
    }
  
  //GHL CALCS

  for_mat_B_st=for_mat_B/0.90718;                   //RETURN TO SHORT TONS
   HR_p=(2+8*for_mat_B_st/Thresh_denom)/100;        //Region-specific harvest rate (make sure all are in tons) CHECK!!!

  if(HR_p>0.2){HR=0.2;}
   else{HR=HR_p;}
  if(for_mat_B_st<Thresh){HR=0;}

  if(HR==0.2){GHL=0.2*for_mat_B_st;}
  else{GHL=(2+8*for_mat_B_st/Thresh_denom)/100*for_mat_B_st;}             //GHL calculation differs for Sitka
   
                                                  

FUNCTION get_FIGDATA
  FIGDATA.initialize();
  for (int i=1;i<=myrs;i++){
  for (int j=1;j<=1;j++){FIGDATA(i,j)=yminusthree(i);} //mature biomass (tons) (Figure 3)
  for (int j=2;j<=2;j++){FIGDATA(i,j)=yminustwo(i);} //mature biomass (tons) (Figure 3)
  for (int j=3;j<=3;j++){FIGDATA(i,j)=yminusone(i);} //mature biomass  (tons)(Figure 3)
  for (int j=4;j<=4;j++){FIGDATA(i,j)=tot_obs_egg(i+md_offset);}//observed total yearly egg deposition in trillions (Figure 1, Figure 2, Figure 3)
  for (int j=5;j<=5;j++){FIGDATA(i,j)=tot_est_egg(i);}//estimated total yearly egg deposition (Figure 1)
  for (int j=6;j<=6;j++){FIGDATA(i,j)=tot_mat_B(i);}// total mature biomass (metric tonnes) (Figure 2, Figure 3, Figure 6, Figure 14)
  for (int j=7;j<=7;j++){FIGDATA(i,j)=tcb(i+md_offset);}////observed catch biomass (tons)(Figure 2, Figure 3)
  for (int j=8;j<=8;j++){FIGDATA(i,j)=Threshold(i);}}//threshold vector (tons)(Figure 2, Figure 3)

  for (int i=1;i<=myrs-3;i++){
  for (int j=9;j<=9;j++){FIGDATA(i,j)=res_SR(i);}}//spawner-recruit residuals (Figure 4)

  for (int i=1;i<=myrs;i++){
  for (int j=10;j<=10;j++){FIGDATA(i,j)=res_tot_egg(i);}//egg deposition residuals (Figure 4)
  for (int j=11;j<=11;j++){FIGDATA(i,j)=init_age_3(i);}//initial age 3 (millions of fish) (Figure 5; Figure 10)
  for (int j=12;j<=12;j++){FIGDATA(i,j)=tot_sp_B(i);}// total spawning biomass (metric tons) (Figure 6; Figure 10, Figure 14)  
  for (int j=13;j<=13;j++){FIGDATA(i,j)=tot_sp_N(i);} // total spawning abundance [millions](Figure 6)
  for (int j=14;j<=14;j++){FIGDATA(i,j)=tot_mat_N(i);} // total mature abundance[millions] (Figure 6)
  for (int j=15;j<=15;j++){FIGDATA(i,j)=tot_post_N(i);}// total population [mature+ immature] - catch [millions](Figure 6)    
  for (int j=16;j<=16;j++){FIGDATA(i,j)=N(i);}// total abundance (mature + immature [millions] (Figure 6)
  for (int j=17;j<=22;j++){FIGDATA(i,j)=sel_naa_prop(i,j-16);}//proportion of N selected by gear (estimated) (Figure 8)
  for (int j=23;j<=28;j++){FIGDATA(i,j)=obs_c_comp(i+md_offset,j-22);}//observed catch compostion (Figure 8)
  for (int j=29;j<=34;j++){FIGDATA(i,j)=est_sp_comp(i,j-28);}//estimated spawner (cast net) age composition (Figure 9)
  for (int j=35;j<=40;j++){FIGDATA(i,j)=obs_sp_comp(i+md_offset,j-34);}}//observed spawner (cast net) age composition (Figure 9)

  for (int i=1;i<=myrs-3;i++){
  for (int j=41;j<=41;j++){FIGDATA(i,j)=SR(i);}}//(Figure 10)

  for (int i=1;i<=myrs;i++){
  for (int j=42;j<=42;j++){FIGDATA(i,j)=mod_yrs(i);}}//model years
FUNCTION output_FIGDATA

 ofstream figdata("FIGDATA.dat");
 figdata<<"yminusthree yminustwo yminusone tot_obs_egg tot_est_egg tot_mat_B tcb Threshold res_SR res_tot_egg init_age_3 tot_sp_B tot_sp_N tot_mat_N tot_post_N N sel_naa_prop3 sel_naa_prop4 sel_naa_prop5 sel_naa_prop6 sel_naa_prop7 sel_naa_prop8 obs_c_comp3 obs_c_comp4 obs_c_comp5 obs_c_comp6 obs_c_comp7 obs_c_comp8 est_sp_comp3 est_sp_comp4 est_sp_comp5 est_sp_comp6 est_sp_comp7 est_sp_comp8 obs_sp_comp3 obs_sp_comp4 obs_sp_comp5 obs_sp_comp6 obs_sp_comp7 obs_sp_comp8 SR Year"<<endl;
 figdata<<FIGDATA<<endl;

FUNCTION get_FIGDATAAGE
  FIGDATAAGE.initialize();
  for (int i=1;i<=nages;i++){
 //Mature biomass at age forecast
  for (int j=1;j<=1;j++){FIGDATAAGE(i,j)=for_mat_baa(i);}//forecasted mature biomass at age (metric tons) (Figure 11)
 //Mature numbers at age forecast (%)
  for (int j=2;j<=2;j++){FIGDATAAGE(i,j)=for_mat_prop(i);}//projected % mature #s at age (Figure 12)
 //forecast weight at age
  for (int j=3;j<=3;j++){FIGDATAAGE(i,j)=fw_a_a(i);}}//forecasted weight-at-age (Figure 13) 
  

FUNCTION output_FIGDATAAGE

 ofstream figdataage("FIGDATAAGE.dat");
 figdataage<<"for_mat_baa for_mat_prop fw_a_a"<<endl;
 figdataage<<FIGDATAAGE<<endl;  


FUNCTION compute_AICc
  n_tot_obs_egg.initialize();
  n_obs_c_comp.initialize();
  n_C.initialize();
  n_obs_sp_comp.initialize();
  n_S.initialize();
  n_R.initialize();
  n_d.initialize();
  n.initialize();
  sig_1.initialize();
  //sig_d.initialize();
  w_d.initialize();
  lnL_d.initialize();
  lnL.initialize();
  AIC.initialize();
  AICc.initialize();
  p.initialize();
 // n_My.initialize();
 // n_M.initialize();

  //Compute sample sizes
  
  //Egg Deposition
  for (int i=1;i<=myrs;i++)
    {
      n_tot_obs_egg(i)=1;
    }

  //SR function
  n_R=myrs-3;

  //Catch age comp
  for (int i=1;i<=myrs;i++)
    {
      for (int j=1;j<=nages;j++)
        {
          if (obs_c_comp(i+md_offset,j)<0)
            {
              n_obs_c_comp(i,j)=0;
            }
          else
            {
              n_obs_c_comp(i,j)=1;
            }
        }
    }

  n_C=sum(n_obs_c_comp);

  //Spawning age comp
  for (int i=1;i<=myrs;i++)
    {
      for (int j=1;j<=nages;j++)
        {
          if (obs_sp_comp(i+md_offset,1)<0)
            {
              n_obs_sp_comp(i,j)=0;
            }
          else
            {
              n_obs_sp_comp(i,j)=1;
            }
        }
     }

  n_S=sum(n_obs_sp_comp);


  //Mile-days
 // for (int i=1;i<=myrs;i++){
  //  if(mile_days(i+md_offset)>0){n_My(i)=1;}
   // else{n_My(i)=0;}}
 // n_M=sum(n_My);



  //Set up sample size vector
  n_d(1)=n_C;
  n_d(2)=n_S;
  n_d(3)=n_R;
 // n_d(3)=n_M;
  for (int i=1;i<=myrs;i++)
   {
     n_d(i+3)=n_tot_obs_egg(i);
   }

  n=sum(n_d);

  //Set up weighting vector
  w_d(1)=lC;
  w_d(2)=lS;
 // w_d(3)=lM;
  w_d(3)=lR;
  for (int i=1;i<=myrs;i++){
  w_d(i+3)=lE(i+md_offset);}

  //Set up sigma vector
  sig_1=f/n;
  //sig_d(1)=sig_1;
  //for (int i=1;i<=myrs+3;i++){
  //sig_d(i)=sig_1/w_d(i+1);}


  //Calculate log likelihood
  for (int i=1;i<=myrs+3;i++){ 
  if(w_d(i)>0){
   lnL_d(i)=-n_d(i)/2*(log(2*3.141593*sig_1/w_d(i))+1);}
  else{
   lnL_d(i)=0;}}
  lnL=sum(lnL_d);


  //Compute AICc
  p=initial_params::nvarcalc();
  AIC=-2*lnL+2*p;
  AICc=AIC+2*p*(p+1)/(n-p-1);


FUNCTION get_sdreports

  for (int i=1;i<=myrs;i++){
  SpawnBio(i)=tot_sp_B(i);}
  SpawnBioFor=for_mat_B_st;
  AICcsd=AICc;
  GHLsd=GHL;

FUNCTION get_report
       ofstream Report("Report.csv");

    int vsize = Year.size();

    Report << "MODEL RESULTS" <<endl;
 
    Report<<"Objective function value"<<","<<","<<","<<f<<endl;
    Report<<"AICc"<<","<<","<<","<<AICc<<endl;
    Report<<"  "<<endl;

    Report<<"Dataset Components:"<<endl;
    Report<<"Catch SSQ"<<","<<","<<","<<SSQC<<endl;
    Report<<"Spawning SSQ"<<","<<","<<","<<SSQSp<<endl;
    Report<<"Egg dep. SSQ"<<","<<","<<","<<WSSQE<<endl;
    Report<<"Ricker SSQ"<<","<<","<<","<<SSQR<<endl;
    Report<<"  "<<endl;
    Report<<"  "<<endl;

    Report<<"FORECAST INFORMATION:"<<endl;
    Report<<"Mature biomass forecast (short tons)"<<","<<","<<","<<SpawnBioFor<<endl;
    Report<<"Threshold (short tons)"<<","<<","<<","<<Thresh<<endl;
    Report<<"GHL (short tons)"<<","<<","<<","<<GHL<<endl;
    Report<<"Harvest Rate"<<","<<","<<","<<HR<<endl;
    Report<<"  "<<endl;
    Report<<"  "<<endl;

    Report<<","<<","<<","<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
    Report<<"Mature and immature numbers-at-age forecast (millions)"<<","<<","<<","<<","<<for_naa[1]<<","<<for_naa[2]<<","<<for_naa[3]<<","<<for_naa[4]<<","<<for_naa[5]<<","<<for_naa[6]<<endl;
    Report<<"Mature numbers-at-age forecast (millions)"<<","<<","<<","<<","<<for_mat_naa[1]<<","<<for_mat_naa[2]<<","<<for_mat_naa[3]<<","<<for_mat_naa[4]<<","<<for_mat_naa[5]<<","<<for_mat_naa[6]<<endl;
    Report<<"Proportion mature numbers-at-age forecast"<<","<<","<<","<<","<<for_mat_prop[1]<<","<<for_mat_prop[2]<<","<<for_mat_prop[3]<<","<<for_mat_prop[4]<<","<<for_mat_prop[5]<<","<<for_mat_prop[6]<<endl;
    Report<<"Mature biomass-at-age forecast (tons)"<<","<<","<<","<<","<<for_mat_baa[1]/0.90718<<","<<for_mat_baa[2]/0.90718<<","<<for_mat_baa[3]/0.90718<<","<<for_mat_baa[4]/0.90718<<","<<for_mat_baa[5]/0.90718<<","<<for_mat_baa[6]/0.90718<<endl;
    Report<<"Proportion mature biomass-at-age forecast"<<","<<","<<","<<","<<for_mat_b_prop[1]<<","<<for_mat_b_prop[2]<<","<<for_mat_b_prop[3]<<","<<for_mat_b_prop[4]<<","<<for_mat_b_prop[5]<<","<<for_mat_b_prop[6]<<endl;
    Report<<"Weight-at-age used in forecast (g)"<<","<<","<<","<<","<<fw_a_a[1]<<","<<fw_a_a[2]<<","<<fw_a_a[3]<<","<<fw_a_a[4]<<","<<fw_a_a[5]<<","<<fw_a_a[6]<<endl;
    Report<<"Maturity used for forecast"<<","<<","<<","<<","<<mat_for[1]<<","<<mat_for[2]<<","<<mat_for[3]<<","<<mat_for[4]<<","<<mat_for[5]<<","<<mat_for[6]<<endl;
    Report<<"Survival forecast"<<","<<","<<","<<","<<S_for<<","<<S_for<<","<<S_for<<","<<S_for<<","<<S_for<<","<<S_for<<endl;
    Report<<"  "<<endl;
    Report<<"  "<<endl;

   Report<<"Year"<<","<<"Mature biomass (tons)"<<","<<"Spawning biomass (tons)"<<","<<"Catch (tons)"<<","<<"Mature and immature age-3 abundance (recruitment in millions)"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<tot_mat_B[n+1]/0.90718<<","<<tot_sp_B[n+1]/0.90718<<","<<tcb[n+1+md_offset]<<","<<init_age_3[n+1]<<endl;
   Report<<"  "<<endl;


   Report<<"Mature and immature numbers-at-age (millions)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<naa(n+1,1)<<","<<naa(n+1,2)<<","<<naa(n+1,3)<<","<<naa(n+1,4)<<","<<naa(n+1,5)<<","<<naa(n+1,6)<<endl;
   Report<<"  "<<endl;

   Report<<"Mature numbers-at-age (millions)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_mat_naa(n+1,1)<<","<<est_mat_naa(n+1,2)<<","<<est_mat_naa(n+1,3)<<","<<est_mat_naa(n+1,4)<<","<<est_mat_naa(n+1,5)<<","<<est_mat_naa(n+1,6)<<endl;
   Report<<"  "<<endl;

   Report<<"Spawning numbers-at-age (millions)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_sp_naa(n+1,1)<<","<<est_sp_naa(n+1,2)<<","<<est_sp_naa(n+1,3)<<","<<est_sp_naa(n+1,4)<<","<<est_sp_naa(n+1,5)<<","<<est_sp_naa(n+1,6)<<endl;
   Report<<"  "<<endl;

   Report<<"Estimated catch-at-age (millions)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_c_naa(n+1,1)<<","<<est_c_naa(n+1,2)<<","<<est_c_naa(n+1,3)<<","<<est_c_naa(n+1,4)<<","<<est_c_naa(n+1,5)<<","<<est_c_naa(n+1,6)<<endl;
   Report<<"  "<<endl;

   Report<<"Estimated proportion of each age class that is caught (e.g. number of age-3 fish caught divided by the number of age-3 mature at age)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_c_naa(n+1,1)/est_sp_naa(n+1,1)<<","<<est_c_naa(n+1,2)/est_sp_naa(n+1,2)<<","<<est_c_naa(n+1,3)/est_sp_naa(n+1,3)<<","<<est_c_naa(n+1,4)/est_sp_naa(n+1,4)<<","<<est_c_naa(n+1,5)/est_sp_naa(n+1,5)<<","<<est_c_naa(n+1,6)/est_sp_naa(n+1,6)<<endl;
   Report<<"  "<<endl;

   Report<<"Mature biomass-at-age"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_mat_baa(n+1,1)/0.90718<<","<<est_mat_baa(n+1,2)/0.90718<<","<<est_mat_baa(n+1,3)/0.90718<<","<<est_mat_baa(n+1,4)/0.90718<<","<<est_mat_baa(n+1,5)/0.90718<<","<<est_mat_baa(n+1,6)/0.90718<<endl;
   Report<<"  "<<endl;

   Report<<"Spawning biomass-at-age (tons)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_sp_baa(n+1,1)/0.90718<<","<<est_sp_baa(n+1,2)/0.90718<<","<<est_sp_baa(n+1,3)/0.90718<<","<<est_sp_baa(n+1,4)/0.90718<<","<<est_sp_baa(n+1,5)/0.90718<<","<<est_sp_baa(n+1,6)/0.90718<<endl;
   Report<<"  "<<endl;

   Report<<"Survival"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<Sur(n+1,1)<<","<<Sur(n+1,2)<<","<<Sur(n+1,3)<<","<<Sur(n+1,4)<<","<<Sur(n+1,5)<<","<<Sur(n+1,6)<<endl;
   Report<<"  "<<endl;

   Report<<"Maturity"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<Mat(n+1,1)<<","<<Mat(n+1,2)<<","<<Mat(n+1,3)<<","<<Mat(n+1,4)<<","<<Mat(n+1,5)<<","<<Mat(n+1,6)<<endl;
   Report<<"  "<<endl;

   Report<<"Gear Selectivity"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<GS_Sc(n+1,1)<<","<<GS_Sc(n+1,2)<<","<<GS_Sc(n+1,3)<<","<<GS_Sc(n+1,4)<<","<<GS_Sc(n+1,5)<<","<<GS_Sc(n+1,6)<<endl;
   Report<<"  "<<endl;
//change based on the stock
   Report<<"Observed spawning age-composition (cast net)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<obs_sp_comp(n+10,1)<<","<<obs_sp_comp(n+10,2)<<","<<obs_sp_comp(n+10,3)<<","<<obs_sp_comp(n+10,4)<<","<<obs_sp_comp(n+10,5)<<","<<obs_sp_comp(n+10,6)<<endl;
   Report<<"  "<<endl;
//change based on the stock
   Report<<"Observed commercial seine age-composition (spring seine)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<obs_c_comp(n+10,1)<<","<<obs_c_comp(n+10,2)<<","<<obs_c_comp(n+10,3)<<","<<obs_c_comp(n+10,4)<<","<<obs_c_comp(n+10,5)<<","<<obs_c_comp(n+10,6)<<endl;
   Report<<"  "<<endl;

   Report.close();
 
REPORT_SECTION
   //The reports below have to included in the output or the R figures will not work properly
  //Note: you can't name the same variable two different names or else the reptoRlist function in program R will give an error
  report<<"Scaled_Gear_Selectivity"<<endl;//needs to be in report for figures (Figure 7)
  report<<GS_Sc<<endl;
  report<<"Maturity"<<endl;//needs to be in report for figures (Figure 7)
  report<<Mat<<endl;
  report<<"Survival"<<endl;//needs to be in report for figures (Figure 7)
  report<<Sur<<endl;
  report<<"GHL"<<endl;//needs to be in report for figures (Figure 14)
  report<<GHL<<endl;
  report<<"yminusthreeFOR"<<endl;//needs to be in report for figures-make sure in tons (Figure 3)
  report<<yminusthreeFOR<<endl;
  report<<"yminustwoFOR"<<endl;//needs to be in report for figures-make sure in tons (Figure 3)
  report<<yminustwoFOR<<endl;
  report<<"yminusoneFOR"<<endl;//needs to be in report for figures-make sure in tons (Figure 3)
  report<<yminusoneFOR<<endl;
  report<<"for_mat_B_st"<<endl;//needs to be in report for figures-make sure in tons (Figure 2, Figure 3, Figure 14)
  report<<for_mat_B_st<<endl;
  report<<"Thresh"<<endl;//current threshold-make sure in tons//needs to be in report for figures (Figure 14)
  report<<Thresh<<endl;
  report<<"res_sp_comp"<<endl;//needs to be in report for figures (Figure 15)
  report<<res_sp_comp<<endl;
  report<<"est_sp_comp"<<endl;//needs to be in report for figures (Figure 16)
  report<<est_sp_comp<<endl;
  report<<"res_c_comp"<<endl;//needs to be in report for figures (Figure 17)
  report<<res_c_comp<<endl;
  report<<""<<endl;


  report<<"MODEL_RESULTS:"<<endl;
  report<<"Objective_function_value"<<endl;
  report<<f<<endl;
  report<<"Dataset_components:"<<endl;
  report<<"Catch SSQ"<<endl;
  report<<SSQC<<endl;
  report<<"Spawning_SSQ"<<endl;
  report<<SSQSp<<endl;
  report<<"Egg_dep_SSQ"<<endl;
  report<<WSSQE<<endl;
  report<<"Mile_Days SSQ"<<endl;;
  report<<WSSQM<<endl;
  report<<"Ricker_SSQ"<<endl;
  report<<SSQR<<endl;


  report<<"FORECAST INFORMATION:"<<endl;
  report<<"Harvest Rate"<<endl;
  report<<HR<<endl;
  report<<"Numbers-at-age forecast"<<endl;
  report<<for_naa<<endl;
  report<<"Mature biomass-at-age forecast"<<endl;
  report<<for_mat_baa<<endl;
  report<<"Mature biomass at age"<<endl;
  report<<est_mat_baa<<endl;
  report<<"spawning biomass at age"<<endl;
  report<<est_sp_baa<<endl;


RUNTIME_SECTION
  maximum_function_evaluations 5000 5000 5000 5000
  convergence_criteria 0.0001


TOP_OF_MAIN_SECTION
  arrmblsize=5000000;
  gradient_structure::set_MAX_NVAR_OFFSET(5000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(800000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(800000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
	#include <admodel.h>
	#include <string.h>
	#include <time.h>
        adstring model_name;
        adstring data_file;


	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

	adstring BaseFileName;
	adstring ReportFileName;
	adstring NewFileName;

	adstring stripExtension(adstring fileName)
	{
		/*
		This function strips the file extension
		from the fileName argument and returns
		the file name without the extension.
		*/
		const int length = fileName.size();
		for (int i=length; i>=0; --i)
		{
			if (fileName(i)=='.')
			{
				return fileName(1,i-1);
			}
		}
		return fileName;
	}

FINAL_SECTION

	//  Print run time statistics to the screen.
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
        cout<<""<<endl;
	cout<<"--Objective function value: "<<f<<endl;
        cout<<""<<endl;
	cout<<"--AIC: "<<AICc<<endl;
        cout<<""<<endl;
	cout<<"--Maximum gradient component: "<<objective_function_value::gmax<<endl;
        cout<<""<<endl;
	cout<<"*******************************************"<<endl;

        cout<<""<<endl;
        cout<< "O frabjous day!"<<endl;
        cout<< "The sheep frolic!"<<endl;
        cout<<""<<endl;
        cout<<"        ...moo..."<<endl;
        cout<<"            | "<<endl;
        cout<<"            | "<<endl;
        cout<<"            | "<<endl;
        cout<<"             _.%%%%%%%%%%%%%             " <<endl;
        cout<<"            //-_%%%%%%%%%%%%%            " <<endl;
        cout<<"           (_ %\\%%%%%%%%%%%%%%~             "<<endl;
        cout<<"               %%%%%%%%%%%%%%             "<<endl;
        cout<<"                 %%%%%*%%%%              "<<endl;
        cout<<"            ,,,,,,||,,,,||,,,,,         "<<endl;
        cout<<""<<endl;


		cout<<"            ▄▄▄▄▄▄▄▄▄▄▄▄▄            "<<endl;
		cout<<"         ▄▀▀═════════════▀▀▄         "<<endl;
		cout<<"        █═══════════════════█        "<<endl;
		cout<<"       █═════════════════════█       "<<endl;
		cout<<"      █═══▄▄▄▄▄▄▄═══▄▄▄▄▄▄▄═══█      "<<endl;
		cout<<"     █═══█████████═█████████═══█     "<<endl;
		cout<<"     █══██▀    ▀█████▀    ▀██══█     "<<endl;
		cout<<"    ██████  █▀█  ███  █▀█  ██████    "<<endl;
		cout<<"    ██████  ▀▀▀  ███  ▀▀▀  ██████    "<<endl;
		cout<<"     █══▀█▄    ▄██ ██▄    ▄█▀══█     "<<endl;
		cout<<"     █════▀█████▀   ▀█████▀════█     "<<endl;
		cout<<"     █═════════════════════════█     "<<endl;
		cout<<"     █═════════════════════════█     "<<endl;
		cout<<"     █═══════▀▄▄▄▄▄▄▄▄▄▀═══════█     "<<endl;
		cout<<"     █═════════════════════════█     "<<endl;
		cout<<"    ▐▓▓▌═════════════════════▐▓▓▌    "<<endl;
		cout<<"    ▐▐▓▓▌▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▐▓▓▌▌    "<<endl;
		cout<<"    █══▐▓▄▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▄▓▌══█    "<<endl;
		cout<<"   █══▌═▐▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▌═▐══█   "<<endl;
		cout<<"   █══█═▐▓▓▓▓▓▓▄▄▄▄▄▄▄▓▓▓▓▓▓▌═█══█   "<<endl;
		cout<<"   █══█═▐▓▓▓▓▓▓▐██▀██▌▓▓▓▓▓▓▌═█══█   "<<endl;
		cout<<"   █══█═▐▓▓▓▓▓▓▓▀▀▀▀▀▓▓▓▓▓▓▓▌═█══█   "<<endl;
		cout<<"   █══█▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓█══█   "<<endl;
		cout<<"  ▄█══█▐▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▌█══█▄  "<<endl;
		cout<<"  █████▐▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓█ █████  "<<endl;
		cout<<"  ██████▐▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▌███████  "<<endl;
		cout<<"   ▀█▀█  ▐▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▌   █▀█▀   "<<endl;
		cout<<"          ▐▓▓▓▓▓▓▌▐▓▓▓▓▓▓▌           "<<endl;
		cout<<"           ▐▓▓▓▓▌  ▐▓▓▓▓▌            "<<endl;
		cout<<"          ▄████▀    ▀████▄           "<<endl;
		cout<<"          ▀▀▀▀        ▀▀▀▀           "<<endl;

	
