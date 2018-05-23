#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVectorT.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>

using namespace std;

//TODO: implement HacR_alignAGL / HacL_alignAGL which seems to have the epics data
//	about the spectrometer angles. This is also stored from the E tree
//	But this is somehow NOT what is in the db_run.dat file?? This is the
//	EPICs value from the Hall A gen tools value??
//	WHICH ONE DOES THE ANALYZER USE?? -- I bet the db_run.dat file... 

// -------------------------------------------------------
// Defining values for cuts
const double TG_Theta_Max =  0.06 ; // Acceptance cut lim theta 40mrad
const double TG_Theta_Min = -0.06 ; // Acceptance cut lim theta 40mrad
const double TG_Phi_Max   =  0.030; // Acceptance cut lim phi   25mrad (in-plane?)
const double TG_Phi_Min   = -0.030; // Acceptance cut lim phi   25mrad (in-plane?)
const double TG_Dp_Max    =  0.045; // Acceptance cut lim delta 4.5%
const double TG_Dp_Min    = -0.045; // Acceptance cut lim delta 4.5%
const double TG_VZ_Max    =  0.10 ; // Vertex cut lim
const double TG_VZ_Min    = -0.10 ; // Vertex cut lim
const double GC_Cut_L     = 2000. ; // LHRS Gas Cherenkov ADC sum cut lim
const double EP_Cut_L     = 0.80  ; // LHRS Energy/Momentum cut lim
const double EC2_Cut_L    = 300   ; // LHRS Energy deposited cut lim
const double GC_Cut_R     = 2000. ; // RHRS Gas Cherenkov ADC sum cut lim
const double EP_Cut_R     = 0.80  ; // RHRS Energy/Momentum cut lim
const double EC2_Cut_R    = 300   ; // RHRS Energy deposited cut lim
const int Main_Trigger    = 1     ; // Main trigger cut value

// Physics constants
const double mP   = 0.938  ;	// Proton
const double mD   = 1.876  ;	// Deuteron
const double mHe3 = 2.80941;	// Helium-3
const double mH3  = 2.80943;	// Tritium
const double mC   = 11.172 ;	// Carbon
const double mAl  = 25.118 ;	// Aluminum

// Forward declaring functions
double get_kin(const char * kin, double &e_th, double &e_p, double &p_th, double &p_p);
double get_mTarg(const char * filename); // Calculates the target mass from the target position
double Calculate_Charge_from_Scaler(const char * filename);
double Calculate_Time(const char * filename);

// ================================================================================================
int main(int argc, char ** argv)
{
	if ( argc<5 )
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "   skimmer [run number] [1H/2H/3H/3He] [fast/slow/mid] /path/to/first/file [/path/to/second/file...]\n\n";
		return -1;
	}

	int runNumber=atoi(argv[1]);

	// Establish the kinematics
	double e_theta_c, e_mom_c, p_theta_c, p_mom_c;
	get_kin(argv[3],e_theta_c,e_mom_c,p_theta_c,p_mom_c);

	TVectorT<double> kinematics(5);
	kinematics[0]=e_theta_c;
	kinematics[1]=e_mom_c;
	kinematics[2]=p_theta_c;
	kinematics[3]=p_mom_c;

	// Pull the target type 
	const double mTarg = get_mTarg(argv[4]);

	// Loop over list of runs and chain the files
	double totalQ_mC=0.;
	cout << "*************************************************\n";
	cout << "Here are the files that will be skimmed:\n";
	TChain *T = new TChain("T");
	TChain *E = new TChain("E");
	for (int i=4 ; i<argc ; i++)
	{
		double thisCharge = Calculate_Charge_from_Scaler(argv[i]);
		if (thisCharge > totalQ_mC)
			totalQ_mC = thisCharge;

		cout << "     " << argv[i] << "\n";
		T->Add(argv[i]);
		E->Add(argv[i]);
	}
	cout << "*************************************************\n";
	cout << "Target mass taken as: " << mTarg << " GeV/c^2\n";

	// Do the charge vector
	TVectorT<double> totalCharge(1);
	totalCharge[0] = totalQ_mC;
	double testCharge = 0;
	// Establish memory for input tree
	const int maxNTracks=100;
	double 	e_cer[maxNTracks]		;
	double 	e_prl1[maxNTracks]		,	e_prl2[maxNTracks]	,
		e_ntrk[maxNTracks]		,	p_ntrk[maxNTracks]	,
		e_ytar[maxNTracks]		,	p_ytar[maxNTracks]	,
		e_yptar[maxNTracks]		,	p_yptar[maxNTracks]	,
		gold_e_yptar[maxNTracks]	,	gold_p_yptar[maxNTracks],
		e_xptar[maxNTracks]		,	p_xptar[maxNTracks]	,
		gold_e_xptar[maxNTracks]	,	gold_p_xptar[maxNTracks],
		e_delta[maxNTracks]		,	p_delta[maxNTracks]	,
		gold_e_delta[maxNTracks]	,	gold_p_delta[maxNTracks],
		e_mom[maxNTracks]		,	p_mom[maxNTracks]	,
		gold_e_mom[maxNTracks]		,	gold_p_mom[maxNTracks]	,
		e_beta[maxNTracks]		,	p_beta[maxNTracks]	,
		gold_e_beta[maxNTracks]		,	gold_p_beta[maxNTracks]	,
		e_x[maxNTracks]			,	p_x[maxNTracks]		,
		e_y[maxNTracks]			,	p_y[maxNTracks]		,
		e_z[maxNTracks]			,	p_z[maxNTracks]		,
		e_px[maxNTracks]		,	p_px[maxNTracks]	,
		e_py[maxNTracks]		,	p_py[maxNTracks]	,
		e_pz[maxNTracks]		,	p_pz[maxNTracks]	,
		gold_e_px[maxNTracks]		,	gold_p_px[maxNTracks]	,
		gold_e_py[maxNTracks]		,	gold_p_py[maxNTracks]	,
		gold_e_pz[maxNTracks]		,	gold_p_pz[maxNTracks]	,
		t1[maxNTracks]			,	t4[maxNTracks]		,
		tcoinc[maxNTracks]		,
		Pmiss[maxNTracks]		,	Emiss[maxNTracks]	,
		Pmiss_x[maxNTracks]		,	Pmiss_y[maxNTracks]	,
		Pmiss_z[maxNTracks]		,	p_thWe[maxNTracks]	,
		ph_bq[maxNTracks]		,	th_bq[maxNTracks]	,
		ph_xq[maxNTracks]		,	th_xq[maxNTracks]	,
		rast_Pmiss[maxNTracks]		,	rast_Emiss[maxNTracks]	,
		rast_Pmiss_x[maxNTracks]	,	rast_Pmiss_y[maxNTracks],
		rast_Pmiss_z[maxNTracks]	,	p_rast_thWe[maxNTracks] ,
		rast_ph_bq[maxNTracks]		,	rast_th_bq[maxNTracks]	,
		rast_ph_xq[maxNTracks]		,	rast_th_xq[maxNTracks]	,
		exRa_Pmiss[maxNTracks]		,	exRa_Emiss[maxNTracks]	,
		exRa_Pmiss_x[maxNTracks]	,	exRa_Pmiss_y[maxNTracks],
		exRa_Pmiss_z[maxNTracks]	,	p_exRa_thWe[maxNTracks] ,
		exRa_ph_bq[maxNTracks]		,	exRa_th_bq[maxNTracks]	,
		exRa_ph_xq[maxNTracks]		,	exRa_th_xq[maxNTracks]	,
		Q2[maxNTracks]			,	W2[maxNTracks]		,
		Nu[maxNTracks]			,	ph_q[maxNTracks]	,
		th_q[maxNTracks]		,	xB[maxNTracks]		,
		q3m[maxNTracks]			,	q_x[maxNTracks]		,
		q_y[maxNTracks]			,	q_z[maxNTracks]		,
		e_th[maxNTracks]		,
		rast_Q2[maxNTracks]		,	rast_W2[maxNTracks]	,
		rast_Nu[maxNTracks]		,	rast_ph_q[maxNTracks]	,
		rast_th_q[maxNTracks]		,	rast_xB[maxNTracks]	,
		rast_q3m[maxNTracks]		,	rast_q_x[maxNTracks]	,
		rast_q_y[maxNTracks]		,	rast_q_z[maxNTracks]	,		
		e_rast_th[maxNTracks]		,
		exRa_Q2[maxNTracks]		,	exRa_W2[maxNTracks]	,
		exRa_Nu[maxNTracks]		,	exRa_ph_q[maxNTracks]	,
		exRa_th_q[maxNTracks]		,	exRa_xB[maxNTracks]	,
		exRa_q3m[maxNTracks]		,	exRa_q_x[maxNTracks]	,
		exRa_q_y[maxNTracks]		,	exRa_q_z[maxNTracks]	,
		e_exRa_th[maxNTracks]		,
		e_ext_deltaDp[maxNTracks]	,	e_ext_deltaP[maxNTracks],
		e_ext_deltaTh[maxNTracks]	,	e_ext_delta[maxNTracks]	,
		e_ext_mom[maxNTracks]		,	e_ext_yptar[maxNTracks]	,
		e_ext_xptar[maxNTracks]		,	e_ext_px[maxNTracks]	,
		e_ext_py[maxNTracks]		,	e_ext_pz[maxNTracks]	,		
		p_ext_deltaDp[maxNTracks]	,	p_ext_deltaP[maxNTracks],
		p_ext_deltaTh[maxNTracks]	,	p_ext_delta[maxNTracks]	,
		p_ext_mom[maxNTracks]		,	p_ext_yptar[maxNTracks]	,
		p_ext_xptar[maxNTracks]		,	p_ext_px[maxNTracks]	,
		p_ext_py[maxNTracks]		,	p_ext_pz[maxNTracks]	,		
		BCMcharge[maxNTracks]		,	BCMcurr[maxNTracks]	,
		BCMrenew[maxNTracks]		;


	// Checked against analyzer output?
	double testE;
	// Pre-shower pion rej.
	T->SetBranchAddress("L.prl1.e"    , e_prl1 );			// checked
	// Show pion rej.
	T->SetBranchAddress("L.prl2.e"    , e_prl2 );			// checked
	// Gas cherenkov
	T->SetBranchAddress("L.cer.asum_c", e_cer  );
	// tr = Track information
	// track number
	T->SetBranchAddress("L.tr.n"      , e_ntrk );
	T->SetBranchAddress("R.tr.n"      , p_ntrk );
	// target y coordinate
	T->SetBranchAddress("L.tr.tg_y"   , e_ytar );			// checked
	T->SetBranchAddress("R.tr.tg_y"   , p_ytar );			// checked
	// tangent of target phi angle
	T->SetBranchAddress("L.tr.tg_ph"  , e_yptar);			// checked
	T->SetBranchAddress("R.tr.tg_ph"  , p_yptar);			// checked
	T->SetBranchAddress("L.gold.ph"	  , gold_e_yptar );
	T->SetBranchAddress("R.gold.ph"	  , gold_p_yptar );
	// tangent of target theta angle
	T->SetBranchAddress("L.tr.tg_th"  , e_xptar);			// checked
	T->SetBranchAddress("R.tr.tg_th"  , p_xptar);			// checked
	T->SetBranchAddress("L.gold.th"	  , gold_e_xptar );
	T->SetBranchAddress("R.gold.th"	  , gold_p_xptar );
	// target delta
	T->SetBranchAddress("L.tr.tg_dp"  , e_delta);			// checked
	T->SetBranchAddress("R.tr.tg_dp"  , p_delta);			// checked
	T->SetBranchAddress("L.gold.dp"	  , gold_e_delta );
	T->SetBranchAddress("R.gold.dp"	  , gold_p_delta );
	// track momentum
	T->SetBranchAddress("L.tr.p"      , e_mom  );			// checked
	T->SetBranchAddress("R.tr.p"      , p_mom  );			// checked
	T->SetBranchAddress("L.gold.p"      , gold_e_mom  );
	T->SetBranchAddress("R.gold.p"      , gold_p_mom  );
	// track beta
	T->SetBranchAddress("L.tr.beta"   , e_beta );			// checked
	T->SetBranchAddress("R.tr.beta"   , p_beta );			// checked
	T->SetBranchAddress("L.gold.beta"      , gold_e_beta  );
	T->SetBranchAddress("R.gold.beta"      , gold_p_beta  );
	// vertex (z)
	T->SetBranchAddress("L.tr.vx"     , e_x    );			// checked
	T->SetBranchAddress("L.tr.vy"     , e_y    );			// checked
	T->SetBranchAddress("L.tr.vz"     , e_z    );			// checked
	T->SetBranchAddress("R.tr.vx"     , p_x    );			// checked
	T->SetBranchAddress("R.tr.vy"     , p_y    );			// checked
	T->SetBranchAddress("R.tr.vz"     , p_z    );			// checked
	// track momentum components
	T->SetBranchAddress("L.tr.px"     , e_px   );			// checked	
	T->SetBranchAddress("L.tr.py"     , e_py   );			// checked
	T->SetBranchAddress("L.tr.pz"     , e_pz   );			// checked
	T->SetBranchAddress("R.tr.px"     , p_px   );			// checked
	T->SetBranchAddress("R.tr.py"     , p_py   );			// checked
	T->SetBranchAddress("R.tr.pz"     , p_pz   );			// checked
	T->SetBranchAddress("L.gold.px"     , gold_e_px   );
	T->SetBranchAddress("L.gold.py"     , gold_e_py   );
	T->SetBranchAddress("L.gold.pz"     , gold_e_pz   );
	T->SetBranchAddress("R.gold.px"     , gold_p_px   );
	T->SetBranchAddress("R.gold.py"     , gold_p_py   );
	T->SetBranchAddress("R.gold.pz"     , gold_p_pz   );
	// There are also variables like 
	// R.tr.th and R.tr.ph  and these
	// ones are at the VDC plane, not at the 
	// target, so I don't think it's so
	// useful...
	// Now these are some time variables
	T->SetBranchAddress("DR.rrawt1"   , t1     );			// checked
	T->SetBranchAddress("DR.rrawt7"   , t4     );			// checked
	T->SetBranchAddress("DR.rrawt4"   , tcoinc );			// checked
	// All of the EKL/EKR variables are IDEAL beam assumed
	// 	you can access idealbeam with branch 'ib'
	// 	some members include x/y/px/py/pz/th/ph
	// 	but it just points beam along (0,0,1)
	// All of the EKLc/EKRc are corrected for rastered beam
	// 	you can access rasterbeam with branch 'Lrb'/'Rrb'
	// 	rastered beam has pointing information of the beam!!!
	// 	I'm not sure if this is actually implemented in the EKLc
	// 	class but I would hope so, otherwise that's really stupid.
	// 	It has things like x/y/dir.x/dir.y so I'm assuming
	// 	you get beam pointing.
	// All of the EKLx/EKRx are corrected for rastered beam and extended target
	// 	you can access the reaction point with branch 'rpl'/'rpr'
	// 	which uses the rastered beam information. then this is
	// 	passed into the exL/exR extended target class for LHRS.
	// 	finally these are both passed to the EKLx/EKRx, which in principle
	// 	will correct variables for the extended target. The problem is
	// 	that you don't have access to theta/phi angles via this
	// 	class. and so you need to grab corrections to track from
	// 	the extended target class. The extended target class has all the
	// 	theta and phi corrections you need to apply to the track class
	// Right arm IDEAL BEAM class
	T->SetBranchAddress("EKR.pmiss"    , Pmiss  );			// checked
	T->SetBranchAddress("EKR.pmiss_x"  , Pmiss_x);			// checked
	T->SetBranchAddress("EKR.pmiss_y"  , Pmiss_y);			// checked
	T->SetBranchAddress("EKR.pmiss_z"  , Pmiss_z);			// checked
	T->SetBranchAddress("EKR.emiss"    , Emiss);			// checked
	T->SetBranchAddress("EKR.ph_bq"    , ph_bq);			// checked
	T->SetBranchAddress("EKR.th_bq"    , th_bq);			// checked
	T->SetBranchAddress("EKR.ph_xq"    , ph_xq);			// checked
	T->SetBranchAddress("EKR.th_xq"    , th_xq);			// checked
	T->SetBranchAddress("EKR.xangle"   , p_thWe); 			// checked
	// Right arm raster-corrected class
	T->SetBranchAddress("EKRc.pmiss"   , rast_Pmiss  );
	T->SetBranchAddress("EKRc.pmiss_x" , rast_Pmiss_x);
	T->SetBranchAddress("EKRc.pmiss_y" , rast_Pmiss_y);
	T->SetBranchAddress("EKRc.pmiss_z" , rast_Pmiss_z);
	T->SetBranchAddress("EKRc.emiss"   , rast_Emiss);
	T->SetBranchAddress("EKRc.xangle"  , p_rast_thWe );
	T->SetBranchAddress("EKRc.ph_bq"   , rast_ph_bq);
	T->SetBranchAddress("EKRc.th_bq"   , rast_th_bq);
	T->SetBranchAddress("EKRc.ph_xq"   , rast_ph_xq);
	T->SetBranchAddress("EKRc.th_xq"   , rast_th_xq);
	// Right arm raster-corrected and extended-target corrected class
	T->SetBranchAddress("EKRx.pmiss"   , exRa_Pmiss  );		// checked
	T->SetBranchAddress("EKRx.pmiss_x" , exRa_Pmiss_x);		// checked
	T->SetBranchAddress("EKRx.pmiss_y" , exRa_Pmiss_y);		// checked
	T->SetBranchAddress("EKRx.pmiss_z" , exRa_Pmiss_z);		// checked
	T->SetBranchAddress("EKRx.emiss"   , exRa_Emiss);		// checked
	T->SetBranchAddress("EKRx.xangle"  , p_exRa_thWe );		// checked
	T->SetBranchAddress("EKRx.ph_bq"   , exRa_ph_bq);		// checked
	T->SetBranchAddress("EKRx.th_bq"   , exRa_th_bq);		// checked
	T->SetBranchAddress("EKRx.ph_xq"   , exRa_ph_xq);		// checked
	T->SetBranchAddress("EKRx.th_xq"   , exRa_th_xq);		// checked
	// Left arm IDEAL BEAM class
	T->SetBranchAddress("EKL.Q2"	   , Q2	      );		// checked
	T->SetBranchAddress("EKL.W2"	   , W2	      );		// checked
	T->SetBranchAddress("EKL.nu"	   , Nu	      );		// checked
	T->SetBranchAddress("EKL.ph_q"	   , ph_q     );		// checked
	T->SetBranchAddress("EKL.th_q"	   , th_q     );		// checked
	T->SetBranchAddress("EKL.x_bj"	   , xB       );		// checked
	T->SetBranchAddress("EKL.q3m"      , q3m      );		// checked
	T->SetBranchAddress("EKL.q_x"	   , q_x      );		// checked
	T->SetBranchAddress("EKL.q_y"	   , q_y      );		// checked
	T->SetBranchAddress("EKL.q_z"	   , q_z      );		// checked
	T->SetBranchAddress("EKL.angle"    , e_th     );		
	// Left arm raster-corrected class
	T->SetBranchAddress("EKLc.Q2"	   , rast_Q2	      );
	T->SetBranchAddress("EKLc.W2"	   , rast_W2	      );
	T->SetBranchAddress("EKLc.nu"	   , rast_Nu	      );
	T->SetBranchAddress("EKLc.ph_q"	   , rast_ph_q     );
	T->SetBranchAddress("EKLc.th_q"	   , rast_th_q     );
	T->SetBranchAddress("EKLc.x_bj"	   , rast_xB       );
	T->SetBranchAddress("EKLc.q3m"     , rast_q3m      );
	T->SetBranchAddress("EKLc.q_x"	   , rast_q_x      );
	T->SetBranchAddress("EKLc.q_y"	   , rast_q_y      );
	T->SetBranchAddress("EKLc.q_z"	   , rast_q_z      );
	T->SetBranchAddress("EKLc.angle"   , e_rast_th	   );		
	// Left arm raster-corrected and extended-target corrected class
	T->SetBranchAddress("EKLx.Q2"	   , exRa_Q2	   );		// checked
	T->SetBranchAddress("EKLx.W2"	   , exRa_W2	   );		// checked
	T->SetBranchAddress("EKLx.nu"	   , exRa_Nu	   );		// checked
	T->SetBranchAddress("EKLx.ph_q"	   , exRa_ph_q     );		// checked
	T->SetBranchAddress("EKLx.th_q"	   , exRa_th_q     );		// checked
	T->SetBranchAddress("EKLx.x_bj"	   , exRa_xB       );		// checked
	T->SetBranchAddress("EKLx.q3m"     , exRa_q3m      );		// checked
	T->SetBranchAddress("EKLx.q_x"	   , exRa_q_x      );		// checked
	T->SetBranchAddress("EKLx.q_y"	   , exRa_q_y      );		// checked
	T->SetBranchAddress("EKLx.q_z"	   , exRa_q_z      );		// checked
	T->SetBranchAddress("EKLx.angle"   , e_exRa_th	   );		// checked		
	// Corrections from extended target -- to be used on:
	// 	exL.delta_th + L.tr.tg_th  == exL.th
	//   and etc...
	T->SetBranchAddress("exL.delta_dp", 	e_ext_deltaDp  ); 	// checked
	T->SetBranchAddress("exL.delta_p" , 	e_ext_deltaP   );	// checked
	T->SetBranchAddress("exL.delta_th" ,	e_ext_deltaTh  );	// checked
	T->SetBranchAddress("exL.dp" , 		e_ext_delta    );	// checked
	T->SetBranchAddress("exL.p" 	,  	e_ext_mom      );	// checked
	T->SetBranchAddress("exL.ph" , 		e_ext_yptar    );	// checked
	T->SetBranchAddress("exL.th" , 		e_ext_xptar    );	// checked
	T->SetBranchAddress("exL.px" , 		e_ext_px       );	// checked
	T->SetBranchAddress("exL.py" , 		e_ext_py       );	// checked
	T->SetBranchAddress("exL.pz" , 		e_ext_pz       );	// checked
	// and for the right arm...
	T->SetBranchAddress("exR.delta_dp", 	p_ext_deltaDp  ); 	// checked
	T->SetBranchAddress("exR.delta_p" , 	p_ext_deltaP   );	// checked
	T->SetBranchAddress("exR.delta_th" ,	p_ext_deltaTh  );	// checked
	T->SetBranchAddress("exR.dp" , 		p_ext_delta    );	// checked
	T->SetBranchAddress("exR.p" , 	 	p_ext_mom      );	// checked
	T->SetBranchAddress("exR.ph" , 		p_ext_yptar    );	// checked
	T->SetBranchAddress("exR.th" , 		p_ext_xptar    );	// checked
	T->SetBranchAddress("exR.px" , 		p_ext_px       );	// checked
	T->SetBranchAddress("exR.py" , 		p_ext_py       );	// checked
	T->SetBranchAddress("exR.pz" , 		p_ext_pz       );	// checked
	// Grab stuff for knowing beam trips
	T->SetBranchAddress("RightBCMev.charge_dnew", BCMcharge);	// checked
	T->SetBranchAddress("RightBCMev.current_dnew",BCMcurr  );	// checked
	T->SetBranchAddress("RightBCMev.isrenewed",   BCMrenew );	// checked
	// Grab beam energy from other tree
	E->SetBranchAddress("HALLA_p"     ,&testE  );

	int nEvents = T->GetEntries();
	// -------------------------------------------------------
	// Make the outfile
	char outfilename[200];
	sprintf(outfilename,"/volatile/halla/triton/eep_Rootfiles/skimmed/%s_%s/skim_%d.root",argv[2],argv[3],runNumber);
	TFile * outfile = new TFile(outfilename,"RECREATE");
	TTree * outtree = new TTree("sk","skimmed tree");


	double 	sk_e_cer		,
		sk_e_prl1		,	sk_e_prl2		,
		sk_e_ntrk		,	sk_p_ntrk		,
		sk_e_ytar		,	sk_p_ytar		,
		sk_e_yptar		,	sk_p_yptar		,
		sk_gold_e_yptar		,	sk_gold_p_yptar		,
		sk_e_xptar		,	sk_p_xptar		,
		sk_gold_e_xptar		,	sk_gold_p_xptar		,
		sk_e_delta		,	sk_p_delta		,
		sk_gold_e_delta		,	sk_gold_p_delta		,
		sk_e_mom		,	sk_p_mom		,
		sk_gold_e_mom		,	sk_gold_p_mom		,
		sk_e_beta		,	sk_p_beta		,
		sk_gold_e_beta		,	sk_gold_p_beta		,
		sk_e_x			,	sk_p_x			,
		sk_e_y			,	sk_p_y			,
		sk_e_z			,	sk_p_z			,
		sk_e_px			,	sk_e_py			,
		sk_e_pz			,	sk_p_px			,
		sk_p_py			,	sk_p_pz			,
		sk_gold_e_px		,	sk_gold_p_px		,
		sk_gold_e_py		,	sk_gold_p_py		,
		sk_gold_e_pz		,	sk_gold_p_pz		,
		sk_t1			,	sk_t4			,
		sk_tcoinc		,
		sk_Pmiss		,	sk_Emiss		,
		sk_Pmiss_x		,	sk_Pmiss_y		,
		sk_Pmiss_z		,	sk_p_thWe		,
		sk_th_bq		,	sk_ph_bq		,
		sk_th_xq		,	sk_ph_xq		,
		sk_rast_Pmiss		,	sk_rast_Emiss		,
		sk_rast_Pmiss_x		,	sk_rast_Pmiss_y		,
		sk_rast_Pmiss_z		,	sk_p_rast_thWe		,
		sk_rast_th_bq		,	sk_rast_ph_bq		,
		sk_rast_th_xq		,	sk_rast_ph_xq		,
		sk_exRa_Pmiss		,	sk_exRa_Emiss		,
		sk_exRa_Pmiss_x		,	sk_exRa_Pmiss_y		,
		sk_exRa_Pmiss_z		,	sk_p_exRa_thWe		,
		sk_exRa_th_bq		,	sk_exRa_ph_bq		,
		sk_exRa_th_xq		,	sk_exRa_ph_xq		,
		sk_Q2			,	sk_W2			,
		sk_Nu			,	sk_ph_q			,
		sk_th_q			,	sk_xB			,
		sk_q3m			,	sk_q_x			,
		sk_q_y			,	sk_q_z			,
		sk_e_th			,
		sk_rast_Q2		,	sk_rast_W2		,
		sk_rast_Nu		,	sk_rast_ph_q		,
		sk_rast_th_q		,	sk_rast_xB		,
		sk_rast_q3m		,	sk_rast_q_x		,
		sk_rast_q_y		,	sk_rast_q_z		,		
		sk_e_rast_th		,
		sk_exRa_Q2		,	sk_exRa_W2		,
		sk_exRa_Nu		,	sk_exRa_ph_q		,
		sk_exRa_th_q		,	sk_exRa_xB		,
		sk_exRa_q3m		,	sk_exRa_q_x		,
		sk_exRa_q_y		,	sk_exRa_q_z		,
		sk_e_exRa_th		,
		sk_e_ext_deltaDp	,	sk_e_ext_deltaP		,
		sk_e_ext_deltaTh	,	sk_e_ext_delta		,
		sk_e_ext_mom		,	sk_e_ext_yptar		,
		sk_e_ext_xptar		,	sk_e_ext_px		,
		sk_e_ext_py		,	sk_e_ext_pz		,		
		sk_p_ext_deltaDp	,	sk_p_ext_deltaP		,
		sk_p_ext_deltaTh	,	sk_p_ext_delta		,
		sk_p_ext_mom		,	sk_p_ext_yptar		,
		sk_p_ext_xptar		,	sk_p_ext_px		,
		sk_p_ext_py		,	sk_p_ext_pz		,		
		sk_kin_eThetaC		,	sk_kin_pThetaC		,
		sk_kin_eMomC		,	sk_kin_pMomC		,
		sk_kin_Ebeam		,	sk_kin_Q		,
		sk_BCMcharge		,	sk_BCMcurr		,
		sk_BCMrenew 		;


	// exL.px/py/pz we can get phi
	// and from EKLx.angle we have theta

	outtree -> Branch("L_prl1"		,&sk_e_prl1		,"L_prl1/D"		);
	outtree -> Branch("L_prl2"		,&sk_e_prl2		,"L_prl2/D"		);

	outtree -> Branch("L_cer"		,&sk_e_cer		,"L_cer/D"		);

	outtree -> Branch("L_ntrk"		,&sk_e_ntrk		,"L_ntrk/D"		);
	outtree -> Branch("R_ntrk"		,&sk_p_ntrk		,"R_ntrk/D"		);

	outtree -> Branch("L_ytar"		,&sk_e_ytar		,"L_ytar/D"		);
	outtree -> Branch("R_ytar"		,&sk_p_ytar		,"R_ytar/D"		);

	outtree -> Branch("L_yptar"		,&sk_e_yptar		,"L_yptar/D"		);
	outtree -> Branch("R_yptar"		,&sk_p_yptar		,"R_yptar/D"		);
	outtree -> Branch("L_gold_yptar"	,&sk_gold_e_yptar	,"L_gold_yptar/D"	);
	outtree -> Branch("R_gold_yptar"	,&sk_gold_p_yptar	,"R_gold_yptar/D"	);

	outtree -> Branch("L_xptar"		,&sk_e_xptar		,"L_xptar/D"		);
	outtree -> Branch("R_xptar"		,&sk_p_xptar		,"R_xptar/D"		);
	outtree -> Branch("L_gold_xptar"	,&sk_gold_e_xptar	,"L_gold_xptar/D"	);
	outtree -> Branch("R_gold_xptar"	,&sk_gold_p_xptar	,"R_gold_xptar/D"	);

	outtree -> Branch("L_dp"		,&sk_e_delta		,"L_dp/D"		);
	outtree -> Branch("R_dp"		,&sk_p_delta		,"R_dp/D"		);
	outtree -> Branch("L_gold_dp"		,&sk_gold_e_delta	,"L_gold_dp/D"		);
	outtree -> Branch("R_gold_dp"		,&sk_gold_p_delta	,"R_gold_dp/D"		);

	outtree -> Branch("L_mom"		,&sk_e_mom		,"L_mom/D"		);
	outtree -> Branch("R_mom"		,&sk_p_mom		,"R_mom/D"		);
	outtree -> Branch("L_gold_mom"		,&sk_gold_e_mom		,"L_gold_mom/D"		);
	outtree -> Branch("R_gold_mom"		,&sk_gold_p_mom		,"R_gold_mom/D"		);

	outtree -> Branch("L_beta"		,&sk_e_beta		,"L_beta/D"		);
	outtree -> Branch("R_beta"		,&sk_p_beta		,"R_beta/D"		);
	outtree -> Branch("L_gold_beta"		,&sk_gold_e_beta	,"L_gold_beta/D"	);
	outtree -> Branch("R_gold_beta"		,&sk_gold_p_beta	,"R_gold_beta/D"	);

	outtree -> Branch("L_vx"		,&sk_e_x		,"L_vx/D"		);
	outtree -> Branch("L_vy"		,&sk_e_y		,"L_vy/D"		);
	outtree -> Branch("L_vz"		,&sk_e_z		,"L_vz/D"		);
	outtree -> Branch("R_vx"		,&sk_p_x		,"R_vx/D"		);
	outtree -> Branch("R_vy"		,&sk_p_y		,"R_vy/D"		);
	outtree -> Branch("R_vz"		,&sk_p_z		,"R_vz/D"		);

	outtree -> Branch("L_px"		,&sk_e_px		,"L_px/D"		);
	outtree -> Branch("L_py"		,&sk_e_py		,"L_py/D"		);
	outtree -> Branch("L_pz"		,&sk_e_pz		,"L_pz/D"		);
	outtree -> Branch("R_px"		,&sk_p_px		,"R_px/D"		);
	outtree -> Branch("R_py"		,&sk_p_py		,"R_py/D"		);
	outtree -> Branch("R_pz"		,&sk_p_pz		,"R_pz/D"		);

	outtree -> Branch("L_gold_px"		,&sk_gold_e_px		,"L_gold_px/D"		);
	outtree -> Branch("L_gold_py"		,&sk_gold_e_py		,"L_gold_py/D"		);
	outtree -> Branch("L_gold_pz"		,&sk_gold_e_pz		,"L_gold_pz/D"		);
	outtree -> Branch("R_gold_px"		,&sk_gold_p_px		,"R_gold_px/D"		);
	outtree -> Branch("R_gold_py"		,&sk_gold_p_py		,"R_gold_py/D"		);
	outtree -> Branch("R_gold_pz"		,&sk_gold_p_pz		,"R_gold_pz/D"		);

	outtree -> Branch("t1"			,&sk_t1			,"t1/D"			);
	outtree -> Branch("t4"			,&sk_t4			,"t4/D"			);
	outtree -> Branch("tcoinc"		,&sk_tcoinc		,"tcoinc/D"		);

	outtree -> Branch("Pmiss"		,&sk_Pmiss		,"Pmiss/D"		);
	outtree -> Branch("Emiss"		,&sk_Emiss		,"Emiss/D"		);
	outtree -> Branch("Pmiss_x"		,&sk_Pmiss_x		,"Pmiss_x/D"		);
	outtree -> Branch("Pmiss_y"		,&sk_Pmiss_y		,"Pmiss_y/D"		);
	outtree -> Branch("Pmiss_z"		,&sk_Pmiss_z		,"Pmiss_z/D"		);
	outtree -> Branch("R_thetaWe"		,&sk_p_thWe		,"R_thetaWe/D"		);
	outtree -> Branch("ph_rq"		,&sk_ph_bq		,"ph_rq/D"		);
	outtree -> Branch("th_rq"		,&sk_th_bq		,"th_rq/D"		);
	outtree -> Branch("ph_xq"		,&sk_ph_xq		,"ph_xq/D"		);
	outtree -> Branch("th_xq"		,&sk_th_xq		,"th_xq/D"		);

	outtree -> Branch("rast_Pmiss"		,&sk_rast_Pmiss		,"rast_Pmiss/D"		);
	outtree -> Branch("rast_Emiss"		,&sk_rast_Emiss		,"rast_Emiss/D"		);
	outtree -> Branch("rast_Pmiss_x"	,&sk_rast_Pmiss_x	,"rast_Pmiss_x/D"	);
	outtree -> Branch("rast_Pmiss_y"	,&sk_rast_Pmiss_y	,"rast_Pmiss_y/D"	);
	outtree -> Branch("rast_Pmiss_z"	,&sk_rast_Pmiss_z	,"rast_Pmiss_z/D"	);
	outtree -> Branch("rast_R_thetaWe"	,&sk_p_rast_thWe	,"rast_R_thetaWe/D"	);
	outtree -> Branch("rast_ph_rq"		,&sk_rast_ph_bq		,"rast_ph_rq/D"		);
	outtree -> Branch("rast_th_rq"		,&sk_rast_th_bq		,"rast_th_rq/D"		);
	outtree -> Branch("rast_ph_xq"		,&sk_rast_ph_xq		,"rast_ph_xq/D"		);
	outtree -> Branch("rast_th_xq"		,&sk_rast_th_xq		,"rast_th_xq/D"		);

	outtree -> Branch("exRa_Pmiss"		,&sk_exRa_Pmiss		,"exRa_Pmiss/D"		);
	outtree -> Branch("exRa_Emiss"		,&sk_exRa_Emiss		,"exRa_Emiss/D"		);
	outtree -> Branch("exRa_Pmiss_x"	,&sk_exRa_Pmiss_x	,"exRa_Pmiss_x/D"	);
	outtree -> Branch("exRa_Pmiss_y"	,&sk_exRa_Pmiss_y	,"exRa_Pmiss_y/D"	);
	outtree -> Branch("exRa_Pmiss_z"	,&sk_exRa_Pmiss_z	,"exRa_Pmiss_z/D"	);
	outtree -> Branch("exRa_R_thetaWe"	,&sk_p_exRa_thWe	,"exRa_R_thetaWe/D"	);
	outtree -> Branch("exRa_ph_rq"		,&sk_exRa_ph_bq		,"exRa_ph_rq/D"		);
	outtree -> Branch("exRa_th_rq"		,&sk_exRa_th_bq		,"exRa_th_rq/D"		);
	outtree -> Branch("exRa_ph_xq"		,&sk_exRa_ph_xq		,"exRa_ph_xq/D"		);
	outtree -> Branch("exRa_th_xq"		,&sk_exRa_th_xq		,"exRa_th_xq/D"		);

	outtree -> Branch("Q2"			,&sk_Q2			,"Q2/D"			);
	outtree -> Branch("W2"			,&sk_W2			,"W2/D"			);
	outtree -> Branch("Nu"			,&sk_Nu			,"Nu/D"			);
	outtree -> Branch("ph_q"		,&sk_ph_q		,"ph_q/D"		);
	outtree -> Branch("th_q"		,&sk_th_q		,"th_q/D"		);
	outtree -> Branch("xB"			,&sk_xB			,"xB/D"			);
	outtree -> Branch("q3m"			,&sk_q3m		,"q3m/D"		);
	outtree -> Branch("q_x"			,&sk_q_x		,"q_x/D"		);
	outtree -> Branch("q_y"			,&sk_q_y		,"q_y/D"		);
	outtree -> Branch("q_z"			,&sk_q_z		,"q_z/D"		);
	outtree -> Branch("L_theta"		,&sk_e_th		,"L_theta/D"		);

	outtree -> Branch("rast_Q2"		,&sk_rast_Q2		,"rast_Q2/D"		);
	outtree -> Branch("rast_W2"		,&sk_rast_W2		,"rast_W2/D"		);
	outtree -> Branch("rast_Nu"		,&sk_rast_Nu		,"rast_Nu/D"		);
	outtree -> Branch("rast_ph_q"		,&sk_rast_ph_q		,"rast_ph_q/D"		);
	outtree -> Branch("rast_th_q"		,&sk_rast_th_q		,"rast_th_q/D"		);
	outtree -> Branch("rast_xB"		,&sk_rast_xB		,"rast_xB/D"		);
	outtree -> Branch("rast_q3m"		,&sk_rast_q3m		,"rast_q3m/D"		);
	outtree -> Branch("rast_q_x"		,&sk_rast_q_x		,"rast_q_x/D"		);
	outtree -> Branch("rast_q_y"		,&sk_rast_q_y		,"rast_q_y/D"		);
	outtree -> Branch("rast_q_z"		,&sk_rast_q_z		,"rast_q_z/D"		);
	outtree -> Branch("rast_L_theta"	,&sk_e_rast_th		,"rast_L_theta/D"	);

	outtree -> Branch("exRa_Q2"		,&sk_exRa_Q2		,"exRa_Q2/D"		);
	outtree -> Branch("exRa_W2"		,&sk_exRa_W2		,"exRa_W2/D"		);
	outtree -> Branch("exRa_Nu"		,&sk_exRa_Nu		,"exRa_Nu/D"		);
	outtree -> Branch("exRa_ph_q"		,&sk_exRa_ph_q		,"exRa_ph_q/D"		);
	outtree -> Branch("exRa_th_q"		,&sk_exRa_th_q		,"exRa_th_q/D"		);
	outtree -> Branch("exRa_xB"		,&sk_exRa_xB		,"exRa_xB/D"		);
	outtree -> Branch("exRa_q3m"		,&sk_exRa_q3m		,"exRa_q3m/D"		);
	outtree -> Branch("exRa_q_x"		,&sk_exRa_q_x		,"exRa_q_x/D"		);
	outtree -> Branch("exRa_q_y"		,&sk_exRa_q_y		,"exRa_q_y/D"		);
	outtree -> Branch("exRa_q_z"		,&sk_exRa_q_z		,"exRa_q_z/D"		);
	outtree -> Branch("exRa_L_theta"	,&sk_e_exRa_th		,"exRa_L_theta/D"	);

	outtree -> Branch("L_ext_delta_dp"	,&sk_e_ext_deltaDp	,"L_ext_delta_dp/D"	);
	outtree -> Branch("L_ext_delta_p"	,&sk_e_ext_deltaP	,"L_ext_delta_p/D"	);
	outtree -> Branch("L_ext_delta_yptar"	,&sk_e_ext_deltaTh	,"L_ext_delta_yptar/D"	);
	outtree -> Branch("L_ext_dp"		,&sk_e_ext_delta	,"L_ext_dp/D"		);
	outtree -> Branch("L_ext_mom"		,&sk_e_ext_mom		,"L_ext_mom/D"		);
	outtree -> Branch("L_ext_yptar"		,&sk_e_ext_yptar	,"L_ext_yptar/D"	);
	outtree -> Branch("L_ext_xptar"		,&sk_e_ext_xptar	,"L_ext_xptar/D"	);
	outtree -> Branch("L_ext_px"		,&sk_e_ext_px		,"L_ext_px/D"		);
	outtree -> Branch("L_ext_py"		,&sk_e_ext_py		,"L_ext_py/D"		);
	outtree -> Branch("L_ext_pz"		,&sk_e_ext_pz		,"L_ext_pz/D"		);	

	outtree -> Branch("R_ext_delta_dp"	,&sk_p_ext_deltaDp	,"R_ext_delta_dp/D"	);
	outtree -> Branch("R_ext_delta_p"	,&sk_p_ext_deltaP	,"R_ext_delta_p/D"	);
	outtree -> Branch("R_ext_delta_yptar"	,&sk_p_ext_deltaTh	,"R_ext_delta_yptar/D"	);
	outtree -> Branch("R_ext_dp"		,&sk_p_ext_delta	,"R_ext_dp/D"		);
	outtree -> Branch("R_ext_mom"		,&sk_p_ext_mom		,"R_ext_mom/D"		);
	outtree -> Branch("R_ext_yptar"		,&sk_p_ext_yptar	,"R_ext_yptar/D"	);
	outtree -> Branch("R_ext_xptar"		,&sk_p_ext_xptar	,"R_ext_xptar/D"	);
	outtree -> Branch("R_ext_px"		,&sk_p_ext_px		,"R_ext_px/D"		);
	outtree -> Branch("R_ext_py"		,&sk_p_ext_py		,"R_ext_py/D"		);
	outtree -> Branch("R_ext_pz"		,&sk_p_ext_pz		,"R_ext_pz/D"		);	

	outtree -> Branch("Kin_L_thetaC"	,&sk_kin_eThetaC	,"Kin_L_thetaC/D"	);
	outtree -> Branch("Kin_R_thetaC"	,&sk_kin_pThetaC	,"Kin_R_thetaC/D"	);
	outtree -> Branch("Kin_L_momC"		,&sk_kin_eMomC		,"Kin_L_momC/D"		);
	outtree -> Branch("Kin_R_momC"		,&sk_kin_pMomC		,"Kin_R_momC/D"		);
	outtree -> Branch("Kin_Ebeam"		,&sk_kin_Ebeam		,"Kin_Ebeam/D"		);
	outtree -> Branch("Kin_Q"		,&sk_kin_Q		,"Kin_Q/D"		);

	outtree -> Branch("BCM_curr"		,&sk_BCMcurr		,"BCM_curr/D"		);
	outtree -> Branch("BCM_charge"		,&sk_BCMcharge		,"BCM_charge/D"		);
	outtree -> Branch("BCM_isrenew"		,&sk_BCMrenew		,"BCM_isrenew/D"	);


	// Loop over events
	cout << "Looping for beam energy \n";
	double Ebeam = 0.;
	int nVal = 0;
	for (int event=0 ; event < E->GetEntries() ; event++){
		E->GetEvent(event);
		if( testE == 4309.12 ) continue;
		if( testE == -1e+32 ) continue;
		Ebeam += testE;
		nVal++;
	}
	Ebeam /= (double) nVal;
	cout << "\tAverage beam energy from BCM with orbit corrections: " << Ebeam << "\n";
	Ebeam *= 1.0025;
	Ebeam /= 1000.;
	cout << "\t**Beam energy that we will use with Doug's black-box correction: "  << Ebeam << "**\n";
	// This is now stored in kinematics, which is a 5-len vector
	kinematics[4] = Ebeam;


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Fill the tree
	cout << nEvents << "\n";
	cout << "Beginning to loop over " << nEvents << " events.\n";
	for (int event=0 ; event < nEvents ; event++)
	{
		if (event%10000 ==0)
			cerr << "Working on event " << event << " out of " << nEvents << "\n";
		/////////////////////////////////////////////////////////////////////////////////////////////
		// Memory set the arrays to avoid weird shit
		memset(	e_cer			,-999		,sizeof(e_cer)		);
		memset(	e_prl1			,-999		,sizeof(e_prl1)		);
		memset(	e_prl2			,-999		,sizeof(e_prl2)		);

		memset( e_ntrk			,-999		,sizeof(e_ntrk) 	);
		memset( p_ntrk			,-999		,sizeof(p_ntrk) 	);

		memset( e_ytar			,-999		,sizeof(e_ytar) 	);
		memset( p_ytar			,-999		,sizeof(p_ytar) 	);

		memset( e_yptar			,-999		,sizeof(e_yptar) 	);
		memset( p_yptar			,-999		,sizeof(p_yptar) 	);
		memset( gold_e_yptar		,-999		,sizeof(gold_e_yptar) 	);
		memset( gold_p_yptar		,-999		,sizeof(gold_p_yptar) 	);
		memset( e_xptar			,-999		,sizeof(e_xptar) 	);
		memset( p_xptar			,-999		,sizeof(p_xptar) 	);
		memset( gold_e_xptar		,-999		,sizeof(gold_e_xptar) 	);
		memset( gold_p_xptar		,-999		,sizeof(gold_p_xptar) 	);

		memset( e_delta			,-999		,sizeof(e_delta) 	);
		memset( p_delta			,-999		,sizeof(p_delta) 	);
		memset( gold_e_delta		,-999		,sizeof(gold_e_delta) 	);
		memset( gold_p_delta		,-999		,sizeof(gold_p_delta) 	);

		memset( e_mom			,-999		,sizeof(e_mom)		);
		memset( p_mom			,-999		,sizeof(p_mom)		);
		memset( gold_e_mom		,-999		,sizeof(gold_e_mom)	);
		memset( gold_p_mom		,-999		,sizeof(gold_p_mom)	);

		memset( e_beta			,-999		,sizeof(e_beta)	);
		memset( p_beta			,-999		,sizeof(p_beta)	);
		memset( gold_e_beta		,-999		,sizeof(gold_e_beta)	);
		memset( gold_p_beta		,-999		,sizeof(gold_p_beta)	);

		memset( e_x			,-999		,sizeof(e_x)		);
		memset( e_y			,-999		,sizeof(e_y)		);
		memset( e_z			,-999		,sizeof(e_z)		);
		memset( p_x			,-999		,sizeof(p_x)		);
		memset( p_y			,-999		,sizeof(p_y)		);
		memset( p_z			,-999		,sizeof(p_z)		);

		memset( e_px			,-999		,sizeof(e_px)		);
		memset( e_py			,-999		,sizeof(e_py)		);
		memset( e_pz			,-999		,sizeof(e_pz)		);
		memset( p_px			,-999		,sizeof(p_px)		);
		memset( p_py			,-999		,sizeof(p_py)		);
		memset( p_pz			,-999		,sizeof(p_pz)		);

		memset( gold_e_px		,-999		,sizeof(gold_e_px)	);
		memset( gold_e_py		,-999		,sizeof(gold_e_py)	);
		memset( gold_e_pz		,-999		,sizeof(gold_e_pz)	);
		memset( gold_p_px		,-999		,sizeof(gold_p_px)	);
		memset( gold_p_py		,-999		,sizeof(gold_p_py)	);
		memset( gold_p_pz		,-999		,sizeof(gold_p_pz)	);

		memset(	t1			,-999		,sizeof(t1)		);
		memset(	t4			,-999		,sizeof(t4)		);

		memset(	tcoinc			,-999		,sizeof(tcoinc)		);

		memset(	Pmiss			,-999		,sizeof(Pmiss)		);
		memset(	Emiss			,-999		,sizeof(Emiss)		);
		memset(	Pmiss_x			,-999		,sizeof(Pmiss_x)	);
		memset(	Pmiss_y			,-999		,sizeof(Pmiss_y)	);
		memset(	Pmiss_z			,-999		,sizeof(Pmiss_z)	);
		memset(	p_thWe			,-999		,sizeof(p_thWe)		);
		memset( ph_bq			,-999		,sizeof(ph_bq)		);
		memset( th_bq			,-999		,sizeof(th_bq)		);
		memset( ph_xq			,-999		,sizeof(ph_xq)		);
		memset( th_xq			,-999		,sizeof(th_xq)		);

		memset(	rast_Pmiss		,-999		,sizeof(rast_Pmiss)	);
		memset(	rast_Emiss		,-999		,sizeof(rast_Emiss)	);
		memset(	rast_Pmiss_x		,-999		,sizeof(rast_Pmiss_x)	);
		memset(	rast_Pmiss_y		,-999		,sizeof(rast_Pmiss_y)	);
		memset(	rast_Pmiss_z		,-999		,sizeof(rast_Pmiss_z)	);
		memset(	p_rast_thWe		,-999		,sizeof(p_rast_thWe)	);
		memset( rast_ph_bq		,-999		,sizeof(rast_ph_bq)	);
		memset( rast_th_bq		,-999		,sizeof(rast_th_bq)	);
		memset( rast_ph_xq		,-999		,sizeof(rast_ph_xq)	);
		memset( rast_th_xq		,-999		,sizeof(rast_th_xq)	);

		memset( exRa_ph_bq		,-999		,sizeof(exRa_ph_bq)	);
		memset( exRa_th_bq		,-999		,sizeof(exRa_th_bq)	);
		memset( exRa_ph_xq		,-999		,sizeof(exRa_ph_xq)	);
		memset( exRa_th_xq		,-999		,sizeof(exRa_th_xq)	);
		memset(	exRa_Pmiss		,-999		,sizeof(exRa_Pmiss)	);
		memset(	exRa_Emiss		,-999		,sizeof(exRa_Emiss)	);
		memset(	exRa_Pmiss_x		,-999		,sizeof(exRa_Pmiss_x)	);
		memset(	exRa_Pmiss_y		,-999		,sizeof(exRa_Pmiss_y)	);
		memset(	exRa_Pmiss_z		,-999		,sizeof(exRa_Pmiss_z)	);
		memset(	p_exRa_thWe		,-999		,sizeof(p_exRa_thWe)	);
		memset( exRa_ph_bq		,-999		,sizeof(exRa_ph_bq)	);
		memset( exRa_th_bq		,-999		,sizeof(exRa_th_bq)	);
		memset( exRa_ph_xq		,-999		,sizeof(exRa_ph_xq)	);
		memset( exRa_th_xq		,-999		,sizeof(exRa_th_xq)	);

		memset( Q2			,-999		,sizeof(Q2)		);
		memset( W2			,-999		,sizeof(W2)		);
		memset( Nu			,-999		,sizeof(Nu)		);
		memset( ph_q			,-999		,sizeof(ph_q)		);
		memset( th_q			,-999		,sizeof(th_q)		);
		memset( xB			,-999		,sizeof(xB)		);
		memset( q3m			,-999		,sizeof(q3m)		);
		memset( q_x			,-999		,sizeof(q_x)		);
		memset( q_y			,-999		,sizeof(q_y)		);
		memset( q_z			,-999		,sizeof(q_z)		);
		memset( e_th			,-999		,sizeof(e_th)		);

		memset( rast_Q2			,-999		,sizeof(rast_Q2)	);
		memset( rast_W2			,-999		,sizeof(rast_W2)	);
		memset( rast_Nu			,-999		,sizeof(rast_Nu)	);
		memset( rast_ph_q		,-999		,sizeof(rast_ph_q)	);
		memset( rast_th_q		,-999		,sizeof(rast_th_q)	);
		memset( rast_xB			,-999		,sizeof(rast_xB)	);
		memset( rast_q3m		,-999		,sizeof(rast_q3m)	);
		memset( rast_q_x		,-999		,sizeof(rast_q_x)	);
		memset( rast_q_y		,-999		,sizeof(rast_q_y)	);
		memset( rast_q_z		,-999		,sizeof(rast_q_z)	);
		memset( e_rast_th		,-999		,sizeof(e_rast_th)	);

		memset( exRa_Q2			,-999		,sizeof(exRa_Q2)	);
		memset( exRa_W2			,-999		,sizeof(exRa_W2)	);
		memset( exRa_Nu			,-999		,sizeof(exRa_Nu)	);
		memset( exRa_ph_q		,-999		,sizeof(exRa_ph_q)	);
		memset( exRa_th_q		,-999		,sizeof(exRa_th_q)	);
		memset( exRa_xB			,-999		,sizeof(exRa_xB)	);
		memset( exRa_q3m		,-999		,sizeof(exRa_q3m)	);
		memset( exRa_q_x		,-999		,sizeof(exRa_q_x)	);
		memset( exRa_q_y		,-999		,sizeof(exRa_q_y)	);
		memset( exRa_q_z		,-999		,sizeof(exRa_q_z)	);
		memset( e_exRa_th		,-999		,sizeof(e_exRa_th)	);

		memset(	e_ext_deltaDp		,-999		,sizeof(e_ext_deltaDp)	);
		memset(	e_ext_deltaP		,-999		,sizeof(e_ext_deltaP)	);
		memset(	e_ext_deltaTh		,-999		,sizeof(e_ext_deltaTh)	);

		memset(	e_ext_delta		,-999		,sizeof(e_ext_delta)	);
		memset( e_ext_mom		,-999		,sizeof(e_ext_mom)	);
		memset( e_ext_yptar		,-999		,sizeof(e_ext_yptar)	);
		memset( e_ext_xptar		,-999		,sizeof(e_ext_xptar)	);
		memset( e_ext_px		,-999		,sizeof(e_ext_px)	);
		memset( e_ext_py		,-999		,sizeof(e_ext_py)	);
		memset( e_ext_pz		,-999		,sizeof(e_ext_pz)	);

		memset(	p_ext_deltaDp		,-999		,sizeof(p_ext_deltaDp)	);
		memset(	p_ext_deltaP		,-999		,sizeof(p_ext_deltaP)	);
		memset(	p_ext_deltaTh		,-999		,sizeof(p_ext_deltaTh)	);

		memset(	p_ext_delta		,-999		,sizeof(p_ext_delta)	);
		memset( p_ext_mom		,-999		,sizeof(p_ext_mom)	);
		memset( p_ext_yptar		,-999		,sizeof(p_ext_yptar)	);
		memset( p_ext_xptar		,-999		,sizeof(p_ext_xptar)	);
		memset( p_ext_px		,-999		,sizeof(p_ext_px)	);
		memset( p_ext_py		,-999		,sizeof(p_ext_py)	);
		memset( p_ext_pz		,-999		,sizeof(p_ext_pz)	);

		memset( BCMcharge		,-999		,sizeof(BCMcharge)	);
		memset( BCMcurr			,-999		,sizeof(BCMcurr)	);
		memset( BCMrenew		,-999		,sizeof(BCMrenew)	);

		/////////////////////////////////////////////////////////////////////////////////////////////
		// Get the event
		T->GetEvent(event);

		/////////////////////////////////////////////////////////////////////////////////////////////
		// Get the charge 
		if( BCMrenew[0] > 0.9){
			if( isnan( (float) BCMcharge[0]) ) continue;
			testCharge += BCMcharge[0];
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		// Do our simple skimming
		// Note that this is different from the "normal" skimmer in that here we
		// require a track on the left OR a track on the right rather than
		// requiring a track on the left AND a track on the right
		if (!(p_ntrk[0]==1)&&!(e_ntrk[0]==1)) continue;

		////////////////////////////////////////////////////////////////////////////////////////////
		// Fill our tree
		//

		// fill kinematic entries
		sk_kin_eThetaC	=	e_theta_c	;
		sk_kin_pThetaC	=	p_theta_c	;
		sk_kin_eMomC	=	e_mom_c		;
		sk_kin_pMomC	=	p_mom_c		;
		sk_kin_Ebeam 	=	Ebeam		;
		sk_kin_Q	=	totalQ_mC	;	


		sk_e_cer		= e_cer		[0];
		sk_e_prl1		= e_prl1	[0];
		sk_e_prl2		= e_prl2	[0];
		sk_e_ntrk		= e_ntrk	[0];
		sk_p_ntrk		= p_ntrk	[0];
		sk_e_ytar		= e_ytar	[0];
		sk_p_ytar		= p_ytar	[0];
		sk_e_yptar		= e_yptar	[0];
		sk_p_yptar		= p_yptar	[0];
		sk_gold_e_yptar		= gold_e_yptar	[0];
		sk_gold_p_yptar		= gold_p_yptar	[0];
		sk_e_xptar		= e_xptar	[0];
		sk_p_xptar		= p_xptar	[0];
		sk_gold_e_xptar		= gold_e_xptar	[0];
		sk_gold_p_xptar		= gold_p_xptar	[0];
		sk_e_delta		= e_delta	[0];
		sk_p_delta		= p_delta	[0];
		sk_gold_e_delta		= gold_e_delta	[0];
		sk_gold_p_delta		= gold_p_delta	[0];
		sk_e_mom		= e_mom		[0];
		sk_p_mom		= p_mom		[0];
		sk_gold_e_mom		= gold_e_mom	[0];
		sk_gold_p_mom		= gold_p_mom	[0];
		sk_e_beta		= e_beta	[0];
		sk_p_beta		= p_beta	[0];
		sk_gold_e_beta		= gold_e_beta	[0];
		sk_gold_p_beta		= gold_p_beta	[0];
		sk_e_x			= e_x		[0];
		sk_p_x			= p_x		[0];
		sk_e_y			= e_y		[0];
		sk_p_y			= p_y		[0];
		sk_e_z			= e_z		[0];
		sk_p_z			= p_z		[0];
		sk_e_px			= e_px		[0];
		sk_p_px			= p_px		[0];
		sk_e_py			= e_py		[0];
		sk_p_py			= p_py		[0];
		sk_e_pz			= e_pz		[0];
		sk_p_pz			= p_pz		[0];
		sk_gold_e_px		= gold_e_px	[0];
		sk_gold_p_px		= gold_p_px	[0];
		sk_gold_e_py		= gold_e_py	[0];
		sk_gold_p_py		= gold_p_py	[0];
		sk_gold_e_pz		= gold_e_pz	[0];
		sk_gold_p_pz		= gold_p_pz	[0];
		sk_t1			= t1		[0];
		sk_t4			= t4		[0];
		sk_tcoinc		= tcoinc	[0];
		sk_Pmiss		= Pmiss		[0];
		sk_Emiss		= Emiss		[0];
		sk_Pmiss_x		= Pmiss_x	[0];
		sk_Pmiss_y		= Pmiss_y	[0];
		sk_Pmiss_z		= Pmiss_z	[0];
		sk_p_thWe		= p_thWe	[0];
		sk_ph_bq		= ph_bq		[0];
		sk_th_bq		= th_bq		[0];
		sk_ph_xq		= ph_xq		[0];
		sk_th_xq		= th_xq		[0];
		sk_rast_Pmiss		= rast_Pmiss	[0];
		sk_rast_Emiss		= rast_Emiss	[0];
		sk_rast_Pmiss_x		= rast_Pmiss_x	[0];
		sk_rast_Pmiss_y		= rast_Pmiss_y	[0];
		sk_rast_Pmiss_z		= rast_Pmiss_z	[0];
		sk_p_rast_thWe		= p_rast_thWe	[0];
		sk_rast_ph_bq		= rast_ph_bq	[0];
		sk_rast_th_bq		= rast_th_bq	[0];
		sk_rast_ph_xq		= rast_ph_xq	[0];
		sk_rast_th_xq		= rast_th_xq	[0];
		sk_exRa_Pmiss		= exRa_Pmiss	[0];
		sk_exRa_Emiss		= exRa_Emiss	[0];
		sk_exRa_Pmiss_x		= exRa_Pmiss_x	[0];
		sk_exRa_Pmiss_y		= exRa_Pmiss_y	[0];
		sk_exRa_Pmiss_z		= exRa_Pmiss_z	[0];
		sk_p_exRa_thWe		= p_exRa_thWe	[0];
		sk_exRa_ph_bq		= exRa_ph_bq	[0];
		sk_exRa_th_bq		= exRa_th_bq	[0];
		sk_exRa_ph_xq		= exRa_ph_xq	[0];
		sk_exRa_th_xq		= exRa_th_xq	[0];
		sk_Q2			= Q2		[0];
		sk_W2			= W2		[0];
		sk_Nu			= Nu		[0];
		sk_ph_q			= ph_q		[0];
		sk_th_q			= th_q		[0];
		sk_xB			= xB		[0];
		sk_q3m			= q3m		[0];
		sk_q_x			= q_x		[0];
		sk_q_y			= q_y		[0];
		sk_q_z			= q_z		[0];
		sk_e_th			= e_th		[0];
		sk_rast_Q2		= rast_Q2	[0];
		sk_rast_W2		= rast_W2	[0];
		sk_rast_Nu		= rast_Nu	[0];
		sk_rast_ph_q		= rast_ph_q	[0];
		sk_rast_th_q		= rast_th_q	[0];
		sk_rast_xB		= rast_xB	[0];
		sk_rast_q3m		= rast_q3m	[0];
		sk_rast_q_x		= rast_q_x	[0];
		sk_rast_q_y		= rast_q_y	[0];
		sk_rast_q_z		= rast_q_z	[0];
		sk_e_rast_th		= e_rast_th	[0];
		sk_exRa_Q2		= exRa_Q2	[0];
		sk_exRa_W2		= exRa_W2	[0];
		sk_exRa_Nu		= exRa_Nu	[0];
		sk_exRa_ph_q		= exRa_ph_q	[0];
		sk_exRa_th_q		= exRa_th_q	[0];
		sk_exRa_xB		= exRa_xB	[0];
		sk_exRa_q3m		= exRa_q3m	[0];
		sk_exRa_q_x		= exRa_q_x	[0];
		sk_exRa_q_y		= exRa_q_y	[0];
		sk_exRa_q_z		= exRa_q_z	[0];
		sk_e_exRa_th		= e_exRa_th	[0];
		sk_e_ext_deltaDp	= e_ext_deltaDp	[0];
		sk_e_ext_deltaP		= e_ext_deltaP	[0];
		sk_e_ext_deltaTh	= e_ext_deltaTh	[0];
		sk_e_ext_delta		= e_ext_delta	[0];
		sk_e_ext_mom		= e_ext_mom	[0];
		sk_e_ext_yptar		= e_ext_yptar	[0];
		sk_e_ext_xptar		= e_ext_xptar	[0];
		sk_e_ext_px		= e_ext_px	[0];
		sk_e_ext_py		= e_ext_py	[0];
		sk_e_ext_pz		= e_ext_pz	[0];
		sk_p_ext_deltaDp	= p_ext_deltaDp	[0];
		sk_p_ext_deltaP		= p_ext_deltaP	[0];
		sk_p_ext_deltaTh	= p_ext_deltaTh	[0];
		sk_p_ext_delta		= p_ext_delta	[0];
		sk_p_ext_mom		= p_ext_mom	[0];
		sk_p_ext_yptar		= p_ext_yptar	[0];
		sk_p_ext_xptar		= p_ext_xptar	[0];
		sk_p_ext_px		= p_ext_px	[0];
		sk_p_ext_py		= p_ext_py	[0];
		sk_p_ext_pz		= p_ext_pz	[0];
		sk_BCMcurr		= BCMcurr	[0];
		sk_BCMcharge		= BCMcharge	[0];
		sk_BCMrenew		= BCMrenew	[0];

		// Calculated variables:
		// I TOOK ALL THESE OUT OF SKIMMER
		// YOU SHOULD DO THESE IN OWN ANALYSIS
		// SCRIPT
		/*
		   sk_eEoverP = (sk_prl1 + sk_prl2)/1000
		   sk_eEoverP = (L_prl1[0] + L_prl2[0])/1000./(sk_e_mom);

		   sk_e_theta = TMath::ACos( (TMath::Cos(e_theta_c) - e_yptar[0]*TMath::Sin(e_theta_c)) / TMath::Sqrt(1. + e_yptar[0]*e_yptar[0] + e_xptar[0]*e_xptar[0]));
		   sk_e_phi   = TMath::ATan2( TMath::Sin(e_theta_c) + e_yptar[0] * TMath::Cos(e_theta_c),e_xptar[0]) - TMath::Pi()/2.;
		   sk_p_theta = TMath::ACos( (TMath::Cos(p_theta_c) + p_yptar[0]*TMath::Sin(p_theta_c)) / TMath::Sqrt(1. + p_yptar[0]*p_yptar[0] + p_xptar[0]*p_xptar[0]));
		   sk_p_phi   = TMath::ATan2(-TMath::Sin(p_theta_c) + p_yptar[0] * TMath::Cos(p_theta_c),p_xptar[0]) + 1.5*TMath::Pi();

		// Calculate the missing mass
		// Went back to calculating pp and q because 
		// these values include shifts to the central momenta
		TLorentzVector q(-sk_e_mom*sin(sk_e_theta)*cos(sk_e_phi),
		-sk_e_mom*sin(sk_e_theta)*sin(sk_e_phi),
		Ebeam - sk_e_mom*cos(sk_e_theta),
		Ebeam - sk_e_mom);
		TLorentzVector pp(sk_p_mom*sin(sk_p_theta)*cos(sk_p_phi),
		sk_p_mom*sin(sk_p_theta)*sin(sk_p_phi),
		sk_p_mom*cos(sk_p_theta),
		TMath::Sqrt(sk_p_mom*sk_p_mom + mP*mP));

		TLorentzVector A(0.,0.,0.,mTarg);
		TVector3 V3pMiss = q.Vect()  - pp.Vect();

		sk_Q2    = -q.Mag2();
		sk_Nu    = Ebeam - sk_e_mom;
		sk_Pm    = V3pMiss.Mag();
		sk_Pmx   = V3pMiss.X();
		sk_Pmy   = V3pMiss.Y();
		sk_Pmz   = V3pMiss.Z();
		sk_thrq  = V3pMiss.Angle(q.Vect());
		sk_xB    = sk_Q2 / (2.*mP*sk_Nu);
		sk_Em_Larry = sk_Nu - TMath::Sqrt(sk_p_mom*sk_p_mom + mP*mP) + mP; 	   
		sk_Em    = sk_Nu - (TMath::Sqrt(sk_p_mom*sk_p_mom + mP*mP) - mP) - (sqrt(sk_Pm*sk_Pm + mD*mD) - mD); 
		*/

		outtree -> Fill();
	}

	//////////////////////////////////////////////////////////////////////////////////////
	// Ending
	cout << "Accumulated charge: " << totalQ_mC << " mC" << endl; 
	cout << "Charge in this file: " << testCharge/1000. << " mC" << endl;
	outfile->cd   ();
	outtree->Write();
	totalCharge.Write("totalCharge");
	kinematics.Write("kinematics");
	outfile->Close();
	delete T;

	return 0;
}
// ====================================================================================================
double get_kin(const char * kin, double &e_th, double &e_p, double &p_th, double &p_p)
{
	if (strcmp(kin,"fast")==0)
	{
		// THESE NUMBERS UPDATED BASED ON EPICS DATA
		e_th = 20.892 * TMath::DegToRad();
		e_p  = 3.54332; // GeV/c
		p_th = 48.8094 * TMath::DegToRad();
		p_p  = 1.4805; // GeV/c
		cout << "Assuming fast kinematics\n(LHRS: p="<<e_p<<"GeV, Th="<< e_th<<"rad. RHRS: p="
			<<p_p<<"GeV, Th="<<p_th<<"rad)"<<endl;
	}
	else if (strcmp(kin,"fast2")==0)
	{
		// THESE NUMBERS UPDATED BASED ON EPICS DATA
		e_th = 20.8798 * TMath::DegToRad();
		e_p  = 3.5433; // GeV/c
		p_th = 48.799 * TMath::DegToRad();
		p_p  = 1.48051; // GeV/c
		cout << "Assuming fast2 kinematics\n(LHRS: p="<<e_p<<"GeV, Th="<< e_th<<"rad. RHRS: p="
			<<p_p<<"GeV, Th="<<p_th<<"rad)"<<endl;
	}
	else if (strcmp(kin,"slow")==0)
	{
		// THESE ANGLES HAVE BEEN UPDATED BASED ON THE SURVEY DONE
		// 	still need to ''move rear jacks along tangent???
		// The momenta updated based on EPICs data
		e_th = 20.892 * TMath::DegToRad(); //20.88 * TMath::DegToRad();
		e_p  = 3.54332; // GeV/c
		p_th = 58.490 * TMath::DegToRad(); //58.50 * TMath::DegToRad();
		p_p  = 1.246; // GeV/c
		cout << "Assuming slow kinematics\n(LHRS: p="<<e_p<<"GeV, Th="<< e_th<<"rad. RHRS: p="
			<<p_p<<"GeV, Th="<<p_th<<"rad)"<<endl;
	}
	else if (strcmp(kin,"mid")==0)
	{
		// THESE NUMBERS UPDATED BASED ON EPICS DATA
		e_th = 17.8018* TMath::DegToRad(); //17.8003 * TMath::DegToRad();
		e_p  = 3.54334;//3.54334; // GeV/c
		p_th = -48.82* TMath::DegToRad(); //48.8094 * TMath::DegToRad(); 
		p_p  = 1.4805; //1.4805; // GeV/c

		// Add corrections found in Efrain toy-analysis for mid
		//e_th += 0.00067;
		//p_th += 0.000445557;
		//e_p   += -0.000268191;
		//p_p   += -0.000390631; 	


		cout << "Assuming mid kinematics\n(LHRS: p="<<e_p<<"GeV, Th="<< e_th<<"rad. RHRS: p="
			<< p_p<<"GeV, Th="<<p_th<<"rad)"<<endl;
	}
	else if (strcmp(kin,"mid2")==0)
	{
		// THESE ANGLES UPDATED BASED ON LVDT DATA
		e_th = 17.8053* TMath::DegToRad(); 
		e_p  = 3.54334;// GeV/c
		p_th = -48.8188* TMath::DegToRad(); 
		p_p  = 1.4805; // GeV/c

		cout << "Assuming mid2 kinematics\n(LHRS: p="<<e_p<<"GeV, Th="<< e_th<<"rad. RHRS: p="
			<< p_p<<"GeV, Th="<<p_th<<"rad)"<<endl;
	}
	else
	{
		cerr << "You asked for the following kinematics: " << kin << "\n";
		cerr << "  This is not a valid setting. Bailing out.\n";
		exit(-2);
	}
}
// ====================================================================================================
double get_mTarg(const char * filename){
	TFile * firstFile = new TFile(filename);
	if (firstFile->IsZombie())
	{
		cerr << "The first rootfile provided cannot be opened. Bailing out.\n";
		return -3;
	}
	TTree * targTree = (TTree*)firstFile->Get("E");
	Double_t targ_pos=-1.E32;
	targTree->SetBranchAddress("haBDSPOS.VAL",&targ_pos);
	int counter=0;
	while (targ_pos < 1E5)
	{
		targTree->GetEvent(counter);
		counter +=1;
	}
	firstFile->Close();

	double fabsrange=1000.;

	if (fabs(targ_pos - 33106235.)<fabsrange)	return mH3;
	if (fabs(targ_pos - 21875520.)<fabsrange)	return mHe3;
	if (fabs(targ_pos - 29367355.)<fabsrange)	return mD;
	if (fabs(targ_pos - 13108119.)<fabsrange)	return mC;
	if (fabs(targ_pos - 25610043.)<fabsrange)	return mP;
	return mAl;
}
// ====================================================================================================
double Calculate_Charge_from_Scaler(const char * filename){

	TFile * infile = new TFile(filename    );
	TTree * TS     = (TTree*)infile->Get("TSRight");

	Double_t Dnew_count, D3_count, CL_count;

	TS -> SetBranchAddress("Rightdnew"  ,&Dnew_count);
	TS -> SetBranchAddress("Rightd3"    ,&D3_count  );
	TS -> SetBranchAddress("RightLclock",&CL_count  );

	Double_t Time, BCM_dnew_count, BCM_d3_count;
	Double_t CL_count_previous  = 0.0;
	Double_t Dnew_count_previous= 0.0;
	Double_t D3_count_previous  = 0.0;
	Double_t CHARGE   = -1.0;
	Double_t CHARGE_d3= -1.0;

	const Double_t fast_clock = 103700.0  ; // Hz - fast clock
	const Double_t slow_clock =   1024.0  ; // Hz - slow clock
	const Double_t clock_rate = fast_clock;

	// Coefficients from BCM calibration
	const Double_t coef1      = 0.0003264;
	const Double_t coef2      = 0.1055   ;
	const Double_t coef3      = 0.0001063;
	const Double_t coef4      = 0.1985   ;

	int Ev_No = TS->GetEntries();

	for(Int_t index=0;index<Ev_No;index++){
		TS->GetEntry(index);
		if(Dnew_count>Dnew_count_previous){
			Time              += (CL_count-CL_count_previous)/clock_rate;
			BCM_dnew_count    += Dnew_count-Dnew_count_previous         ;
			BCM_d3_count      += D3_count-D3_count_previous             ;
			CL_count_previous  = CL_count                               ;
			Dnew_count_previous= Dnew_count                             ;
			D3_count_previous  = D3_count                               ;
		}
	}

	CHARGE   = coef1*BCM_dnew_count + coef2*Time*1e-6;
	CHARGE_d3= coef3*BCM_d3_count   + coef4*Time*1e-6;

	infile->Close();

	if(CHARGE>0) return CHARGE/1000.;

	cerr << "Calculated charge is negative. What's wrong with the world?! Bailing out\n";
	exit(1);
}

// ====================================================================================================
double Calculate_Time(const char * filename){

	TFile * infile = new TFile(filename    );
	TTree * TS     = (TTree*)infile->Get("TSRight");

	Double_t Dnew_count, CL_count;

	TS -> SetBranchAddress("RightLclock",&CL_count  );

	Double_t Time, BCM_dnew_count;
	Double_t CL_count_previous  = 0.0;
	Double_t Dnew_count_previous= 0.0;

	const Double_t fast_clock = 103700.0  ; // Hz - fast clock
	const Double_t slow_clock =   1024.0  ; // Hz - slow clock
	const Double_t clock_rate = fast_clock;

	int Ev_No = TS->GetEntries();

	for(Int_t index=0;index<Ev_No;index++){
		TS->GetEntry(index);
		if(Dnew_count>Dnew_count_previous){
			Time += (CL_count-CL_count_previous)/clock_rate;
		}
	}

	infile->Close();
	return Time;
}

