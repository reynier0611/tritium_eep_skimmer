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

// ================================================================================================
int main(int argc, char ** argv)
{
	if ( argc<2 )
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "   combiner [list of runs...]\n\n";
		return -1;
	}

	TVectorT<double> kinematics(5);
	TVectorT<double> totalCharge(1);
	totalCharge[0] = 0.;

	// Establish the output file
	TFile * outfile = new TFile("combine_out.root","RECREATE");
	TTree * outtree = new TTree("sk","skimmed tree");
	
	ofstream outtxt ("eventCounts.txt");

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

	int totEv = 0; // number of events for all files
	int totEv_0to50 = 0;
	int totEv_50to100 = 0;
	int totEv_100to150 = 0;
	int totEv_150to200 = 0;
	int totEv_200to250 = 0;
	int totEv_250to300 = 0;
	int totEv_300to350 = 0;
	int totEv_350to400 = 0;
	int totEv_400to450 = 0;
	int totEv_450to600 = 0;
	int totEv_above600 = 0;
	// Loop over input files
	for (int i=1 ; i<argc ; i++)
	{
		//cout << "Working on file " << argv[i] << "\n";
		TFile * thisFile = new TFile(argv[i]);

		if (i==1) // Copy over the kinematics
		{
			TVectorT<double> * kin = (TVectorT<double>*)thisFile->Get("kinematics");
			kinematics[0] = (*kin)[0];
			kinematics[1] = (*kin)[1];
			kinematics[2] = (*kin)[2];
			kinematics[3] = (*kin)[3];
			kinematics[4] = (*kin)[4];
		}

		// Copy over the charge
		TVectorT<double> * charge = (TVectorT<double>*)thisFile->Get("totalCharge");
		totalCharge[0] += (*charge)[0];
		// Set up the tree and branches
		TTree * thisTree = (TTree*) thisFile->Get("sk");


		thisTree -> SetBranchAddress("L_prl1"		,&sk_e_prl1		);
		thisTree -> SetBranchAddress("L_prl2"		,&sk_e_prl2		);

		thisTree -> SetBranchAddress("L_cer"		,&sk_e_cer		);

		thisTree -> SetBranchAddress("L_ntrk"		,&sk_e_ntrk		);
		thisTree -> SetBranchAddress("R_ntrk"		,&sk_p_ntrk		);

		thisTree -> SetBranchAddress("L_ytar"		,&sk_e_ytar		);
		thisTree -> SetBranchAddress("R_ytar"		,&sk_p_ytar		);

		thisTree -> SetBranchAddress("L_yptar"		,&sk_e_yptar		);
		thisTree -> SetBranchAddress("R_yptar"		,&sk_p_yptar		);
		thisTree -> SetBranchAddress("L_gold_yptar"	,&sk_gold_e_yptar	);
		thisTree -> SetBranchAddress("R_gold_yptar"	,&sk_gold_p_yptar	);

		thisTree -> SetBranchAddress("L_xptar"		,&sk_e_xptar		);
		thisTree -> SetBranchAddress("R_xptar"		,&sk_p_xptar		);
		thisTree -> SetBranchAddress("L_gold_xptar"	,&sk_gold_e_xptar	);
		thisTree -> SetBranchAddress("R_gold_xptar"	,&sk_gold_p_xptar	);
		
		thisTree -> SetBranchAddress("L_dp"		,&sk_e_delta		);
		thisTree -> SetBranchAddress("R_dp"		,&sk_p_delta		);
		thisTree -> SetBranchAddress("L_gold_dp"	,&sk_gold_e_delta	);
		thisTree -> SetBranchAddress("R_gold_dp"	,&sk_gold_p_delta	);
		
		thisTree -> SetBranchAddress("L_mom"		,&sk_e_mom		);
		thisTree -> SetBranchAddress("R_mom"		,&sk_p_mom		);
		thisTree -> SetBranchAddress("L_gold_mom"	,&sk_gold_e_mom		);
		thisTree -> SetBranchAddress("R_gold_mom"	,&sk_gold_p_mom		);

		thisTree -> SetBranchAddress("L_beta"		,&sk_e_beta		);
		thisTree -> SetBranchAddress("R_beta"		,&sk_p_beta		);
		thisTree -> SetBranchAddress("L_gold_beta"	,&sk_gold_e_beta	);
		thisTree -> SetBranchAddress("R_gold_beta"	,&sk_gold_p_beta	);

		thisTree -> SetBranchAddress("L_vx"		,&sk_e_x		);
		thisTree -> SetBranchAddress("L_vy"		,&sk_e_y		);
		thisTree -> SetBranchAddress("L_vz"		,&sk_e_z		);
		thisTree -> SetBranchAddress("R_vx"		,&sk_p_x		);
		thisTree -> SetBranchAddress("R_vy"		,&sk_p_y		);
		thisTree -> SetBranchAddress("R_vz"		,&sk_p_z		);
		
		thisTree -> SetBranchAddress("L_px"		,&sk_e_px		);
		thisTree -> SetBranchAddress("L_py"		,&sk_e_py		);
		thisTree -> SetBranchAddress("L_pz"		,&sk_e_pz		);
		thisTree -> SetBranchAddress("R_px"		,&sk_p_px		);
		thisTree -> SetBranchAddress("R_py"		,&sk_p_py		);
		thisTree -> SetBranchAddress("R_pz"		,&sk_p_pz		);
		
		thisTree -> SetBranchAddress("L_gold_px"	,&sk_gold_e_px		);
		thisTree -> SetBranchAddress("L_gold_py"	,&sk_gold_e_py		);
		thisTree -> SetBranchAddress("L_gold_pz"	,&sk_gold_e_pz		);
		thisTree -> SetBranchAddress("R_gold_px"	,&sk_gold_p_px		);
		thisTree -> SetBranchAddress("R_gold_py"	,&sk_gold_p_py		);
		thisTree -> SetBranchAddress("R_gold_pz"	,&sk_gold_p_pz		);

		thisTree -> SetBranchAddress("t1"		,&sk_t1			);
		thisTree -> SetBranchAddress("t4"		,&sk_t4			);
		thisTree -> SetBranchAddress("tcoinc"		,&sk_tcoinc		);
		
		thisTree -> SetBranchAddress("Pmiss"		,&sk_Pmiss		);
		thisTree -> SetBranchAddress("Emiss"		,&sk_Emiss		);
		thisTree -> SetBranchAddress("Pmiss_x"		,&sk_Pmiss_x		);
		thisTree -> SetBranchAddress("Pmiss_y"		,&sk_Pmiss_y		);
		thisTree -> SetBranchAddress("Pmiss_z"		,&sk_Pmiss_z		);
		thisTree -> SetBranchAddress("R_thetaWe"	,&sk_p_thWe		);
		thisTree -> SetBranchAddress("ph_rq"		,&sk_ph_bq		);
		thisTree -> SetBranchAddress("th_rq"		,&sk_th_bq		);
		thisTree -> SetBranchAddress("ph_xq"		,&sk_ph_xq		);
		thisTree -> SetBranchAddress("th_xq"		,&sk_th_xq		);

		thisTree -> SetBranchAddress("rast_Pmiss"	,&sk_rast_Pmiss		);
		thisTree -> SetBranchAddress("rast_Emiss"	,&sk_rast_Emiss		);
		thisTree -> SetBranchAddress("rast_Pmiss_x"	,&sk_rast_Pmiss_x	);
		thisTree -> SetBranchAddress("rast_Pmiss_y"	,&sk_rast_Pmiss_y	);
		thisTree -> SetBranchAddress("rast_Pmiss_z"	,&sk_rast_Pmiss_z	);
		thisTree -> SetBranchAddress("rast_R_thetaWe"	,&sk_p_rast_thWe	);
		thisTree -> SetBranchAddress("rast_ph_rq"	,&sk_rast_ph_bq		);
		thisTree -> SetBranchAddress("rast_th_rq"	,&sk_rast_th_bq		);
		thisTree -> SetBranchAddress("rast_ph_xq"	,&sk_rast_ph_xq		);
		thisTree -> SetBranchAddress("rast_th_xq"	,&sk_rast_th_xq		);

		thisTree -> SetBranchAddress("exRa_Pmiss"	,&sk_exRa_Pmiss		);
		thisTree -> SetBranchAddress("exRa_Emiss"	,&sk_exRa_Emiss		);
		thisTree -> SetBranchAddress("exRa_Pmiss_x"	,&sk_exRa_Pmiss_x	);
		thisTree -> SetBranchAddress("exRa_Pmiss_y"	,&sk_exRa_Pmiss_y	);
		thisTree -> SetBranchAddress("exRa_Pmiss_z"	,&sk_exRa_Pmiss_z	);
		thisTree -> SetBranchAddress("exRa_R_thetaWe"	,&sk_p_exRa_thWe	);
		thisTree -> SetBranchAddress("exRa_ph_rq"	,&sk_exRa_ph_bq		);
		thisTree -> SetBranchAddress("exRa_th_rq"	,&sk_exRa_th_bq		);
		thisTree -> SetBranchAddress("exRa_ph_xq"	,&sk_exRa_ph_xq		);
		thisTree -> SetBranchAddress("exRa_th_xq"	,&sk_exRa_th_xq		);
		
		thisTree -> SetBranchAddress("Q2"		,&sk_Q2			);
		thisTree -> SetBranchAddress("W2"		,&sk_W2			);
		thisTree -> SetBranchAddress("Nu"		,&sk_Nu			);
		thisTree -> SetBranchAddress("ph_q"		,&sk_ph_q		);
		thisTree -> SetBranchAddress("th_q"		,&sk_th_q		);
		thisTree -> SetBranchAddress("xB"		,&sk_xB			);
		thisTree -> SetBranchAddress("q3m"		,&sk_q3m		);
		thisTree -> SetBranchAddress("q_x"		,&sk_q_x		);
		thisTree -> SetBranchAddress("q_y"		,&sk_q_y		);
		thisTree -> SetBranchAddress("q_z"		,&sk_q_z		);
		thisTree -> SetBranchAddress("L_theta"		,&sk_e_th		);

		thisTree -> SetBranchAddress("rast_Q2"		,&sk_rast_Q2		);
		thisTree -> SetBranchAddress("rast_W2"		,&sk_rast_W2		);
		thisTree -> SetBranchAddress("rast_Nu"		,&sk_rast_Nu		);
		thisTree -> SetBranchAddress("rast_ph_q"	,&sk_rast_ph_q		);
		thisTree -> SetBranchAddress("rast_th_q"	,&sk_rast_th_q		);
		thisTree -> SetBranchAddress("rast_xB"		,&sk_rast_xB		);
		thisTree -> SetBranchAddress("rast_q3m"		,&sk_rast_q3m		);
		thisTree -> SetBranchAddress("rast_q_x"		,&sk_rast_q_x		);
		thisTree -> SetBranchAddress("rast_q_y"		,&sk_rast_q_y		);
		thisTree -> SetBranchAddress("rast_q_z"		,&sk_rast_q_z		);
		thisTree -> SetBranchAddress("rast_L_theta"	,&sk_e_rast_th		);

		thisTree -> SetBranchAddress("exRa_Q2"		,&sk_exRa_Q2		);
		thisTree -> SetBranchAddress("exRa_W2"		,&sk_exRa_W2		);
		thisTree -> SetBranchAddress("exRa_Nu"		,&sk_exRa_Nu		);
		thisTree -> SetBranchAddress("exRa_ph_q"	,&sk_exRa_ph_q		);
		thisTree -> SetBranchAddress("exRa_th_q"	,&sk_exRa_th_q		);
		thisTree -> SetBranchAddress("exRa_xB"		,&sk_exRa_xB		);
		thisTree -> SetBranchAddress("exRa_q3m"		,&sk_exRa_q3m		);
		thisTree -> SetBranchAddress("exRa_q_x"		,&sk_exRa_q_x		);
		thisTree -> SetBranchAddress("exRa_q_y"		,&sk_exRa_q_y		);
		thisTree -> SetBranchAddress("exRa_q_z"		,&sk_exRa_q_z		);
		thisTree -> SetBranchAddress("exRa_L_theta"	,&sk_e_exRa_th		);

		thisTree -> SetBranchAddress("L_ext_delta_dp"	,&sk_e_ext_deltaDp	);
		thisTree -> SetBranchAddress("L_ext_delta_p"	,&sk_e_ext_deltaP	);
		thisTree -> SetBranchAddress("L_ext_delta_yptar",&sk_e_ext_deltaTh	);
		thisTree -> SetBranchAddress("L_ext_dp"		,&sk_e_ext_delta	);
		thisTree -> SetBranchAddress("L_ext_mom"	,&sk_e_ext_mom		);
		thisTree -> SetBranchAddress("L_ext_yptar"	,&sk_e_ext_yptar	);
		thisTree -> SetBranchAddress("L_ext_xptar"	,&sk_e_ext_xptar	);
		thisTree -> SetBranchAddress("L_ext_px"		,&sk_e_ext_px		);
		thisTree -> SetBranchAddress("L_ext_py"		,&sk_e_ext_py		);
		thisTree -> SetBranchAddress("L_ext_pz"		,&sk_e_ext_pz		);	

		thisTree -> SetBranchAddress("R_ext_delta_dp"	,&sk_p_ext_deltaDp	);
		thisTree -> SetBranchAddress("R_ext_delta_p"	,&sk_p_ext_deltaP	);
		thisTree -> SetBranchAddress("R_ext_delta_yptar",&sk_p_ext_deltaTh	);
		thisTree -> SetBranchAddress("R_ext_dp"		,&sk_p_ext_delta	);
		thisTree -> SetBranchAddress("R_ext_mom"	,&sk_p_ext_mom		);
		thisTree -> SetBranchAddress("R_ext_yptar"	,&sk_p_ext_yptar	);
		thisTree -> SetBranchAddress("R_ext_xptar"	,&sk_p_ext_xptar	);
		thisTree -> SetBranchAddress("R_ext_px"		,&sk_p_ext_px		);
		thisTree -> SetBranchAddress("R_ext_py"		,&sk_p_ext_py		);
		thisTree -> SetBranchAddress("R_ext_pz"		,&sk_p_ext_pz		);	

		thisTree -> SetBranchAddress("Kin_L_thetaC"	,&sk_kin_eThetaC	);
		thisTree -> SetBranchAddress("Kin_R_thetaC"	,&sk_kin_pThetaC	);
		thisTree -> SetBranchAddress("Kin_L_momC"	,&sk_kin_eMomC		);
		thisTree -> SetBranchAddress("Kin_R_momC"	,&sk_kin_pMomC		);
		thisTree -> SetBranchAddress("Kin_Ebeam"	,&sk_kin_Ebeam		);
		thisTree -> SetBranchAddress("Kin_Q"		,&sk_kin_Q		);
		
		thisTree -> SetBranchAddress("BCM_curr"		,&sk_BCMcurr		);
		thisTree -> SetBranchAddress("BCM_charge"	,&sk_BCMcharge		);
		thisTree -> SetBranchAddress("BCM_isrenew"	,&sk_BCMrenew		);


		// Loop over events
		int indEv = 0; // number of events for individual file
		int indEv_0to50 = 0;
		int indEv_50to100 = 0;
		int indEv_100to150 = 0;
		int indEv_150to200 = 0;
		int indEv_200to250 = 0;
		int indEv_250to300 = 0;
		int indEv_300to350 = 0;
		int indEv_350to400 = 0;
		int indEv_400to450 = 0;
		int indEv_450to600 = 0;
		int indEv_above600 = 0;

		double nEvents = thisTree->GetEntries();
		//cout << "Beginning to loop over " << nEvents << " events.\n";
		for (int event=0 ; event < nEvents ; event++)
		{
			thisTree->GetEvent(event);
			outtree->Fill();
			
			// Now let's calculate our current status on numbers
			double xbCut, thrqCut_hi, thrqCut_lo;
			if( (abs(sk_kin_eMomC-3.54332)<0.1) && (abs(sk_kin_pMomC-1.4805)<0.1 ) && (abs(sk_kin_eThetaC*180./M_PI- 20.892)<0.1) && (abs(sk_kin_pThetaC*180./M_PI- 48.8094)<0.1) ){
				// fast cuts
				xbCut = 0.;	//1.2;
				thrqCut_hi = 40.;
				thrqCut_lo = 0.;
			}
			if( (abs(sk_kin_eMomC-3.54332)<0.1) && (abs(sk_kin_pMomC-1.246)<0.1  ) && (abs(sk_kin_eThetaC*180./M_PI- 20.892)<0.1) && (abs(sk_kin_pThetaC*180./M_PI- 58.49)<0.1 ) ){
				// slow cuts	
				xbCut = 1.2;
				thrqCut_hi = 40.;
				thrqCut_lo = 0.;
			}
			if( (abs(sk_kin_eMomC-3.54332)<0.1) && (abs(sk_kin_pMomC-1.4805)<0.1 ) && (abs(sk_kin_eThetaC*180./M_PI-17.8018)<0.1) && (abs(sk_kin_pThetaC*180./M_PI+ 48.82)<0.1)  ){
				// mid
				xbCut = 0.;
				thrqCut_hi = 360.;
				thrqCut_lo = -360.;
			}

			// Basic cuts for all kinematics
			if( pow(((sk_p_ytar+1.5*sk_p_ext_yptar)/0.08),2) + pow(((1.5*sk_p_ext_xptar)/0.08),2) <= 1 ){
			if( pow(((sk_e_ytar+1.5*sk_e_ext_yptar)/0.08),2) + pow(((1.5*sk_e_ext_xptar)/0.08),2) <= 1 ){
			if ( abs(sk_e_ext_delta)<0.045 &&
			     abs(sk_p_ext_delta)<0.045 ){
			if ( ((sk_e_prl1 + sk_e_prl2)/1000.)/(sk_e_ext_mom) > 0.8 && 
			     ((sk_e_prl1 + sk_e_prl2)/1000.)/(sk_e_ext_mom) < 1.3 ){
			if( sk_e_z > -0.1 && sk_e_z < 0.11 && 
			    sk_p_z > -0.1 && sk_p_z < 0.11 &&
			    TMath::Abs( sk_e_z - sk_p_z - 0.005) < 0.02 ){
			if( sk_tcoinc > 10) {
			if( sk_exRa_xB > xbCut ){
			if( (TMath::RadToDeg() * sk_exRa_th_bq < thrqCut_hi) && (TMath::RadToDeg() * sk_exRa_th_bq > thrqCut_lo ) ){
			
				indEv ++;
				if( (sk_exRa_Pmiss > 0.00) && (sk_exRa_Pmiss < 0.05) ) indEv_0to50++;
				if( (sk_exRa_Pmiss > 0.05) && (sk_exRa_Pmiss < 0.10) ) indEv_50to100++;
				if( (sk_exRa_Pmiss > 0.10) && (sk_exRa_Pmiss < 0.15) ) indEv_100to150++;
				if( (sk_exRa_Pmiss > 0.15) && (sk_exRa_Pmiss < 0.20) ) indEv_150to200++;
				if( (sk_exRa_Pmiss > 0.20) && (sk_exRa_Pmiss < 0.25) ) indEv_200to250++;
				if( (sk_exRa_Pmiss > 0.25) && (sk_exRa_Pmiss < 0.30) ) indEv_250to300++;
				if( (sk_exRa_Pmiss > 0.30) && (sk_exRa_Pmiss < 0.35) ) indEv_300to350++;
				if( (sk_exRa_Pmiss > 0.35) && (sk_exRa_Pmiss < 0.40) ) indEv_350to400++;
				if( (sk_exRa_Pmiss > 0.40) && (sk_exRa_Pmiss < 0.45) ) indEv_400to450++;
				if( (sk_exRa_Pmiss > 0.45) && (sk_exRa_Pmiss < 0.60) ) indEv_450to600++;
				if( (sk_exRa_Pmiss > 0.60) 			     ) indEv_above600++;
			}
			}
			}
			}
			}
			}
			}
			}
		}
		outtxt        << "\n\n\tFile: " << argv[i] << "\n\tCharge       : " << (*charge)[0] 
					      << "\n\tIndEv        : " << indEv  
					      << "\n\tIndEv0to50   : " << indEv_0to50
					      << "\n\tIndEv50to100 : " << indEv_50to100
					      << "\n\tIndEv100to150: " << indEv_100to150
					      << "\n\tIndEv150to200: " << indEv_150to200
					      << "\n\tIndEv200to250: " << indEv_200to250
					      << "\n\tIndEv250to300: " << indEv_250to300
					      << "\n\tIndEv300to350: " << indEv_300to350
					      << "\n\tIndEv350to400: " << indEv_350to400
					      << "\n\tIndEv400to450: " << indEv_400to450
					      << "\n\tIndEv450to600: " << indEv_450to600
					      << "\n\tIndEvabove600: " << indEv_above600;
		totEv 	       += indEv;                                  
		totEv_0to50  += indEv_0to50;
		totEv_50to100  += indEv_50to100;
		totEv_100to150 += indEv_100to150;
		totEv_150to200 += indEv_150to200;
		totEv_200to250 += indEv_200to250;
		totEv_250to300 += indEv_250to300;
		totEv_300to350 += indEv_300to350;
		totEv_350to400 += indEv_350to400;
		totEv_400to450 += indEv_400to450;
		totEv_450to600 += indEv_450to600;
		totEv_above600 += indEv_above600;
		thisFile->Close();
	}

			//if( TMath::RadToDeg() * 		// this is how I was calculating theta_rq before but 
			//    TMath::ACos(			// somehow it's almost correct but not exactly...?
			//	( sk_exRa_q3m * sk_exRa_q3m - 
			//		( sk_p_ext_px*sk_exRa_q_x + sk_p_ext_py*sk_exRa_q_y + sk_p_ext_pz*sk_exRa_q_z )) 
			//	/ ( sk_exRa_q3m * sk_exRa_Pmiss  ) 
			//    ) 
			//    < 40.){

	cout << "\n=================================================== \n";
	cout 			      << "\n\tTotal Charge    : " << totalCharge[0] 
        			      << "\n\tTotal Ev        : " << totEv  
        			      << "\n\tTotal Ev0to50   : " << totEv_0to50
        			      << "\n\tTotal Ev50to100 : " << totEv_50to100
        			      << "\n\tTotal Ev100to150: " << totEv_100to150
        			      << "\n\tTotal Ev150to200: " << totEv_150to200
        			      << "\n\tTotal Ev200to250: " << totEv_200to250
        			      << "\n\tTotal Ev250to300: " << totEv_250to300
        			      << "\n\tTotal Ev300to350: " << totEv_300to350
        			      << "\n\tTotal Ev350to400: " << totEv_350to400
        			      << "\n\tTotal Ev400to450: " << totEv_400to450
        			      << "\n\tTotal Ev450to600: " << totEv_450to600
        			      << "\n\tTotal Evabove600: " << totEv_above600;
	cout << "\n=================================================== \n";


	outfile->cd   ();
	outtree->Write();
	totalCharge.Write("totalCharge"); 
	kinematics.Write("kinematics");
	outfile->Close();
	outtxt.close();
	return 0;
}
