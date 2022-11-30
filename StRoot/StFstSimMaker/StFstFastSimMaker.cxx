#include "StFstFastSimMaker.h"

#include "St_base/StMessMgr.h"

#include "StEvent/StEvent.h"
#include "StEvent/StRnDHit.h"
#include "StEvent/StRnDHitCollection.h"

#include "tables/St_g2t_fts_hit_Table.h"
#include "tables/St_g2t_track_Table.h"

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarGenerator/UTIL/StarRandom.h"

#include "TCanvas.h"
#include "TCernLib.h"
#include "TH2F.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TString.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorT.h"

#include <array>
#define _USE_MATH_DEFINES
#include <cmath>
#include <map>
#include <algorithm>

#include "St_db_Maker/St_db_Maker.h"
#include "tables/St_Survey_Table.h"

// lets not polute the global scope
namespace FstGlobal{
	// converts the position error to a cov mat
    StMatrixF Hack1to6(const StHit *stHit);
	// used to convert between uniform and 1 sigma error
    constexpr float SQRT12 = sqrt(12.0);
    
    double PI = TMath::Pi();

    double PHI[] = {5.*PI/12.,  3.*PI/12.,  1.*PI/12.,  -1.*PI/12.,
                    -3.*PI/12., -5.*PI/12., -7.*PI/12., -9.*PI/12.,
                    -11.*PI/12.,-13.*PI/12.,-15.*PI/12.,-17.*PI/12.};

    // Disk segmentation
    //
    float RMIN[] = {0.95 * 4.3, 0.95 * 4.3, 0.95 * 4.3, 0.95 * 5.0, 0.95 * 5.0, 0.95 * 5.0};
    float RMAX[] = {1.05 * 15.0, 1.05 * 25.0, 1.05 * 25.0, 1.05 * 28.0, 1.05 * 28.0, 1.05 * 28.0};

    //NEXT IS only for disk ARRAY 456 with the radius from 5 to 28.
    float RSegment[] = {5., 7.875, 10.75, 13.625, 16.5, 19.375, 22.25, 25.125, 28.};

    // controls some extra output
    const bool verbose = false;

	// key type for lookup tables on disk r-phi strip
	typedef std::tuple<int, int, int> FstKeyTriple;


}

StFstFastSimMaker::StFstFastSimMaker(const Char_t *name)
	: StMaker{name},
    mNumR{8},
    mNumPHI{128},
    mNumSEC{12},
    mRaster{0},
    mInEff{0},
    mHist{false},
    mQAFileName(0),
    hTrutHitYXDisk(0),
    hTrutHitRDisk(0),
    hTrutHitRShower{0},
    hTrutHitPhiDisk(0),
    hTrutHitPhiZ(0),
    hRecoHitYXDisk(0),
    hRecoHitRDisk(0),
    hRecoHitPhiDisk(0),
    hRecoHitPhiZ(0),
    hGlobalDRDisk(0),
    hGlobalZ(0),
    h2GlobalXY(0),
    h2GlobalSmearedXY(0),
    h2GlobalDeltaXY(0),
    h3GlobalDeltaXYDisk(0),
    h3GlobalDeltaXYR(0) { }

int StFstFastSimMaker::Init() {

        //-----  LOAD ALIGNMENT MATRICES  -----//
        St_db_Maker *dbMk=new St_db_Maker("db", "MySQL:StarDb", "$STAR/StarDb");
        dbMk->SetDebug();
        dbMk->SetDateTime(20211026,7); // event or run start time, set to your liking
        dbMk->SetFlavor("ofl");

        dbMk->Init();
        dbMk->Make();
       
        Survey_st *fstOnTpc;
        Survey_st *hssOnFst;
        Survey_st *fstWedgeOnHss;
        Survey_st *fstSensorOnWedge;

        TMatrixD MfstOnTpc(4,4);
        TMatrixD MhssOnFst(4,4);
        TMatrixD MfstWedgeOnHss(4,4);
        TMatrixD MfstSensorOnWedge(4,4);

        MfstOnTpc.Zero();
        MhssOnFst.Zero();
        MfstWedgeOnHss.Zero();
        MfstSensorOnWedge.Zero();

        TDataSet *DB = 0;
        DB = dbMk->GetDataBase("Geometry/fst/fstOnTpc");
        if (!DB) {
          std::cout << "ERROR: no fstOnTpc table found in db, or malformed local db config" << std::endl;
        }

        St_Survey *dataset = 0;
        dataset = (St_Survey*) DB->Find("fstOnTpc");
        
        if (dataset) {

            Int_t rows = dataset->GetNRows();
            if (rows > 1) {
              std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
            }

            TDatime val[2];
            dbMk->GetValidity((TTable*)dataset,val);
            std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - " 
              << val[1].GetDate() << "." << val[1].GetTime() << " ] "
              << std::endl;

            fstOnTpc = dataset->GetTable();

            MfstOnTpc(0,0) = fstOnTpc[0].r00;
            MfstOnTpc(0,1) = fstOnTpc[0].r01;
            MfstOnTpc(0,2) = fstOnTpc[0].r02;
            MfstOnTpc(1,0) = fstOnTpc[0].r10;
            MfstOnTpc(1,1) = fstOnTpc[0].r11;
            MfstOnTpc(1,2) = fstOnTpc[0].r12;
            MfstOnTpc(2,0) = fstOnTpc[0].r20;
            MfstOnTpc(2,1) = fstOnTpc[0].r21;
            MfstOnTpc(2,2) = fstOnTpc[0].r22;
            MfstOnTpc(0,3) = fstOnTpc[0].t0 ;
            MfstOnTpc(1,3) = fstOnTpc[0].t1 ;
            MfstOnTpc(2,3) = fstOnTpc[0].t2 ;
            MfstOnTpc(3,3) = 1.0            ;

        } else {
            std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        }

        DB = dbMk->GetDataBase("Geometry/fst/hssOnFst");
        if (!DB) {
          std::cout << "ERROR: no hssOnFst table found in db, or malformed local db config" << std::endl;
        }

        dataset = (St_Survey*) DB->Find("hssOnFst");
        
        if (dataset) {

            Int_t rows = dataset->GetNRows();
            if (rows > 1) {
              std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
            }

            TDatime val[2];
            dbMk->GetValidity((TTable*)dataset,val);
            std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - " 
              << val[1].GetDate() << "." << val[1].GetTime() << " ] "
              << std::endl;

            hssOnFst = dataset->GetTable();
           
        } else {
            std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        }

        DB = dbMk->GetDataBase("Geometry/fst/fstWedgeOnHss");
        if (!DB) {
          std::cout << "ERROR: no fstWedgeOnHss table found in db, or malformed local db config" << std::endl;
        }

        dataset = (St_Survey*) DB->Find("fstWedgeOnHss");

        if (dataset) {

            Int_t rows = dataset->GetNRows();
            if (rows > 1) {
              std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
            }

            TDatime val[2];
            dbMk->GetValidity((TTable*)dataset,val);
            std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - " 
              << val[1].GetDate() << "." << val[1].GetTime() << " ] "
              << std::endl;

            fstWedgeOnHss = dataset->GetTable();

        } else {
            std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        }


        DB = dbMk->GetDataBase("Geometry/fst/fstSensorOnWedge");
        if (!DB) {
          std::cout << "ERROR: no fstSensorOnWedge table found in db, or malformed local db config" << std::endl;
        }

        dataset = (St_Survey*) DB->Find("fstSensorOnWedge");
        
        if (dataset) {

            Int_t rows = dataset->GetNRows();
            if (rows > 1) {
              std::cout << "INFO: found INDEXED table with " << rows << " rows" << std::endl;
            }

            TDatime val[2];
            dbMk->GetValidity((TTable*)dataset,val);
            std::cout << "Dataset validity range: [ " << val[0].GetDate() << "." << val[0].GetTime() << " - " 
              << val[1].GetDate() << "." << val[1].GetTime() << " ] "
              << std::endl;

            fstSensorOnWedge = dataset->GetTable();
        
        } else {
            std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        }


        // default FST sensor z locations, found externally from ideal geometry
        for (int is = 0; is < 108; is++) {
            
            // FST half
            int h = (is / 18) % 2; // 0 (left +x half), 1 (right -x half) 
        
            // FST wedge
            int w = is / 3; // 0-35

            // create matrices from alignment tables
            MhssOnFst(0,0) = hssOnFst[h].r00; 
            MhssOnFst(0,1) = hssOnFst[h].r01; 
            MhssOnFst(0,2) = hssOnFst[h].r02;
            MhssOnFst(1,0) = hssOnFst[h].r10; 
            MhssOnFst(1,1) = hssOnFst[h].r11; 
            MhssOnFst(1,2) = hssOnFst[h].r12;
            MhssOnFst(2,0) = hssOnFst[h].r20; 
            MhssOnFst(2,1) = hssOnFst[h].r21; 
            MhssOnFst(2,2) = hssOnFst[h].r22;
            MhssOnFst(0,3) = hssOnFst[h].t0 ;
            MhssOnFst(1,3) = hssOnFst[h].t1 ;
            MhssOnFst(2,3) = hssOnFst[h].t2 ;
            MhssOnFst(3,3) = 1.0            ;

            MfstWedgeOnHss(0,0) = fstWedgeOnHss[w].r00; 
            MfstWedgeOnHss(0,1) = fstWedgeOnHss[w].r01; 
            MfstWedgeOnHss(0,2) = fstWedgeOnHss[w].r02;
            MfstWedgeOnHss(1,0) = fstWedgeOnHss[w].r10; 
            MfstWedgeOnHss(1,1) = fstWedgeOnHss[w].r11; 
            MfstWedgeOnHss(1,2) = fstWedgeOnHss[w].r12;
            MfstWedgeOnHss(2,0) = fstWedgeOnHss[w].r20; 
            MfstWedgeOnHss(2,1) = fstWedgeOnHss[w].r21; 
            MfstWedgeOnHss(2,2) = fstWedgeOnHss[w].r22;
            MfstWedgeOnHss(0,3) = fstWedgeOnHss[w].t0 ;
            MfstWedgeOnHss(1,3) = fstWedgeOnHss[w].t1 ;
            MfstWedgeOnHss(2,3) = fstWedgeOnHss[w].t2 ;
            MfstWedgeOnHss(3,3) = 1.0                 ;

            if(is != 36)
            {
              MfstSensorOnWedge(0,0) = fstSensorOnWedge[is].r00; 
              MfstSensorOnWedge(0,1) = fstSensorOnWedge[is].r01; 
              MfstSensorOnWedge(0,2) = fstSensorOnWedge[is].r02;
              MfstSensorOnWedge(1,0) = fstSensorOnWedge[is].r10; 
              MfstSensorOnWedge(1,1) = fstSensorOnWedge[is].r11; 
              MfstSensorOnWedge(1,2) = fstSensorOnWedge[is].r12;
              MfstSensorOnWedge(2,0) = fstSensorOnWedge[is].r20; 
              MfstSensorOnWedge(2,1) = fstSensorOnWedge[is].r21; 
              MfstSensorOnWedge(2,2) = fstSensorOnWedge[is].r22;
              MfstSensorOnWedge(0,3) = fstSensorOnWedge[is].t0 ;
              MfstSensorOnWedge(1,3) = fstSensorOnWedge[is].t1 ;
              MfstSensorOnWedge(2,3) = fstSensorOnWedge[is].t2 ;
              MfstSensorOnWedge(3,3) = 1.0                     ;
            }
            else
            {       
              double cost = TMath::Cos(0.002);    
              double sint = TMath::Sin(0.002);    
              double cos75 = TMath::Cos(1.309);    
              double sin75 = TMath::Sin(1.309);  
              double dxu =  cos75*0.005; // dv = 100 um   
              double dyu =  sin75*0.005; // dv = 100 um   
              double dxv = -sin75*0.005; // dv = 100 um   
              double dyv =  cos75*0.005; // dv = 100 um   
              double dx = dxu + dxv; 
              double dy = dyu + dyv; 
              MfstSensorOnWedge(0,0) = fstSensorOnWedge[is].r00; 
              MfstSensorOnWedge(0,1) = fstSensorOnWedge[is].r01; 
              MfstSensorOnWedge(0,2) = fstSensorOnWedge[is].r02;
              MfstSensorOnWedge(1,0) = fstSensorOnWedge[is].r10; 
              MfstSensorOnWedge(1,1) = fstSensorOnWedge[is].r11; 
              MfstSensorOnWedge(1,2) = fstSensorOnWedge[is].r12;
              MfstSensorOnWedge(2,0) = fstSensorOnWedge[is].r20; 
              MfstSensorOnWedge(2,1) = fstSensorOnWedge[is].r21; 
              MfstSensorOnWedge(2,2) = fstSensorOnWedge[is].r22;
              MfstSensorOnWedge(0,3) = fstSensorOnWedge[is].t0 ;
              MfstSensorOnWedge(1,3) = fstSensorOnWedge[is].t1 ;
              MfstSensorOnWedge(2,3) = fstSensorOnWedge[is].t2 ;
              MfstSensorOnWedge(3,3) = 1.0                     ;
            }

            // Rotate and Translate plane normal vectors and origin
            TMatrixD M = MfstOnTpc * MhssOnFst * MfstWedgeOnHss * MfstSensorOnWedge;
            //M.Print();

            // save this inverse matrix for use with misaligned simulated data
            mInverseM[is].ResizeTo(4,4);
            mInverseM[is] = M.Invert();
            //mInverseM[is].Print();

            MhssOnFst.Zero();
            MfstWedgeOnHss.Zero();
            MfstSensorOnWedge.Zero();
        }


	if(mHist){
		fOut = new TFile(mQAFileName.Data(), "RECREATE");
		AddHist(hTrutHitYXDisk = new TH3F("hTrutHitYXDisk", "Global hits before segmentation", 151, -75.5, 75.5, 151, -75.5, 75.5, 10, 0, 10));
		AddHist(hTrutHitRDisk = new TH2F("hTrutHitRDisk", "Global hits before segmentation", 400, 0, 40, 10, 0, 10));
		AddHist(hTrutHitRShower[0] = new TH2F("hTrutHitRShower_4", "Global hits before segmentation", 400, 0, 40, 20, -10, 10));
		AddHist(hTrutHitRShower[1] = new TH2F("hTrutHitRShower_5", "Global hits before segmentation", 400, 0, 40, 20, -10, 10));	
           	AddHist(hTrutHitRShower[2] = new TH2F("hTrutHitRShower_6", "Global hits before segmentation", 400, 0, 40, 20, -10, 10));
                AddHist(hMCHit[0] = new TH2F("hMCHitDisk_4", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hMCHit[1] = new TH2F("hMCHitDisk_5", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hMCHit[2] = new TH2F("hMCHitDisk_6", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hMCPhiZOut[0] = new TH2F("hMCPhiZOutDisk_4", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hMCPhiZOut[1] = new TH2F("hMCPhiZOutDisk_5", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hMCPhiZOut[2] = new TH2F("hMCPhiZOutDisk_6", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hMCPhiZIn[0] = new TH2F("hMCPhiZInDisk_4", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hMCPhiZIn[1] = new TH2F("hMCPhiZInDisk_5", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hMCPhiZIn[2] = new TH2F("hMCPhiZInDisk_6", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCHit[0] = new TH2F("hRCHitDisk_4", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hRCHit[1] = new TH2F("hRCHitDisk_5", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hRCHit[2] = new TH2F("hRCHitDisk_6", "Global hits before segmentation", 150, -30.0, 30.0, 150, -30.0, 30.0));
                AddHist(hRCPhiZOut[0] = new TH2F("hRCPhiZOutDisk_4", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCPhiZOut[1] = new TH2F("hRCPhiZOutDisk_5", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCPhiZOut[2] = new TH2F("hRCPhiZOutDisk_6", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCPhiZIn[0] = new TH2F("hRCPhiZInDisk_4", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCPhiZIn[1] = new TH2F("hRCPhiZInDisk_5", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
                AddHist(hRCPhiZIn[2] = new TH2F("hRCPhiZInDisk_6", "Global hits before segmentation", 360, 0, 360, 2000, 125, 200));
		AddHist(hTrutHitPhiDisk = new TH2F("hTrutHitPhiDisk", "Global hits before segmentation", 360, 0, 360, 10, 0, 10));
		AddHist(hTrutHitPhiZ = new TH2F("hTrutHitPhiZ", "Global hits before segmentation", 360, 0, 360, 6000, 0, 600));
		AddHist(hRecoHitYXDisk = new TH3F("hRecoHitYXDisk", "Global hits after segmentation", 151, -75.5, 75.5, 151, -75.5, 75.5, 10, 0, 10));
		AddHist(hRecoHitRDisk = new TH2F("hRecoHitRDisk", "Global hits after segmentation", 400, 0, 40, 10, 0, 10));
		AddHist(hRecoHitPhiDisk = new TH2F("hRecoHitPhiDisk", "Global hits after segmentation", 360, 0, 360, 10, 0, 10));
		AddHist(hRecoHitPhiZ = new TH2F("hRecoHitPhiZ", "Global hits after segmentation", 360, 0, 360, 6000, 0, 600));
		AddHist(hGlobalDRDisk = new TH2F("hGlobalDRDisk", "; Reco. r - MC r [cm]; Events;", 1000, -50, 50, 10, 0, 10));
		AddHist(hGlobalZ = new TH1F("hGlobalZ", "; Z [cm]; Events;", 6000, 0, 600));
		AddHist(h3GlobalDeltaXYR = new TH3F("h3GlobalDeltaXYR", ";globalDeltaX; globalDeltaY; R", 300, -0.3, 0.3, 300, -3, 3, 100, 0, 30));
		AddHist(h2GlobalXY = new TH2F("h2GlobalXY", ";globalX; globalY", 1510, -75.5, 75.5, 1510, -75.5, 75.5));
		AddHist(h2GlobalSmearedXY = new TH2F("h2GlobalSmearedXY", ";globalSmearedX; globalSmearedY", 1510, -75.5, 75.5, 1510, -75.5, 75.5));
		AddHist(h2GlobalDeltaXY = new TH2F("h2GlobalDeltaXY", ";globalDeltaX; globalDeltaY", 151, -75.5, 75.5, 151, -75.5, 75.5));
		AddHist(h3GlobalDeltaXYDisk = new TH3F("h3GlobalDeltaXYDisk", ";globalDeltaX; globalDeltaY; Disk", 151, -75.5, 75.5, 151, -75.5, 75.5, 10, 0, 10));
	}
	return StMaker::Init();
}

void StFstFastSimMaker::SetDisk(const int i, const float rmn, const float rmx) {
    FstGlobal::RMIN[i] = rmn;
    FstGlobal::RMAX[i] = rmx;
}

Int_t StFstFastSimMaker::Make() {
	LOG_DEBUG << "StFstFastSimMaker::Make" << endm;

	// Get the existing StEvent, or add one if it doesn't exist.
	StEvent *event = static_cast<StEvent *>(GetDataSet("StEvent"));
	if (!event) {
		event = new StEvent;
		AddData(event);
		LOG_DEBUG << "Creating StEvent" << endm;
	}

	if (0 == event->rndHitCollection()) {
		event->setRnDHitCollection(new StRnDHitCollection());
		LOG_DEBUG << "Creating StRnDHitCollection for FTS" << endm;
	}

	// Digitize GEANT FTS hits
	FillSilicon(event);

	return kStOk;
}

/* Fill an event with StFtsHits. */
/* This should fill StFtsStrip for realistic simulator and let clustering fill StFtsHit */
/* For now skipping StFtsStrip and clustering, and fill StFtsHits directly here*/

void StFstFastSimMaker::FillSilicon(StEvent *event) {

	StRnDHitCollection *fsicollection = event->rndHitCollection();

	const int NDISC = 6;
	const int MAXR = mNumR;
	const int MAXPHI = mNumPHI * mNumSEC;

	float X0[] = {0, 0, 0, 0, 0, 0};
	float Y0[] = {0, 0, 0, 0, 0, 0};
	
	if (mRaster > 0)
		for (int i = 0; i < 6; i++) {
			X0[i] = mRaster * TMath::Cos(i * 60 * TMath::DegToRad());
			Y0[i] = mRaster * TMath::Sin(i * 60 * TMath::DegToRad());
		}

	
	// maps for hit and energy for each disk's r-phi strip
	std::map< FstGlobal::FstKeyTriple, StRnDHit* > hitMap;
	std::map< FstGlobal::FstKeyTriple, double > energySum;
	std::map< FstGlobal::FstKeyTriple, double > energyMax;

	// Read the g2t table
	St_g2t_fts_hit *hitTable = static_cast<St_g2t_fts_hit *>(GetDataSet("g2t_fsi_hit"));
	if (!hitTable) {
		LOG_INFO << "g2t_fsi_hit table is empty" << endm;
		return; // Nothing to do
	}

	const Int_t nHits = hitTable->GetNRows();
	LOG_DEBUG << "g2t_fsi_hit table has " << nHits << " hits" << endm;
	const g2t_fts_hit_st *hit = hitTable->GetTable();
	
	StPtrVecRnDHit hits;

	// track table
	St_g2t_track *trkTable = static_cast<St_g2t_track *>(GetDataSet("g2t_track"));
	if (!trkTable) {
		LOG_INFO << "g2t_track table is empty" << endm;
		return; // Nothing to do
	}

	const Int_t nTrks = trkTable->GetNRows();
	
	LOG_DEBUG << "g2t_track table has " << nTrks << " tracks" << endm;
	
	const g2t_track_st *trk = trkTable->GetTable();

	
	int count = 0;
	for (Int_t i = 0; i < nHits; ++i) {
		
		hit = (g2t_fts_hit_st *)hitTable->At(i);
		
		// skip bad hits
		if (!hit) 
			continue;


		int volume_id = hit->volume_id;
                //in misaligned geometry, this is perfectly labeled to match all mapping throught FST framework
                //in original geometry the wedges alternate. We are going to assume misaligned geometry from here on out
		LOG_DEBUG << "volume_id = " << volume_id << endm;
		int disk = volume_id / 1000;         // disk id    range = [4,6]
		int wedge = (volume_id % 1000) / 10; // wedge id   range = [1,12]
		int sensor = volume_id % 10;         // sensor id  range = [1,3]
                int globalSensorId = (disk-4)*36 + (wedge-1)*3 + (sensor-1);
		
		// used as an index for various arrays
		size_t disk_index = disk - 1;

		LOG_DEBUG << "disk = " << disk << ", wedge = " << wedge << ", sensor = " << sensor << endm;
		//cout << "disk = " << disk << ", wedge = " << wedge << ", sensor = " << sensor << endl;

		// skip non-FST hits
		if ( disk > 6 ) continue;


		double energy = hit->de;
		int track = hit->track_p;

		trk = (g2t_track_st *)trkTable->At(track);
		int isShower = false;
		if (trk)
			isShower = trk->is_shower;

		// raster coordinate offsets
		double xc = X0[disk_index];
		double yc = Y0[disk_index];

		// hit coordinates
		double x = hit->x[0];
		double y = hit->x[1];
		double z = hit->x[2];

		if (z > 200)
			continue; // skip large disks

		// rastered
		double rastered_x = x - xc;
		double rastered_y = y - yc;

                // Apply inverse misalignment matrices to find x,y in sensor ideal coordinate frame
                // ideal local sensor coordinates rotated into ideal global position
                // used to find the overall r and phi index for disk
                double gpos[4] = {x,y,z,1.};
                double rgpos[4] = {rastered_x,rastered_y,z,1.};
                TMatrixD gPosition(4,1,gpos); 
                TMatrixD r_gPosition(4,1,rgpos);
      
                //cout << "globalSensorId = " << globalSensorId << endl;
 
                //cout << "Original global position" << endl;
                //gPosition.Print();
                //r_gPosition.Print();

                //cout << "Misalignment matrix to apply" << endl;
                //mInverseM[globalSensorId].Print();

                gPosition = mInverseM[globalSensorId] * gPosition;
                r_gPosition = mInverseM[globalSensorId] * r_gPosition;             

                //cout << "Ideal global position after inverse misalignment" << endl;
                //gPosition.Print();
                //r_gPosition.Print();


		double r = sqrt(gPosition(0,0)*gPosition(0,0) + gPosition(1,0)*gPosition(1,0));
		double p = atan2(gPosition(1,0), gPosition(0,0));

		// rastered
		double rr = sqrt(r_gPosition(0,0)*r_gPosition(0,0) + r_gPosition(1,0)*r_gPosition(1,0));
		double pp = atan2(r_gPosition(1,0), r_gPosition(0,0));
                
                //cout << "r =  " << r     <<  "  phi  = " << p     << endl;
               
		// wrap an angle between 0 and 2pi
		auto wrapAngle = [&]( double angle ) {
			angle = fmod( angle, 2.0 * M_PI );
			if ( angle < 0 )
				angle += 2.0 * M_PI;
			return angle;
		};

		p = wrapAngle( p );
		pp = wrapAngle( pp );

		LOG_DEBUG << "rr = " << rr << " pp=" << pp << endm;
		LOG_DEBUG << "RMIN = " << FstGlobal::RMIN[disk_index] << " RMAX= " << FstGlobal::RMAX[disk_index] << endm;

		// Cuts made on rastered value to require the r value is within limits
		if (rr < FstGlobal::RMIN[disk_index] || rr > FstGlobal::RMAX[disk_index])
			continue;

		LOG_DEBUG << "rr = " << rr << endm;

		// Strip numbers on rastered value
		int r_index = floor(MAXR * (rr - FstGlobal::RMIN[disk_index]) / (FstGlobal::RMAX[disk_index] - FstGlobal::RMIN[disk_index]));
		
		// this gives a different conflicting answer for r_index and does not handle r outside of range
		for (int ii = 0; ii < MAXR; ii++)
			if (rr > FstGlobal::RSegment[ii] && rr <= FstGlobal::RSegment[ii + 1])
				r_index = ii;
		
                double center_phi = 5. * M_PI / 12. -  (wedge - 1) * M_PI / 6.; // radians
                double phiGap = 1.0 * M_PI / 180.; // radians               

		// Phi number
		int phi_index;
                if(sensor == 1)
                {
                  phi_index = int(MAXPHI * pp / 2.0 / M_PI);
                }
                if(sensor == 2 || sensor == 3)
                {
                  double local_phi = pp - center_phi;
                  double gapless_phi;
                  if(local_phi > 0.0) gapless_phi = pp - phiGap / 2.;
                  if(local_phi < 0.0) gapless_phi = pp + phiGap / 2.;
                  phi_index = int(MAXPHI * gapless_phi / 2.0 / M_PI);
                }
            
		if (r_index >= 8)
			continue;

		if (MAXR)
			assert(r_index < MAXR);
		if (MAXPHI)
			assert(phi_index < MAXPHI);

		StRnDHit *fsihit = nullptr;

		// key of this disk's r & phi strip
		auto threeKey = std::tie( disk_index, r_index, phi_index );
		
		if (hitMap.count( threeKey ) == 0) { // New hit

			if (FstGlobal::verbose){
				LOG_INFO << Form("NEW d=%1d xyz=%8.4f %8.4f %8.4f r=%8.4f phi=%8.4f iR=%2d iPhi=%4d dE=%8.4f[MeV] truth=%d",
						disk, x, y, z, r, p, r_index, phi_index, energy * 1000.0, track)
					<< endm;
			}

			count++;
			fsihit = new StRnDHit();
			fsihit->setDetectorId(kFtsId);
			fsihit->setLayer(disk);
                        fsihit->setVolumeId(globalSensorId); // range = [0,127]             

			//
			// Set position and position error based on radius-constant bins
			//
			double p0  = (phi_index + 0.5) * 2.0 * M_PI / double(MAXPHI);
                        if(sensor == 2 || sensor == 3)
                        {
                          double local_phi = pp - center_phi;
                          if(local_phi > 0.0) p0 += phiGap / 2.;
                          if(local_phi < 0.0) p0 -= phiGap / 2.;
                        }
			double dp = 2.0 * M_PI / double(MAXPHI) / FstGlobal::SQRT12;
			
			// ONLY valid for the disk array 456, no difference for each disk
			double r0 = (FstGlobal::RSegment[r_index] + FstGlobal::RSegment[r_index + 1]) * 0.5;
			double dr = FstGlobal::RSegment[r_index + 1] - FstGlobal::RSegment[r_index];
			
			double xtemp = r0 * cos(p0) + xc;
			double ytemp = r0 * sin(p0) + yc;

                        //cout << "Rasterized position" << endl;
                        //cout << "x =  " << xtemp <<  "    y  = " << ytemp << endl;
                        //cout << "r =  " << r     <<  "  phi  = " << p     << endl;
                        //cout << "r0 = " << r0    <<  "  phi0 = " << p0    << endl;

                        // Apply last rotation to the hit to place it in the sensor's local coordinate frame
                        //double angle = (TMath::Pi() * 3./8.) - wedge * (TMath::Pi() * 1./4.);
                        // Inverse rotation
                
                        double phiAxisShift = 0.0;
                        if((disk-4) == 0 || (disk-4) == 2)
                        {
                          if((wedge-1)%2 == 0)
                          {
                            if((sensor-1)%3 == 1)      phiAxisShift = +8.0 * M_PI / 180.0;
                            else if((sensor-1)%3 == 2) phiAxisShift = -8.0 * M_PI / 180.0;
                          }
                          else if((wedge-1)%2 == 1)
                          {
                            if((sensor-1)%3 == 1)      phiAxisShift = -8.0 * M_PI / 180.0;
                            else if((sensor-1)%3 == 2) phiAxisShift = +8.0 * M_PI / 180.0;
                          } 
                        }
                        else if((disk-4) == 1)
                        {
                          if((wedge-1)%2 == 0)
                          {
                            if((sensor-1)%3 == 1)      phiAxisShift = -8.0 * M_PI / 180.0;
                            else if((sensor-1)%3 == 2) phiAxisShift = +8.0 * M_PI / 180.0;
                          }
                          else if((wedge-1)%2 == 1)
                          {
                            if((sensor-1)%3 == 1)      phiAxisShift = +8.0 * M_PI / 180.0;
                            else if((sensor-1)%3 == 2) phiAxisShift = -8.0 * M_PI / 180.0;
                          } 
                        }
                         
                        phiAxisShift = 0.0;

                        double x0 = xtemp *       TMath::Cos(FstGlobal::PHI[wedge-1] + phiAxisShift) + ytemp * TMath::Sin(FstGlobal::PHI[wedge-1] + phiAxisShift);
                        double y0 = xtemp * -1. * TMath::Sin(FstGlobal::PHI[wedge-1] + phiAxisShift) + ytemp * TMath::Cos(FstGlobal::PHI[wedge-1] + phiAxisShift);
            
                        //cout << "Final Position of hit after rotation" << endl;
                        //cout << "x = " << x0 << "   y = " << y0 << endl;

			assert(TMath::Abs(x0) + TMath::Abs(y0) > 0);
			double dz = 0.03 / FstGlobal::SQRT12;
			double er = dr / FstGlobal::SQRT12;
			fsihit->setPosition(StThreeVectorF(x0, y0, z));
			
			fsihit->setPositionError(StThreeVectorF(er, dp, dz));
			// set covariance matrix
			fsihit->setErrorMatrix(&FstGlobal::Hack1to6(fsihit)[0][0]);

			fsihit->setCharge(energy);
			fsihit->setIdTruth(track, 100);
			hits.push_back(fsihit);
			
			hitMap[ threeKey ] = fsihit;
			energySum[ threeKey ] = energy;
			energyMax[ threeKey ] = energy;

			if (FstGlobal::verbose){
				LOG_INFO << Form("NEW d=%1d xyz=%8.4f %8.4f %8.4f ", disk, x, y, z) << endm;
				LOG_INFO << Form("smeared xyz=%8.4f %8.4f %8.4f ", fsihit->position().x(), fsihit->position().y(), fsihit->position().z()) << endm;
			}

			if(mHist){
				TVector2 hitpos_mc(x, y);
				TVector2 hitpos_rc(fsihit->position().x(), fsihit->position().y());

				hTrutHitYXDisk->Fill(x, y, disk);
				hTrutHitRDisk->Fill(hitpos_mc.Mod(), disk);
				
				if (disk == 4)
					hTrutHitRShower[0]->Fill(hitpos_mc.Mod(), isShower);
				if (disk == 5)
					hTrutHitRShower[1]->Fill(hitpos_mc.Mod(), isShower);
				if (disk == 6)
					hTrutHitRShower[2]->Fill(hitpos_mc.Mod(), isShower);

				hTrutHitPhiDisk->Fill(hitpos_mc.Phi() * 180.0 / TMath::Pi(), disk);
				hTrutHitPhiZ->Fill(hitpos_mc.Phi() * 180.0 / TMath::Pi(), z);
				hRecoHitYXDisk->Fill(fsihit->position().x(), fsihit->position().y(), disk);
				hRecoHitRDisk->Fill(hitpos_rc.Mod(), disk);
				hRecoHitPhiDisk->Fill(hitpos_rc.Phi() * 180.0 / TMath::Pi(), disk);
				hRecoHitPhiZ->Fill(hitpos_rc.Phi() * 180.0 / TMath::Pi(), z);
				hGlobalDRDisk->Fill(hitpos_rc.Mod() - hitpos_mc.Mod(), disk);
				hGlobalZ->Fill(fsihit->position().z());

				// cout << "CHECK : " << fsihit->position().x()-x << " |||  "<<  fsihit->position().y()-y << endl;
				h2GlobalXY->Fill(x, y);
				h2GlobalSmearedXY->Fill(fsihit->position().x(), fsihit->position().y());
				h2GlobalDeltaXY->Fill(fsihit->position().x() - x, fsihit->position().y() - y);
				h3GlobalDeltaXYDisk->Fill(fsihit->position().x() - x, fsihit->position().y() - y, disk);

				h3GlobalDeltaXYR->Fill(fsihit->position().x() - x, fsihit->position().y() - y, sqrt(pow(fsihit->position().x(), 2) + pow(fsihit->position().y(), 2)));
                          
                                double rrc = hitpos_rc.Mod();
                                double prc = hitpos_rc.Phi();
                                double rmc = hitpos_mc.Mod();
                                double pmc = hitpos_mc.Phi();

                                hMCHit[disk-4]->Fill(x,y);
                                hRCHit[disk-4]->Fill(fsihit->position().x(), fsihit->position().y());

                                if (sensor == 1){
                                  hMCPhiZIn[disk-4]->Fill(pmc * 180.0 / TMath::Pi(), z);
                                  hRCPhiZIn[disk-4]->Fill(prc * 180.0 / TMath::Pi(), z); 
                                }
                                else {
                                  hMCPhiZOut[disk-4]->Fill(pmc * 180.0 / TMath::Pi(), z);
                                  hRCPhiZOut[disk-4]->Fill(prc * 180.0 / TMath::Pi(), z);
                                }
			}
		}
		else { // Hit on this strip already exists, adding energy to old hit
			// get hit from the map
			fsihit = hitMap[ threeKey ] ;
			fsihit->setCharge(fsihit->charge() + energy);

			// Add energy to running sum
			energySum[ threeKey ] = energySum[ threeKey ] + energy;

			if (energy> energyMax[ threeKey ])
				energyMax[ threeKey ] = energy;

			// keep idtruth but dilute it...
			track = fsihit->idTruth();

			fsihit->setIdTruth(track, 100 * (energyMax[ threeKey ] / energySum[ threeKey ]));
		}
	} // loop on hits
	int nfsihit = hits.size();

	StarRandom &rand = StarRandom::Instance();

	// NOW run back through the hits and add them if they pass an efficiency roll
	for (int i = 0; i < nfsihit; i++) {
		double rnd_save = rand.flat();
		if (rnd_save > mInEff){
			fsicollection->addHit(hits[i]);
		}
	}
	if (FstGlobal::verbose) {
		LOG_DEBUG << Form("Found %d/%d g2t hits in %d cells, created %d hits with ADC>0", count, nHits, nfsihit, fsicollection->numberOfHits()) << endm;
	}

}
//

int StFstFastSimMaker::Finish() {
	if(mHist){
		fOut->cd();
		hTrutHitYXDisk->Write();
		hTrutHitRDisk->Write();
		hTrutHitRShower[0]->Write();
		hTrutHitRShower[1]->Write();
		hTrutHitRShower[2]->Write();
		hTrutHitPhiDisk->Write();
		hTrutHitPhiZ->Write();
		hRecoHitYXDisk->Write();
		hRecoHitRDisk->Write();
		hRecoHitPhiDisk->Write();
		hRecoHitPhiZ->Write();
		hGlobalDRDisk->Write();
		hGlobalZ->Write();
		h3GlobalDeltaXYR->Write();
		h2GlobalXY->Write();
		h2GlobalSmearedXY->Write();
		h2GlobalDeltaXY->Write();
		h3GlobalDeltaXYDisk->Write();
                hMCHit[0]->Write();
                hMCHit[1]->Write();
                hMCHit[2]->Write();
                hMCPhiZOut[0]->Write();
                hMCPhiZOut[1]->Write();
                hMCPhiZOut[2]->Write();
                hMCPhiZIn[0]->Write();
                hMCPhiZIn[1]->Write();
                hMCPhiZIn[2]->Write();
                hRCHit[0]->Write();
                hRCHit[1]->Write();
                hRCHit[2]->Write();
                hRCPhiZOut[0]->Write();
                hRCPhiZOut[1]->Write();
                hRCPhiZOut[2]->Write();
                hRCPhiZIn[0]->Write();
                hRCPhiZIn[1]->Write();
                hRCPhiZIn[2]->Write();
		fOut->Close();
	}
	return kStOK;
}

//_____________________________________________________________________________
StMatrixF FstGlobal::Hack1to6(const StHit *stHit) {
	//   X = R*cos(Fi), Y=R*sin(Fi), Z = z
	//   dX/dR  = (    cos(Fi)  ,sin(Fi),0)
	//   dX/dFi = (-R*sin(Fi), R*cos(Fi),0)
	//   dX/dZ  = (         0,         0,1)

	auto hiPos = stHit->position();
	auto hiErr = stHit->positionError();
	double Rxy = sqrt(hiPos[0] * hiPos[0] + hiPos[1] * hiPos[1]);
	double cosFi = hiPos[0] / Rxy;
	double sinFi = hiPos[1] / Rxy;
	double T[3][3] = {{cosFi, -Rxy * sinFi, 0}, {sinFi, Rxy * cosFi, 0}, {0, 0, 1}};
	double Ginp[6] = {hiErr[0] * hiErr[0], 0, hiErr[1] * hiErr[1], 0, 0, hiErr[2] * hiErr[2]};
	double Gout[6];

	TCL::trasat(T[0], Ginp, Gout, 3, 3);
	StMatrixF mtxF(3, 3);

	for (int i = 0, li = 0; i < 3; li += ++i) {
		for (int j = 0; j <= i; j++) {
			mtxF[i][j] = Gout[li + j];
			mtxF[j][i] = mtxF[i][j];
		}
	}

	return mtxF;
}
