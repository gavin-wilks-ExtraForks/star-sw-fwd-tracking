#ifndef TRACK_FITTER_H
#define TRACK_FITTER_H

#include "GenFit/ConstField.h"
#include "GenFit/EventDisplay.h"
#include "GenFit/Exception.h"
#include "GenFit/FieldManager.h"
#include "GenFit/KalmanFitStatus.h"
#include "GenFit/GblFitter.h"
#include "GenFit/GFGbl.h"
#include "StFwdTrackMaker/StFwdGbl.h"

#include "GenFit/KalmanFitter.h"
#include "GenFit/KalmanFitterInfo.h"
#include "GenFit/KalmanFitterRefTrack.h"
#include "GenFit/MaterialEffects.h"
#include "GenFit/PlanarMeasurement.h"
#include "GenFit/RKTrackRep.h"
#include "GenFit/SpacepointMeasurement.h"
#include "GenFit/StateOnPlane.h"
#include "GenFit/TGeoMaterialInterface.h"
#include "GenFit/Track.h"
#include "GenFit/TrackPoint.h"

#include "Criteria/SimpleCircle.h"

#include "TDatabasePDG.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TProfile.h"

#include <vector>
#include <memory>

#include "StFwdTrackMaker/Common.h"

#include "StFwdTrackMaker/include/Tracker/FwdHit.h"
#include "StFwdTrackMaker/include/Tracker/TrackFitter.h"
#include "StFwdTrackMaker/include/Tracker/STARField.h"
#include "StFwdTrackMaker/include/Tracker/FwdGeomUtils.h"

#include "StarGenerator/UTIL/StarRandom.h"

#include "St_base/StMessMgr.h"
#include "St_db_Maker/St_db_Maker.h"

#include "tables/St_Survey_Table.h"
#include "StMaker.h"
//#include "TGeoHMatrix.h"

/* Cass for fitting tracks(space points) with GenFit
 *
 */
class TrackFitter {

// Accessors and options
  public:
    genfit::FitStatus getStatus() { return mFitStatus; }
    genfit::AbsTrackRep *getTrackRep() { return mTrackRep; }
    genfit::Track *getTrack() { return mFitTrack; }
    void setGenerateHistograms( bool gen) { mGenHistograms = gen;}

  public:
    // ctor 
    // provide the main configuration object
    TrackFitter(FwdTrackerConfig _mConfig) : mConfig(_mConfig) {
        mTrackRep = 0;
        mFitTrack = 0;
    }


    void setup() {

        // the geometry manager that GenFit will use
        TGeoManager * gMan = nullptr;
        mMisaligned = mConfig.get<bool>("Geometry:misaligned",false);

        // Setup the Geometry used by GENFIT
        TGeoManager::Import(mConfig.get<string>("Geometry", "fGeom.root").c_str());
        gMan = gGeoManager;
        // Set up the material interface and set material effects on/off from the config
        genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
        genfit::MaterialEffects::getInstance()->setNoEffects(mConfig.get<bool>("TrackFitter::noMaterialEffects", false)); // false means defaul ON

        // Determine which Magnetic field to use
        // Either constant field or real field from StarFieldAdaptor
        if (mConfig.get<bool>("TrackFitter:constB", false)) {
            mBField = std::unique_ptr<genfit::AbsBField>(new genfit::ConstField(0., 0., 5.)); // 0.5 T Bz
            // mBField = std::unique_ptr<genfit::AbsBField>(new genfit::ConstField(0., 0., 0.)); // ZERO FIELD
            LOG_INFO << "StFwdTrackMaker: Tracking with constant magnetic field" << endl;
        } else if (mConfig.get<bool>("TrackFitter:zeroB", false)) {
            mBField = std::unique_ptr<genfit::AbsBField>(new genfit::ConstField(0., 0., 0.)); // ZERO FIELD
            LOG_INFO << "StFwdTrackMaker: Tracking with ZERO magnetic field" << endl;
        } else {
            mBField = std::unique_ptr<genfit::AbsBField>(new StarFieldAdaptor());
            LOG_INFO << "StFwdTrackMaker: Tracking with StarFieldAdapter" << endl;
        }
        // we must have one of the two available fields at this point
        // note, the pointer is still bound to the lifetime of the TackFitter
        genfit::FieldManager::getInstance()->init(mBField.get()); 

        // initialize the main mFitter using a KalmanFitter with reference tracks
        mFitter = std::unique_ptr<genfit::AbsKalmanFitter>(new genfit::KalmanFitterRefTrack());

        // Here we load several options from the config, 
        // to customize the mFitter behavior
        mFitter->setMaxFailedHits(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MaxFailedHits", -1)); // default -1, no limit
        mFitter->setDebugLvl(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:DebugLvl", 0)); // default 0, no output
        mFitter->setMaxIterations(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MaxIterations", 4)); // default 4 iterations
        mFitter->setMinIterations(mConfig.get<int>("TrackFitter.KalmanFitterRefTrack:MinIterations", 0)); // default 0 iterations

        if(mConfig.get<bool>("TrackFitter:refitGBL",false)) {
            mRefitGBL = true;
            mGblFitter = std::unique_ptr<StFwdGbl>(new StFwdGbl());
            mGblFitter->setGBLOptions("HH",true,true);
            mGblFitter->setDebugLvl(1);
            mGblFitter->beginRun();
        }
        // FwdGeomUtils looks into the loaded geometry and gets detector z locations if present
        FwdGeomUtils fwdGeoUtils( gMan );

        // these default values are the default if the detector is 
        // a) not found in the geometry 
        // b) not provided in config

        // NOTE: these defaults are needed since the geometry file might not include FST (bug being worked on separately)
        //mFSTZLocations = fwdGeoUtils.fstZ(
        //    mConfig.getVector<double>("TrackFitter.Geometry:fst", 
        //        {140.286011, 154.286011, 168.286011 }
        //        // 144.633,158.204,171.271
        //    )
        //);
        
        mFSTZLocations.push_back(151.750);
        mFSTZLocations.push_back(165.248);
        mFSTZLocations.push_back(178.781);

        if ( fwdGeoUtils.fstZ( 0 ) < 1.0 ) { // returns 0.0 on failure
            LOG_WARN << "Using FST z-locations from config or defautl, may not match hits" << endm;
        }

        const double dzInnerFst = 1.715 + 0.04; // cm relative to "center" of disk + residual...
        const double dzOuterFst = 0.240 + 0.04; // cm relative to "center" of disk

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

            //std::cout << "Id:        " << fstOnTpc[0].Id        << std::endl;
            //std::cout << "r00:       " << fstOnTpc[0].r00       << std::endl;
            //std::cout << "r01:       " << fstOnTpc[0].r01       << std::endl;
            //std::cout << "r02:       " << fstOnTpc[0].r02       << std::endl;
            //std::cout << "r10:       " << fstOnTpc[0].r10       << std::endl;
            //std::cout << "r11:       " << fstOnTpc[0].r11       << std::endl;
            //std::cout << "r12:       " << fstOnTpc[0].r12       << std::endl;
            //std::cout << "r20:       " << fstOnTpc[0].r20       << std::endl;
            //std::cout << "r21:       " << fstOnTpc[0].r21       << std::endl;
            //std::cout << "r22:       " << fstOnTpc[0].r22       << std::endl;
            //std::cout << "t0:        " << fstOnTpc[0].t0        << std::endl;
            //std::cout << "t1:        " << fstOnTpc[0].t1        << std::endl;
            //std::cout << "t2:        " << fstOnTpc[0].t2        << std::endl;
            //std::cout << "sigmaRotX: " << fstOnTpc[0].sigmaRotX << std::endl;
            //std::cout << "sigmaRotY: " << fstOnTpc[0].sigmaRotY << std::endl;
            //std::cout << "sigmaRotZ: " << fstOnTpc[0].sigmaRotZ << std::endl;
            //std::cout << "sigmaTrX:  " << fstOnTpc[0].sigmaTrX  << std::endl;
            //std::cout << "sigmaTrY:  " << fstOnTpc[0].sigmaTrY  << std::endl;
            //std::cout << "sigmaTrZ:  " << fstOnTpc[0].sigmaTrZ  << std::endl;
            //std::cout << "comment:   " << fstOnTpc[0].comment   << std::endl;

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
           
            //for(int i = 0; i < 2; i++) {
            //    std::cout << "Id:        " << hssOnFst[i].Id        << std::endl;
            //    std::cout << "r00:       " << hssOnFst[i].r00       << std::endl;
            //    std::cout << "r01:       " << hssOnFst[i].r01       << std::endl;
            //    std::cout << "r02:       " << hssOnFst[i].r02       << std::endl;
            //    std::cout << "r10:       " << hssOnFst[i].r10       << std::endl;
            //    std::cout << "r11:       " << hssOnFst[i].r11       << std::endl;
            //    std::cout << "r12:       " << hssOnFst[i].r12       << std::endl;
            //    std::cout << "r20:       " << hssOnFst[i].r20       << std::endl;
            //    std::cout << "r21:       " << hssOnFst[i].r21       << std::endl;
            //    std::cout << "r22:       " << hssOnFst[i].r22       << std::endl;
            //    std::cout << "t0:        " << hssOnFst[i].t0        << std::endl;
            //    std::cout << "t1:        " << hssOnFst[i].t1        << std::endl;
            //    std::cout << "t2:        " << hssOnFst[i].t2        << std::endl;
            //    std::cout << "sigmaRotX: " << hssOnFst[i].sigmaRotX << std::endl;
            //    std::cout << "sigmaRotY: " << hssOnFst[i].sigmaRotY << std::endl;
            //    std::cout << "sigmaRotZ: " << hssOnFst[i].sigmaRotZ << std::endl;
            //    std::cout << "sigmaTrX:  " << hssOnFst[i].sigmaTrX  << std::endl;
            //    std::cout << "sigmaTrY:  " << hssOnFst[i].sigmaTrY  << std::endl;
            //    std::cout << "sigmaTrZ:  " << hssOnFst[i].sigmaTrZ  << std::endl;
            //    std::cout << "comment:   " << hssOnFst[i].comment   << std::endl;
            //}
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

            //for(int i = 0; i < 36; i++) {
            //    std::cout << "Id:        " << fstWedgeOnHss[i].Id        << std::endl;
            //    std::cout << "r00:       " << fstWedgeOnHss[i].r00       << std::endl;
            //    std::cout << "r01:       " << fstWedgeOnHss[i].r01       << std::endl;
            //    std::cout << "r02:       " << fstWedgeOnHss[i].r02       << std::endl;
            //    std::cout << "r10:       " << fstWedgeOnHss[i].r10       << std::endl;
            //    std::cout << "r11:       " << fstWedgeOnHss[i].r11       << std::endl;
            //    std::cout << "r12:       " << fstWedgeOnHss[i].r12       << std::endl;
            //    std::cout << "r20:       " << fstWedgeOnHss[i].r20       << std::endl;
            //    std::cout << "r21:       " << fstWedgeOnHss[i].r21       << std::endl;
            //    std::cout << "r22:       " << fstWedgeOnHss[i].r22       << std::endl;
            //    std::cout << "t0:        " << fstWedgeOnHss[i].t0        << std::endl;
            //    std::cout << "t1:        " << fstWedgeOnHss[i].t1        << std::endl;
            //    std::cout << "t2:        " << fstWedgeOnHss[i].t2        << std::endl;
            //    std::cout << "sigmaRotX: " << fstWedgeOnHss[i].sigmaRotX << std::endl;
            //    std::cout << "sigmaRotY: " << fstWedgeOnHss[i].sigmaRotY << std::endl;
            //    std::cout << "sigmaRotZ: " << fstWedgeOnHss[i].sigmaRotZ << std::endl;
            //    std::cout << "sigmaTrX:  " << fstWedgeOnHss[i].sigmaTrX  << std::endl;
            //    std::cout << "sigmaTrY:  " << fstWedgeOnHss[i].sigmaTrY  << std::endl;
            //    std::cout << "sigmaTrZ:  " << fstWedgeOnHss[i].sigmaTrZ  << std::endl;
            //    std::cout << "comment:   " << fstWedgeOnHss[i].comment   << std::endl;
            //}
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
            
            for(int i = 0; i < 108; i++) {
                std::cout << "Id:        " << fstSensorOnWedge[i].Id        << std::endl;
                std::cout << "r00:       " << fstSensorOnWedge[i].r00       << std::endl;
                std::cout << "r01:       " << fstSensorOnWedge[i].r01       << std::endl;
                std::cout << "r02:       " << fstSensorOnWedge[i].r02       << std::endl;
                std::cout << "r10:       " << fstSensorOnWedge[i].r10       << std::endl;
                std::cout << "r11:       " << fstSensorOnWedge[i].r11       << std::endl;
                std::cout << "r12:       " << fstSensorOnWedge[i].r12       << std::endl;
                std::cout << "r20:       " << fstSensorOnWedge[i].r20       << std::endl;
                std::cout << "r21:       " << fstSensorOnWedge[i].r21       << std::endl;
                std::cout << "r22:       " << fstSensorOnWedge[i].r22       << std::endl;
                std::cout << "t0:        " << fstSensorOnWedge[i].t0        << std::endl;
                std::cout << "t1:        " << fstSensorOnWedge[i].t1        << std::endl;
                std::cout << "t2:        " << fstSensorOnWedge[i].t2        << std::endl;
                std::cout << "sigmaRotX: " << fstSensorOnWedge[i].sigmaRotX << std::endl;
                std::cout << "sigmaRotY: " << fstSensorOnWedge[i].sigmaRotY << std::endl;
                std::cout << "sigmaRotZ: " << fstSensorOnWedge[i].sigmaRotZ << std::endl;
                std::cout << "sigmaTrX:  " << fstSensorOnWedge[i].sigmaTrX  << std::endl;
                std::cout << "sigmaTrY:  " << fstSensorOnWedge[i].sigmaTrY  << std::endl;
                std::cout << "sigmaTrZ:  " << fstSensorOnWedge[i].sigmaTrZ  << std::endl;
                std::cout << "comment:   " << fstSensorOnWedge[i].comment   << std::endl;
            }
        } else {
            std::cout << "ERROR: dataset does not contain requested table" << std::endl;
        }

        //delete DB;
        //delete dataset;

        // Now add the Si detector planes at the desired location
        std::stringstream sstr;
        sstr << "Adding FST Planes at: ";
        string delim = "";
        for (auto z : mFSTZLocations) {
            mFSTPlanes.push_back(
                genfit::SharedPlanePtr(
                    // these normals make the planes face along z-axis
                    new genfit::DetPlane(TVector3(0, 0, z), TVector3(1, 0, 0), TVector3(0, 1, 0) )
                )
            );

            //// Inner Module FST planes
            //mFSTPlanesInner.push_back(
            //    genfit::SharedPlanePtr(
            //        // these normals make the planes face along z-axis
            //        new genfit::DetPlane(TVector3(0, 0, z - dzInnerFst), TVector3(1, 0, 0), TVector3(0, 1, 0) )
            //    )
            //);
            //mFSTPlanesInner.push_back(
            //    genfit::SharedPlanePtr(
            //        // these normals make the planes face along z-axis
            //        new genfit::DetPlane(TVector3(0, 0, z + dzInnerFst), TVector3(1, 0, 0), TVector3(0, 1, 0) )
            //    )
            //);
            //// Outer Module FST planes
            //mFSTPlanesOuter.push_back(
            //    genfit::SharedPlanePtr(
            //        // these normals make the planes face along z-axis
            //        new genfit::DetPlane(TVector3(0, 0, z - dzOuterFst), TVector3(1, 0, 0), TVector3(0, 1, 0) )
            //    )
            //);
            //mFSTPlanesOuter.push_back(
            //    genfit::SharedPlanePtr(
            //        // these normals make the planes face along z-axis
            //        new genfit::DetPlane(TVector3(0, 0, z + dzOuterFst), TVector3(1, 0, 0), TVector3(0, 1, 0) )
            //    )
            //);

            //sstr << delim << z << " (-dzInner=" << z - dzInnerFst << ", +dzInner=" << z+dzInnerFst << ", -dzOuter=" << z - dzOuterFst << ", +dzOuter=" << z + dzOuterFst << ")";
            //delim = ", ";
        }
 
        // default FST sensor z locations, found externally from ideal geometry
        double fstDefaultZ[12] = {150.008101,151.403100,153.491899,152.096900,166.989901,165.594901,163.506101,164.901101,177.039106,178.434106,180.522905,179.127906};
        for (int is = 0; is < 108; is++) {
            
            double z;

            // FST half
            int h = (is / 18) % 2; // 0 (left +x half), 1 (right -x half) 
        
            // FST disk (integer division rounds down)
            int d = is / 36; // 0-2

            // FST wedge
            int w = is / 3; // 0-35

            // FST sensor
            int s = is % 3; // 0 (inner), 1 (outer), 2 (outer)
            int ds = (s == 0)? 0 : 1; // +0 for inner, +1 for outer
       
            int defaultZidx = d * 4 + 2 * (w % 2) + ds;
            z = fstDefaultZ[defaultZidx];


            double angle = TMath::Pi()*5./12. - double(w)*TMath::Pi()*1./6.;

            double phiAxisShift = 0.0;
            if(d == 0 || d == 2)
            {
              if(w%2 == 0)
              {
                if(s == 1)      phiAxisShift = +8.0 * TMath::Pi() / 180.0;
                else if(s == 2) phiAxisShift = -8.0 * TMath::Pi() / 180.0;
              }
              else if(w%2 == 1)
              {
                if(s == 1)      phiAxisShift = -8.0 * TMath::Pi() / 180.0;
                else if(s == 2) phiAxisShift = +8.0 * TMath::Pi() / 180.0;
              } 
            }
            else if(d == 1)
            {
              if(w%2 == 0)
              {
                if(s == 1)      phiAxisShift = -8.0 * TMath::Pi() / 180.0;
                else if(s == 2) phiAxisShift = +8.0 * TMath::Pi() / 180.0;
              }
              else if(w%2 == 1)
              {
                if(s == 1)      phiAxisShift = +8.0 * TMath::Pi() / 180.0;
                else if(s == 2) phiAxisShift = -8.0 * TMath::Pi() / 180.0;
              } 
            }
            
            //angle += phiAxisShift;

            TMatrixD RotZ(4,4);
            RotZ(0,0) =  TMath::Cos(angle);          
            RotZ(0,1) = -TMath::Sin(angle);          
            RotZ(1,0) =  TMath::Sin(angle);          
            RotZ(1,1) =  TMath::Cos(angle);          
            RotZ(2,2) =  1.0              ;          
            RotZ(3,3) =  1.0              ;          
  
            double ua[4] = {1.,0.,0.,1.};
            double va[4] = {0.,1.,0.,1.};
            double oa[4] = {0.,0.,0.,1.};
            TMatrixD u4(4,1,ua); // default u (corresponds to x on the sensor, (x,0,0) = (r,0,0) along x-axis) 
            TMatrixD v4(4,1,va); // default v (corresponds to y on the sensor)
            TMatrixD o4(4,1,oa); // default origin

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

            //if(is != 36)
            //{
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
            //}
            //else
            //{
            //  double cost = TMath::Cos(0.001);    
            //  double sint = TMath::Sin(0.001);    
            //  double cos75 = TMath::Cos(1.309);    
            //  double sin75 = TMath::Sin(1.309);  
            //  double dxu =  cos75*0.01; // dv = 100 um   
            //  double dyu =  sin75*0.01; // dv = 100 um   
            //  double dxv = -sin75*0.01; // dv = 100 um   
            //  double dyv =  cos75*0.01; // dv = 100 um   
            //  double dx = dxu + dxv; 
            //  double dy = dyu + dyv; 
            //  MfstSensorOnWedge(0,0) = cost;//fstSensorOnWedge[is].r00; 
            //  MfstSensorOnWedge(0,1) = -sint;//fstSensorOnWedge[is].r01; 
            //  MfstSensorOnWedge(0,2) = fstSensorOnWedge[is].r02;
            //  MfstSensorOnWedge(1,0) = sint;//fstSensorOnWedge[is].r10; 
            //  MfstSensorOnWedge(1,1) = cost;//fstSensorOnWedge[is].r11; 
            //  MfstSensorOnWedge(1,2) = fstSensorOnWedge[is].r12;
            //  MfstSensorOnWedge(2,0) = fstSensorOnWedge[is].r20; 
            //  MfstSensorOnWedge(2,1) = fstSensorOnWedge[is].r21; 
            //  MfstSensorOnWedge(2,2) = fstSensorOnWedge[is].r22;
            //  MfstSensorOnWedge(0,3) = dxv;//fstSensorOnWedge[is].t0 ;
            //  MfstSensorOnWedge(1,3) = dyv;//fstSensorOnWedge[is].t1 ;
            //  MfstSensorOnWedge(2,3) = fstSensorOnWedge[is].t2 ;
            //  MfstSensorOnWedge(3,3) = 1.0                     ;
            //}

            // Rotate and Translate plane normal vectors and origin
            TMatrixD M = MfstOnTpc * MhssOnFst * MfstWedgeOnHss * MfstSensorOnWedge * RotZ;
            u4 = M * u4;           
            v4 = M * v4;            
            o4 = M * o4;            

            // save this inverse matrix for use with misaligned simulated data
            mInverseM[is].ResizeTo(4,4);
            mInverseM[is] = M.Invert();

            TVector3 u(u4(0,0),u4(1,0),u4(2,0));
            TVector3 v(v4(0,0),v4(1,0),v4(2,0));
            TVector3 o(o4(0,0),o4(1,0),o4(2,0));

            mFSTSensorPlanes.push_back(
                genfit::SharedPlanePtr(
                    // these normals make the planes face along z-axis
                    new genfit::DetPlane(TVector3(0,0,z)+o, u-o, v-o)
                )
            );
            mFSTSensorPlanes[is]->Print();

            //cout << "sensor = " << is << endl;
            //cout << "u vector" << endl;
            //u.Print();
            //cout << "v vector" << endl;
            //v.Print();
            //cout << "o vector" << endl;
            //o.Print();
               
            MhssOnFst.Zero();
            MfstWedgeOnHss.Zero();
            MfstSensorOnWedge.Zero();
        }

        //LOG_INFO  << sstr.str() << endm;

        // Now load FTT
        // mConfig.getVector<>(...) requires a default, hence the 
        mFTTZLocations = fwdGeoUtils.fttZ(
            mConfig.getVector<double>("TrackFitter.Geometry:ftt", {0.0f, 0.0f, 0.0f, 0.0f})
            );

        if ( fwdGeoUtils.fttZ( 0 ) < 1.0 ) { // returns 0.0 on failure
            LOG_WARN << "Using FTT z-locations from config or default, may not match hits" << endm;
        }

        if ( mFTTZLocations.size() != 4 ){
            LOG_ERROR << "Wrong number of FTT layers, got " << mFTTZLocations.size() << " but expected 4" << endm;
        }

        sstr.str("");
        sstr.clear();
        sstr << "Adding FTT Planes at: ";
        delim = "";
        for (auto z : mFTTZLocations) {
            mFTTPlanes.push_back(
                genfit::SharedPlanePtr(
                    // these normals make the planes face along z-axis
                    new genfit::DetPlane(TVector3(0, 0, z), TVector3(1, 0, 0), TVector3(0, 1, 0))
                )
            );
            sstr << delim << z;
            delim = ", ";
        }
        LOG_INFO << sstr.str() << endm;

        // get default vertex values used in simulation from the config
        mVertexSigmaXY = mConfig.get<double>("TrackFitter.Vertex:sigmaXY", 1.0);
        mVertexSigmaZ = mConfig.get<double>("TrackFitter.Vertex:sigmaZ", 30.0);
        mVertexPos = mConfig.getVector<double>("TrackFitter.Vertex:pos", {0.0,0.0,0.0});
        mIncludeVertexInFit = mConfig.get<bool>("TrackFitter.Vertex:includeInFit", false);
        mSmearMcVertex = mConfig.get<bool>("TrackFitter.Vertex:smearMcVertex", false);

        if ( mGenHistograms )
            makeHistograms();

        //delete tempMaker;
    }

    void finish() {
      mGblFitter->endRun();
    }

    void makeHistograms() {
        std::string n = "";
        mHist["ECalProjPosXY"] = new TH2F("ECalProjPosXY", ";X;Y", 1000, -500, 500, 1000, -500, 500);
        mHist["ECalProjSigmaXY"] = new TH2F("ECalProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 50, 0, 0.5, 50, 0, 0.5);
        mHist["ECalProjSigmaR"] = new TH1F("ECalProjSigmaR", ";#sigma_{XY} (cm) at ECAL", 50, 0, 0.5);

        mHist["SiProjPosXY"] = new TH2F("SiProjPosXY", ";X;Y", 1000, -500, 500, 1000, -500, 500);
        mHist["SiProjSigmaXY"] = new TH2F("SiProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 150, 0, 15, 150, 0, 15);

        mHist["VertexProjPosXY"] = new TH2F("VertexProjPosXY", ";X;Y", 100, -5, 5, 100, -5, 5);
        mHist["VertexProjSigmaXY"] = new TH2F("VertexProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 150, 0, 20, 150, 0, 20);

        mHist["VertexProjPosZ"] = new TH1F("VertexProjPosZ", ";Z;", 100, -50, 50);
        mHist["VertexProjSigmaZ"] = new TH1F("VertexProjSigmaZ", ";#sigma_{Z};", 100, 0, 10);

        mHist["SiWrongProjPosXY"] = new TH2F("SiWrongProjPosXY", ";X;Y", 1000, -500, 500, 1000, -500, 500);
        mHist["SiWrongProjSigmaXY"] = new TH2F("SiWrongProjSigmaXY", ";#sigma_{X};#sigma_{Y}", 50, 0, 0.5, 50, 0, 0.5);

        mHist["SiDeltaProjPosXY"] = new TH2F("SiDeltaProjPosXY", ";X;Y", 1000, 0, 20, 1000, 0, 20);

        mHist["FstDiffZVsR"] = new TH2F( "FstDiffZVsR", ";R;dz", 400, 0, 40, 500, -5, 5 );
        mHist["FstDiffZVsPhiSliceInner"] = new TH2F( "FstDiffZVsPhiSliceInner", ";slice;dz", 15, 0, 15, 500, -5, 5 );
        mHist["FstDiffZVsPhiSliceOuter"] = new TH2F( "FstDiffZVsPhiSliceOuter", ";slice;dz", 15, 0, 15, 500, -5, 5 );

        mHist["FstDiffZVsPhiOuter"] = new TH2F( "FstDiffZVsPhiOuter", ";slice;dz", 628, 0, TMath::Pi()*2, 500, -5, 5 );

        mHist["CorrFstDiffZVsPhiSliceInner"] = new TH2F( "CorrFstDiffZVsPhiSliceInner", ";slice;dz", 15, 0, 15, 500, -5, 5 );
        mHist["CorrFstDiffZVsPhiSliceOuter"] = new TH2F( "CorrFstDiffZVsPhiSliceOuter", ";slice;dz", 15, 0, 15, 500, -5, 5 );

        //for(int in = 1; in < 4; in++) 
        //{
        //  for(int id = 0; id < 3; id++)
        //  {
        //    for(int is = 0; is < 8; is++)
        //    {
        //      n = Form("TruthPhi_l%d_s%d_n%d",id,is,in);               mHist[n] = new TH1F(n.c_str(), "", 150, -TMath::Pi()/11.5, TMath::Pi()/11.5);
        //      n = Form("TruthR_l%d_s%d_n%d",id,is,in);                 mHist[n] = new TH1F(n.c_str(), "", 100, 5.0+is*2.875-0.25, 5.0+(is+1)*2.875+0.25);
        //      n = Form("MeasuredPhi_l%d_s%d_n%d",id,is,in);            mHist[n] = new TH1F(n.c_str(), "", 150, -TMath::Pi()/11.5, TMath::Pi()/11.5);
        //      n = Form("MeasuredR_l%d_s%d_n%d",id,is,in);              mHist[n] = new TH1F(n.c_str(), "", 100, 5.0+is*2.875-0.25, 5.0+(is+1)*2.875+0.25);
        //      n = Form("PredictedPhi_l%d_s%d_n%d",id,is,in);           mHist[n] = new TH1F(n.c_str(), "", 150, -TMath::Pi()/11.5, TMath::Pi()/11.5);
        //      n = Form("PredictedR_l%d_s%d_n%d",id,is,in);             mHist[n] = new TH1F(n.c_str(), "", 100, 5.0+is*2.875-0.25, 5.0+(is+1)*2.875+0.25);

        //      n = Form("TruthMinusMeasuredPhi_l%d_s%d_n%d",id,is,in);  mHist[n] = new TH1F(n.c_str(), "", 100, -0.004, 0.004);
        //      n = Form("TruthMinusMeasuredR_l%d_s%d_n%d",id,is,in);    mHist[n] = new TH1F(n.c_str(), "", 100, -2.0, 2.0);
        //      n = Form("TruthMinusMeasuredY_l%d_s%d_n%d",id,is,in);    mHist[n] = new TH1F(n.c_str(), "", 100, -0.4, .4);
        //      n = Form("TruthMinusMeasuredX_l%d_s%d_n%d",id,is,in);    mHist[n] = new TH1F(n.c_str(), "", 100, -2.0, 2.0);

        //      n = Form("PredictedMinusTruthPhi_l%d_s%d_n%d",id,is,in);  mHist[n] = new TH1F(n.c_str(), "", 100, -0.004, 0.004);
        //      n = Form("PredictedMinusTruthR_l%d_s%d_n%d",id,is,in);    mHist[n] = new TH1F(n.c_str(), "", 100, -2.0, 2.0);
        //      n = Form("PredictedMinusTruthY_l%d_s%d_n%d",id,is,in);    mHist[n] = new TH1F(n.c_str(), "", 100, -0.4, .4);
        //      n = Form("PredictedMinusTruthX_l%d_s%d_n%d",id,is,in);    mHist[n] = new TH1F(n.c_str(), "", 100, -2.0, 2.0);

        //      n = Form("PredictedMinusMeasuredPhi_l%d_s%d_n%d",id,is,in);  mHist[n] = new TH1F(n.c_str(), "", 100, -0.004, 0.004);
        //      n = Form("PredictedMinusMeasuredR_l%d_s%d_n%d",id,is,in);    mHist[n] = new TH1F(n.c_str(), "", 100, -2.0, 2.0);
        //      n = Form("PredictedMinusMeasuredY_l%d_s%d_n%d",id,is,in);    mHist[n] = new TH1F(n.c_str(), "", 100, -0.4, .4);
        //      n = Form("PredictedMinusMeasuredX_l%d_s%d_n%d",id,is,in);    mHist[n] = new TH1F(n.c_str(), "", 100, -2.0, 2.0);

        //      n = Form("PredictedMinusMeasuredVsTruthPhi_l%d_s%d_n%d",id,is,in);  mHist[n] = new TProfile(n.c_str(), "", 150, -TMath::Pi()/11., TMath::Pi()/11.);
        //      n = Form("PredictedMinusMeasuredVsTruthR_l%d_s%d_n%d",id,is,in);    mHist[n] = new TProfile(n.c_str(), "", 100, 5.0+is*2.875, 5.0+(is+1)*2.875);
        //      //n = Form("PredictedMinusMeasuredVsTruthY_l%d_s%d",id,is);    mHist[n] = new TProfile(n.c_str(), "", 100, -0.15, .15);
        //      //n = Form("PredictedMinusMeasuredVsTruthX_l%d_s%d",id,is);    mHist[n] = new TProfile(n.c_str(), "", 100, -2.0, 2.0);
  
        //      n = Form("PredictedMinusMeasuredVsTruthPhiStrip_l%d_s%d_n%d",id,is,in);  mHist[n] = new TProfile(n.c_str(), "", 100, -0.0025, 0.0025);

        //      n = Form("PredictedMinusMeasuredVsTruthPhi2D_l%d_s%d_n%d",id,is,in);  mHist2D[n] = new TH2F(n.c_str(), "", 150, -TMath::Pi()/11., TMath::Pi()/11., 100, -0.004, 0.004);
        //      n = Form("PredictedMinusMeasuredVsTruthR2D_l%d_s%d_n%d",id,is,in);    mHist2D[n] = new TH2F(n.c_str(), "", 100, 5.0+is*2.875, 5.0+(is+1)*2.875, 100, -2.0, 2.0);


        //      n = Form("PredictedMinusMeasuredVsTruthPhi2DStrip_l%d_s%d_n%d",id,is,in);  mHist2D[n] = new TH2F(n.c_str(), "", 100, -0.0025, 0.0025, 100, -0.004, 0.004);
        //      //n = Form("PredictedMinusMeasuredVsTruthY2D_l%d_s%d",id,is);    mHist2D[n] = new TH2F(n.c_str(), "", 100, -0.15, .15);
        //      //n = Form("PredictedMinusMeasuredVsTruthX2D_l%d_s%d",id,is);    mHist2D[n] = new TH2F(n.c_str(), "", 100, -2.0, 2.0);
        //    }
        //  }
        //}

        for(int id = 0; id < 3; id++)
        {
          for(int is = 0; is < 8; is++)
          {
            n = Form("TruthPhi_l%d_s%d",id,is);               mHist[n] = new TH1F(n.c_str(), "", 640, -TMath::Pi()/12., TMath::Pi()/12.);
            n = Form("TruthPhiStrip_l%d_s%d",id,is);          mHist[n] = new TH1F(n.c_str(), "", 100, -0.004, 0.004);
            n = Form("TruthR_l%d_s%d",id,is);                 mHist[n] = new TH1F(n.c_str(), "", 100, 5.0+is*2.875-0.25, 5.0+(is+1)*2.875+0.25);
            n = Form("MeasuredPhi_l%d_s%d",id,is);            mHist[n] = new TH1F(n.c_str(), "", 640, -TMath::Pi()/12., TMath::Pi()/12.);
            n = Form("MeasuredPhiStrip_l%d_s%d",id,is);       mHist[n] = new TH1F(n.c_str(), "", 100, -0.004, 0.004);
            n = Form("MeasuredR_l%d_s%d",id,is);              mHist[n] = new TH1F(n.c_str(), "", 100, 5.0+is*2.875-0.25, 5.0+(is+1)*2.875+0.25);
            n = Form("PredictedPhi_l%d_s%d",id,is);           mHist[n] = new TH1F(n.c_str(), "", 640, -TMath::Pi()/12., TMath::Pi()/12.);
            n = Form("PredictedPhiStrip_l%d_s%d",id,is);      mHist[n] = new TH1F(n.c_str(), "", 100, -0.004, 0.004);
            n = Form("PredictedR_l%d_s%d",id,is);             mHist[n] = new TH1F(n.c_str(), "", 100, 5.0+is*2.875-0.25, 5.0+(is+1)*2.875+0.25);

            n = Form("TruthMinusMeasuredPhi_l%d_s%d",id,is);  mHist[n] = new TH1F(n.c_str(), "", 100, -0.004, 0.004);
            n = Form("TruthMinusMeasuredR_l%d_s%d",id,is);    mHist[n] = new TH1F(n.c_str(), "", 100, -2.0, 2.0);
            n = Form("TruthMinusMeasuredY_l%d_s%d",id,is);    mHist[n] = new TH1F(n.c_str(), "", 100, -0.4, .4);
            n = Form("TruthMinusMeasuredX_l%d_s%d",id,is);    mHist[n] = new TH1F(n.c_str(), "", 100, -2.0, 2.0);

            n = Form("PredictedMinusTruthPhi_l%d_s%d",id,is);  mHist[n] = new TH1F(n.c_str(), "", 100, -0.004, 0.004);
            n = Form("PredictedMinusTruthR_l%d_s%d",id,is);    mHist[n] = new TH1F(n.c_str(), "", 100, -0.5, 0.5);
            n = Form("PredictedMinusTruthY_l%d_s%d",id,is);    mHist[n] = new TH1F(n.c_str(), "", 100, -0.08, 0.08);
            n = Form("PredictedMinusTruthX_l%d_s%d",id,is);    mHist[n] = new TH1F(n.c_str(), "", 100, -0.5, 0.5);

            n = Form("PredictedMinusMeasuredPhi_l%d_s%d",id,is);  mHist[n] = new TH1F(n.c_str(), "", 100, -0.004, 0.004);
            n = Form("PredictedMinusMeasuredR_l%d_s%d",id,is);    mHist[n] = new TH1F(n.c_str(), "", 100, -2.0, 2.0);
            n = Form("PredictedMinusMeasuredY_l%d_s%d",id,is);    mHist[n] = new TH1F(n.c_str(), "", 100, -0.4, .4);
            n = Form("PredictedMinusMeasuredX_l%d_s%d",id,is);    mHist[n] = new TH1F(n.c_str(), "", 100, -2.0, 2.0);

            n = Form("PredictedMinusMeasuredVsTruthPhi_l%d_s%d",id,is);  mHist[n] = new TProfile(n.c_str(), "", 256, -TMath::Pi()/12., TMath::Pi()/12.);
            n = Form("PredictedMinusMeasuredVsTruthR_l%d_s%d",id,is);    mHist[n] = new TProfile(n.c_str(), "", 100, 5.0+is*2.875, 5.0+(is+1)*2.875);
            //n = Form("PredictedMinusMeasuredVsTruthY_l%d_s%d",id,is);    mHist[n] = new TProfile(n.c_str(), "", 100, -0.15, .15);
            //n = Form("PredictedMinusMeasuredVsTruthX_l%d_s%d",id,is);    mHist[n] = new TProfile(n.c_str(), "", 100, -2.0, 2.0);
  
            n = Form("PredictedMinusMeasuredVsTruthPhiStrip_l%d_s%d",id,is);  mHist[n] = new TProfile(n.c_str(), "", 100, -0.0025, 0.0025);

            n = Form("PredictedMinusMeasuredVsTruthPhi2D_l%d_s%d",id,is);  mHist2D[n] = new TH2F(n.c_str(), "", 256, -TMath::Pi()/12., TMath::Pi()/12., 100, -0.004, 0.004);
            n = Form("PredictedMinusMeasuredVsTruthR2D_l%d_s%d",id,is);    mHist2D[n] = new TH2F(n.c_str(), "", 100, 5.0+is*2.875, 5.0+(is+1)*2.875, 100, -2.0, 2.0);

            n = Form("PredictedMinusMeasuredVsTruthPhi2DStrip_l%d_s%d",id,is);  mHist2D[n] = new TH2F(n.c_str(), "", 100, -0.0025, 0.0025, 100, -0.004, 0.004);

            //n = Form("PredictedMinusMeasuredVsTruthY2D_l%d_s%d",id,is);    mHist2D[n] = new TH2F(n.c_str(), "", 100, -0.15, .15);
            //n = Form("PredictedMinusMeasuredVsTruthX2D_l%d_s%d",id,is);    mHist2D[n] = new TH2F(n.c_str(), "", 100, -2.0, 2.0);
          }
        }



        n = "seed_curv";
        mHist[n] = new TH1F(n.c_str(), ";curv", 1000, 0, 10000);
        n = "seed_pT";
        mHist[n] = new TH1F(n.c_str(), ";pT (GeV/c)", 500, 0, 10);
        n = "seed_eta";
        mHist[n] = new TH1F(n.c_str(), ";eta", 500, 0, 5);

        n = "delta_fit_seed_pT";
        mHist[n] = new TH1F(n.c_str(), ";#Delta( fit, seed ) pT (GeV/c)", 500, -5, 5);
        n = "delta_fit_seed_eta";
        mHist[n] = new TH1F(n.c_str(), ";#Delta( fit, seed ) eta", 500, 0, 5);
        n = "delta_fit_seed_phi";
        mHist[n] = new TH1F(n.c_str(), ";#Delta( fit, seed ) phi", 500, -5, 5);

        n = "FitStatus";
        mHist[n] = new TH1F(n.c_str(), ";", 5, 0, 5);
        FwdTrackerUtils::labelAxis(mHist[n]->GetXaxis(), {"Total", "Pass", "Fail", "GoodCardinal", "Exception"});

        n = "FitDuration";
        mHist[n] = new TH1F(n.c_str(), "; Duraton (ms)", 5000, 0, 50000);

        n = "FailedFitDuration";
        mHist[n] = new TH1F(n.c_str(), "; Duraton (ms)", 500, 0, 50000);
    }

    // writes mHistograms stored in map only if mGenHistograms is true
    void writeHistograms() {
        if ( !mGenHistograms )
            return;
        for (auto nh : mHist) {
            nh.second->SetDirectory(gDirectory);
            nh.second->Write();
        }
        for (auto nh : mHist2D) {
            nh.second->SetDirectory(gDirectory);
            nh.second->Write();
        }
    }

    /* Convert the 3x3 covmat to 2x2 by dropping z
    *
    */
    TMatrixDSym CovMatPlane(KiTrack::IHit *h){
        TMatrixDSym cm(2);
        cm(0, 0) = static_cast<FwdHit*>(h)->_covmat(0, 0);
        cm(1, 1) = static_cast<FwdHit*>(h)->_covmat(1, 1);
        cm(0, 1) = static_cast<FwdHit*>(h)->_covmat(0, 1);
        cm(1, 0) = static_cast<FwdHit*>(h)->_covmat(1, 0);
        //cm.Print();
        return cm;
    }


    /* FitSimpleCircle
     * Used to determine a seed transverse momentum based on space points
     * Takes a list of space points KiTrack::IHit *
     * Takes three indecise used to lookup three of the possible hits within the list
     */ 
    float fitSimpleCircle(Seed_t trackCand, size_t i0, size_t i1, size_t i2) {
        float curv = 0;

        // ensure that no index is outside of range for FST or FTT volumes
        if (i0 > 12 || i1 > 12 || i2 > 12)
            return 0;

        try {
            KiTrack::SimpleCircle sc(trackCand[i0]->getX(), trackCand[i0]->getY(), trackCand[i1]->getX(), trackCand[i1]->getY(), trackCand[i2]->getX(), trackCand[i2]->getY());
            curv = sc.getRadius();
        } catch (KiTrack::InvalidParameter &e) {
            // if we got here we failed to get  a valid seed. We will still try to move forward but the fit will probably fail
        }

        //  make sure the curv is valid
        if (isinf(curv))
            curv = 999999.9;

        return curv;
    }

    /* seedState
     * Determines the seed position and momentum for a list of space points
     */
    float seedState(Seed_t trackCand, TVector3 &seedPos, TVector3 &seedMom) {
        // we require at least 4 hits,  so this should be gauranteed
        if(trackCand.size() < 3){
            // failure
            return 0.0;
        }
            

        // we want to use the LAST 3 hits, since silicon doesnt have R information
        TVector3 p0, p1, p2;
        // use the closest hit to the interaction point for the seed pos
        FwdHit *hit_closest_to_IP = static_cast<FwdHit *>(trackCand[0]);

        // maps from <key=vol_id> to <value=index in trackCand>
        std::map<size_t, size_t> vol_map; 

        // init the map
        for (size_t i = 0; i < 13; i++)
            vol_map[i] = -1;

        for (size_t i = 0; i < trackCand.size(); i++) {
            auto fwdHit = static_cast<FwdHit *>(trackCand[i]);
            vol_map[abs(fwdHit->_vid)] = i;
            LOG_INFO << "fwdHit->_vid = " << fwdHit->_vid << endm; 
            // find the hit closest to IP for the initial position seed
            if (hit_closest_to_IP->getZ() > fwdHit->getZ())
                hit_closest_to_IP = fwdHit;
        }

        // now get an estimate of the pT from several overlapping simple circle fits
        // enumerate the available partitions
        // 12 11 10
        // 12 11 9
        // 12 10 9
        // 11 10 9
        vector<float> curvs;
        curvs.push_back(fitSimpleCircle(trackCand, vol_map[12], vol_map[11], vol_map[10]));
        curvs.push_back(fitSimpleCircle(trackCand, vol_map[12], vol_map[11], vol_map[9]));
        curvs.push_back(fitSimpleCircle(trackCand, vol_map[12], vol_map[10], vol_map[9]));
        curvs.push_back(fitSimpleCircle(trackCand, vol_map[11], vol_map[10], vol_map[9]));

        // average them and exclude failed fits
        float mcurv = 0;
        float nmeas = 0;

        for (size_t i = 0; i < curvs.size(); i++) {
            if (mGenHistograms)
                this->mHist["seed_curv"]->Fill(curvs[i]);
            if (curvs[i] > 10) {
                mcurv += curvs[i];
                nmeas += 1.0;
            }
        }

        if (nmeas >= 1)
            mcurv = mcurv / nmeas;
        else
            mcurv = 10;

        // Now lets get eta information
        // simpler, use farthest points from IP
        if (vol_map[9] < 13)
            p0.SetXYZ(trackCand[vol_map[9]]->getX(), trackCand[vol_map[9]]->getY(), trackCand[vol_map[9]]->getZ());

        if (vol_map[10] < 13)
            p1.SetXYZ(trackCand[vol_map[10]]->getX(), trackCand[vol_map[10]]->getY(), trackCand[vol_map[10]]->getZ());

        const double K = 0.00029979; //K depends on the units used for Bfield
        double pt = mcurv * K * 5; // pT from average measured curv
        double dx = (p1.X() - p0.X());
        double dy = (p1.Y() - p0.Y());
        double dz = (p1.Z() - p0.Z());
        double phi = TMath::ATan2(dy, dx);
        double Rxy = sqrt(dx * dx + dy * dy);
        double theta = TMath::ATan2(Rxy, dz);
        // double eta = -log( tantheta / 2.0 );
        // these starting conditions can probably be improvd, good study for student

        seedMom.SetPtThetaPhi(pt, theta, phi);
        seedPos.SetXYZ(hit_closest_to_IP->getX(), hit_closest_to_IP->getY(), hit_closest_to_IP->getZ());

        if (mGenHistograms) {
            this->mHist["seed_pT"]->Fill(seedMom.Pt());
            this->mHist["seed_eta"]->Fill(seedMom.Eta());
        }

        return mcurv;
    }


    /*genfit::MeasuredStateOnPlane projectToFst(size_t si_plane, genfit::Track *fitTrack) {
        if (si_plane > 2) {
            genfit::MeasuredStateOnPlane nil;
            return nil;
        }

        auto detSi = mFSTSensorPlanes[si_plane];
        genfit::MeasuredStateOnPlane tst = fitTrack->getFittedState(1);
        auto TCM = fitTrack->getCardinalRep()->get6DCov(tst);
        //  can get the track length if needed
        // double len = fitTrack->getCardinalRep()->extrapolateToPlane(tst, detSi, false, true);
        double len = fitTrack->getCardinalRep()->extrapolateToPlane(tst, detSi);

        TCM = fitTrack->getCardinalRep()->get6DCov(tst);

        // can get the projected positions if needed
        float x = tst.getPos().X();
        float y = tst.getPos().Y();
        float z = tst.getPos().Z();
        // and the uncertainties
        LOG_INFO << "Track Uncertainty at FST (plane=" << si_plane << ") @ x= " << x << ", y= " << y << ", z= " << z << " : " << sqrt(TCM(0, 0)) << ", " << sqrt(TCM(1, 1)) << endm;

        return tst;
    }*/

    TVector2 projectToFst(size_t si_plane, genfit::Track *fitTrack) {
        if (si_plane > 107) {
            TVector2 nil;
            return nil;
        }

        auto detSi = mFSTSensorPlanes[si_plane];
        genfit::MeasuredStateOnPlane tst = fitTrack->getFittedState(1);
        auto TCM = fitTrack->getCardinalRep()->get6DCov(tst);

        //  can get the track length if needed
        // double len = fitTrack->getCardinalRep()->extrapolateToPlane(tst, detSi, false, true);
        double len = fitTrack->getCardinalRep()->extrapolateToPlane(tst, detSi);

        TCM = fitTrack->getCardinalRep()->get6DCov(tst);

        // can get the projected positions if needed
        double x = tst.getPos().X();
        double y = tst.getPos().Y();
        double z = tst.getPos().Z();
        // and the uncertainties
        LOG_INFO << "Track Uncertainty at FST (plane=" << si_plane << ") @ x= " << x << ", y= " << y << ", z= " << z << " : " << sqrt(TCM(0, 0)) << ", " << sqrt(TCM(1, 1)) << endm;

        //double vecPos[4] = {x,y,z,1.};
        //TMatrixD labPos(4,1,vecPos);
        //TMatrixD sensorPos = mInverseM[si_plane] * labPos;
        //return TVector2(sensorPos(0,0), sensorPos(1,0));
        return detSi->LabToPlane(tst.getPos());
    }

    genfit::SharedPlanePtr getFstPlane( FwdHit * h ){

        size_t planeId  = h->getSector();
        size_t sensorId = h->getSensor(); 
        
        int sensorIdx = sensorId;        

        if(!mMisaligned){
          //size_t sensorIdx = planeId*36 + moduleId*3 + sensorId;
          //LOG_INFO << "plane: " << planeId << "| module: " << moduleId << "| sensor: " << sensorId << "| sensorIdx: " << sensorIdx << "   PREVIOUS MODULE MAPPING" << endm;
          size_t moduleId = mModuleMap[planeId][(int(sensorId)/3)%12];
          sensorId = mSensorMap[sensorId%3];

          sensorIdx = planeId*36 + moduleId*3 + sensorId;
        }

     
        //size_t sensorIdx = planeId*36 + moduleId*3 + sensorI;

        //LOG_INFO << "plane: " << planeId << "| module: " << moduleId << "| sensor: " << sensorId << "| sensorIdx: " << sensorIdx << endm;
       
        auto planeCorr = mFSTSensorPlanes[sensorIdx];

        // Old method for determining the sensor z position
        //TVector3 hitXYZ( h->getX(), h->getY(), h->getZ() );

        //double phi = hitXYZ.Phi();
        //if ( phi < 0 ) phi = TMath::Pi() * 2 + phi;
        //const double phi_slice = phi / (TMath::Pi() / 6.0); // 2pi/12
        //const int phi_index = ((int)phi_slice);
        //const double r  =sqrt( pow(hitXYZ.x(), 2) + pow(hitXYZ.y(), 2) );

        //const size_t idx = phi_index % 2;
        //auto planeCorr = mFSTPlanesInner[planeId*2 + idx];
        //if ( r > 16 ){
        //    planeCorr = mFSTPlanesOuter[planeId*2 + idx];
        //}
        //double cdz = (h->getZ() - planeCorr->getO().Z());

        //if ( cdz > 0.010 ) {
        //    LOG_WARN << "FST Z =" << h->getZ() << " vs CORR Plane Z = " << planeCorr->getO().Z() << " DIFF: " << cdz << " phi_slice = " << phi_slice << ", phi_index = " << phi_index << " R=" << hitXYZ.Pt() << " idx=" << idx << endm;
        //}

        return planeCorr;

    }

    TVector3 refitTrackWithGBL( genfit::Track *originalTrack ) {
        // mem leak, global track is overwritten without delete.
        TVector3 pOrig = originalTrack->getCardinalRep()->getMom(originalTrack->getFittedState(1, originalTrack->getCardinalRep()));
        
        // auto cardinalStatus = originalTrack->getFitStatus(originalTrack->getCardinalRep());

        if (originalTrack->getFitStatus(originalTrack->getCardinalRep())->isFitConverged() == false) {
            // in this case the original track did not converge so we should not refit. 
            // probably never get here due to previous checks
            return pOrig;
        }

        // Setup the Track Reps
       // auto trackRepPos = new genfit::RKTrackRep(mPdgPositron);
       // auto trackRepNeg = new genfit::RKTrackRep(mPdgElectron);

        // get the space points on the original track
      //  auto trackPoints = originalTrack->getPointsWithMeasurement();
        

     //   TVectorD rawCoords = trackPoints[0]->getRawMeasurement()->getRawHitCoords();
     //   TVector3 seedPos(rawCoords(0), rawCoords(1), rawCoords(2));
     //   TVector3 seedMom = pOrig;

        // Create the ref track using the seed state
     //   auto pFitTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
     //   pFitTrack->addTrackRep(trackRepNeg);
        //cout << "Before copying track" << endl;
        //genfit::Track fitTrack = genfit::Track(*originalTrack);
        //cout << "After copying track" << endl;
        //for (size_t i = 0; i < trackPoints.size(); i++) {
            // clone the track points into this track
        //    fitTrack.insertPoint(new genfit::TrackPoint(trackPoints[i]->getRawMeasurement(), &fitTrack));
       // }

        //auto gblFitter = std::unique_ptr<genfit::GblFitter>(new genfit::GblFitter());
        try {
            // check consistency of all points
            //fitTrack.checkConsistency();
            originalTrack->checkConsistency();

            // do the actual track fit
            //cout << "Before fitting" << endl;
            mGblFitter->processTrackWithRep(originalTrack, originalTrack->getCardinalRep());
            //cout << "After fitting" << endl;
            //mGblFitter->processTrackWithRep(&fitTrack, trackRepNeg);
            //mFitter->processTrack(&fitTrack);

            //fitTrack.checkConsistency();
            originalTrack->checkConsistency();

            // this chooses the lowest chi2 fit result as cardinal
            //fitTrack.determineCardinalRep(); 
            originalTrack->determineCardinalRep(); 
           

        } catch (genfit::Exception &e) {
            // will be caught below by converge check
            LOG_WARN << "Track fit exception : " << e.what() << endm;
        }

        if (originalTrack->getFitStatus(originalTrack->getCardinalRep())->isFitConverged() == false) {
            LOG_WARN << "GBL fit did not converge" << endm;
            return pOrig;
        } else { // we did converge, return new momentum
            
            try {
                // causes seg fault
                auto cardinalRep = originalTrack->getCardinalRep();
                auto cardinalStatus = originalTrack->getFitStatus(cardinalRep);
                mFitStatus = *cardinalStatus; // save the status of last fit
            } catch (genfit::Exception &e) {
                LOG_WARN << "Failed to get cardinal status from converged fit" << endm;
            }

            return originalTrack->getCardinalRep()->getMom(originalTrack->getFittedState(1, originalTrack->getCardinalRep()));
        }
        return pOrig;
    } //refitwith GBL


    /* RefitTracksWithSiHits
     * Takes a previously fit track re-fits it with the newly added silicon hits 
     * 
     */
    TVector3 refitTrackWithSiHits(genfit::Track *originalTrack, Seed_t si_hits) {
        // mem leak, global track is overwritten without delete.
        TVector3 pOrig = originalTrack->getCardinalRep()->getMom(originalTrack->getFittedState(1, originalTrack->getCardinalRep()));
        
        // auto cardinalStatus = originalTrack->getFitStatus(originalTrack->getCardinalRep());

        if (originalTrack->getFitStatus(originalTrack->getCardinalRep())->isFitConverged() == false) {
            // in this case the original track did not converge so we should not refit. 
            // probably never get here due to previous checks
            return pOrig;
        }

        // Setup the Track Reps
        auto trackRepPos = new genfit::RKTrackRep(mPdgPositron);
        auto trackRepNeg = new genfit::RKTrackRep(mPdgElectron);

        // get the space points on the original track
        auto trackPoints = originalTrack->getPointsWithMeasurement();
        
        if ((trackPoints.size() < (mFTTZLocations.size() +1) && mIncludeVertexInFit) || trackPoints.size() < mFTTZLocations.size() ) {
            // we didnt get enough points for a refit
            return pOrig;
        }

        TVectorD rawCoords = trackPoints[0]->getRawMeasurement()->getRawHitCoords();
        double z = mFSTZLocations[0]; //first FTT plane, used if we dont have PV in fit
        if (mIncludeVertexInFit)
            z = rawCoords(2);

        TVector3 seedPos(rawCoords(0), rawCoords(1), z);
        TVector3 seedMom = pOrig;

        // Create the ref track using the seed state
        auto pFitTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
        pFitTrack->addTrackRep(trackRepNeg);

        genfit::Track &fitTrack = *pFitTrack;

        size_t firstFTTIndex = 0;
        if (mIncludeVertexInFit) {
            // clone the PRIMARY VERTEX into this track
            fitTrack.insertPoint(new genfit::TrackPoint(trackPoints[0]->getRawMeasurement(), &fitTrack));
            firstFTTIndex = 1; // start on hit index 1 below
        }

        // initialize the hit coords on plane
        TVectorD hitCoords(2);
        hitCoords[0] = 0;
        hitCoords[1] = 0;

        size_t planeId(0);
        int hitId(5);

        //LOG_INFO << "Start adding silicon hits to track after cloning" << endm;

        // add the hits to the track
        for (auto h : si_hits) {
            if ( nullptr == h ) continue; // if no Si hit in this plane, skip

            hitCoords[0] = h->getX();
            hitCoords[1] = h->getY();
          
            //LOG_INFO << "Hit position: X = " << hitCoords[0] << ", Y = " << hitCoords[1] << endm;

            planeId = h->getSector();
            //LOG_INFO << "Plane ID: " << planeId << endm;
            auto plane = getFstPlane( static_cast<FwdHit*>(h) );
            int sensorId = static_cast<FwdHit*>(h)->getSensor();

            genfit::PlanarMeasurement *measurement = new genfit::PlanarMeasurement(hitCoords, CovMatPlane(h), sensorId+1000, ++hitId, nullptr);

            if (mFSTPlanes.size() <= planeId) {
                LOG_WARN << "invalid VolumId -> out of bounds DetPlane, vid = " << planeId << endm;
                return pOrig;
            }

            //LOG_INFO << "Created the genfit::PlanarMeasurement" << endm;
            // auto plane = mFSTPlanes[planeId];
            //LOG_INFO << "Sensor ID: " << sensorId << endm;           

            //measurement->setPlane(plane, planeId);
            measurement->setPlane(plane, sensorId+1000);
            fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

            //TVector3 hitXYZ( h->getX(), h->getY(), h->getZ() );
            //float phi = hitXYZ.Phi();
            //if ( phi < 0 ) phi = TMath::Pi() * 2 + phi;
            //double phi_slice = phi / (TMath::Pi() / 6.0); // 2pi/12
            //int phi_index = ((int)phi_slice);
            //double dz = (h->getZ() - plane->getO().Z());

            //double r  =sqrt( pow(hitXYZ.x(), 2) + pow(hitXYZ.y(), 2) );

            //size_t idx = phi_index % 2;
            //auto planeCorr = mFSTPlanesInner[planeId + idx];
            //if ( r > 16 ){
            //    planeCorr = mFSTPlanesOuter[planeId + idx];
            //}
            //double cdz = (h->getZ() - planeCorr->getO().Z());

            //if (mGenHistograms){
            //    
            //    ((TH2*)mHist[ "FstDiffZVsR" ])->Fill( r, dz );

            //    if ( r < 16 ) {// inner
            //        mHist["FstDiffZVsPhiSliceInner"]->Fill( phi_slice, dz );
            //        mHist["CorrFstDiffZVsPhiSliceInner"]->Fill( phi_slice, cdz );
            //    } else {
            //        mHist["FstDiffZVsPhiSliceOuter"]->Fill( phi_slice, dz );
            //        mHist["CorrFstDiffZVsPhiSliceOuter"]->Fill( phi_slice, cdz );
            //        mHist["FstDiffZVsPhiOuter"]->Fill( phi, dz );
            //    }
            //}
            // mHist[ "FstDiffZVsPhiSliceInner" ]->Fill( sqrt( pow(hitXYZ.x(), 2), pow(hitXYZ.y(), 2) ), dz );

        }
        //LOG_INFO << "Finish adding silicon hits to track after cloning" << endm;
        // start at 0 if PV not included, 1 otherwise 
        
        //LOG_INFO << "Start adding old vertex and FTT hits" << endm;
        for (size_t i = firstFTTIndex; i < trackPoints.size(); i++) {
            // clone the track points into this track
            fitTrack.insertPoint(new genfit::TrackPoint(trackPoints[i]->getRawMeasurement(), &fitTrack));
        }
        //LOG_INFO << "Finish adding old vertex and FTT hits" << endm;

        //LOG_INFO << "Try refitting track with all hits present" << endm;
        try {
            //Track RE-Fit with GENFIT2
            // check consistency of all points
            
            //LOG_INFO << "fitTrack.checkConsistency();" << endm;
            fitTrack.checkConsistency();

            // do the actual track fit
            //LOG_INFO << "Process Track" << endm;
            mFitter->processTrack(&fitTrack);

            //LOG_INFO << "fitTrack.checkConsistency();" << endm;
            fitTrack.checkConsistency();

            // this chooses the lowest chi2 fit result as cardinal
            //LOG_INFO << "fitTrack.determineCardinalRep();" << endm;
            fitTrack.determineCardinalRep(); 

        } catch (genfit::Exception &e) {
            // will be caught below by converge check
            LOG_WARN << "Track fit exception : " << e.what() << endm;
        }

        if (fitTrack.getFitStatus(fitTrack.getCardinalRep())->isFitConverged() == false) {
            // Did not converge
            return pOrig;
        } else { // we did converge, return new momentum
            
            try {
                // causes seg fault
                auto cardinalRep = fitTrack.getCardinalRep();
                auto cardinalStatus = fitTrack.getFitStatus(cardinalRep);
                mFitStatus = *cardinalStatus; // save the status of last fit
            } catch (genfit::Exception &e) {
            }

            TVector3 p = fitTrack.getCardinalRep()->getMom(fitTrack.getFittedState(1, fitTrack.getCardinalRep()));
            // get status if needed later
            // auto newStatus = fitTrack.getFitStatus(fitTrack.getCardinalRep());

            // try {
            //     LOG_INFO << "Rechecking projected errors" << endm;
            //     auto amsp2 = projectToFst(2, pFitTrack);
            //     auto amsp1 = projectToFst(1, pFitTrack);
            //     auto amsp0 = projectToFst(0, pFitTrack);
            // }catch (genfit::Exception &e) {
            // }
            //
            

            if(mRefitGBL) 
            {
              mGblFitter->setSuccessfulFitFlag(false);
              p = refitTrackWithGBL(&fitTrack);
            }

            if(mGblFitter->getSuccessfulFitFlag())
            {
              int nhits = 0;
              for (auto h : si_hits) 
              {
                if ( nullptr == h ) continue; // if no Si hit in this plane, skip
                nhits++;
              }
              for (auto h : si_hits) 
              {
                if ( nullptr == h ) continue; // if no Si hit in this plane, skip

                int disk = h->getSector();
                int sid = static_cast<FwdHit*>(h)->getSensor();
                cout << "sid = " << sid << endl;

                auto plane = getFstPlane( static_cast<FwdHit*>(h) );
                TVector3 origin = plane->getO();
                double mcZ;
                double rcZ = origin.Z();

                TVectorD rcCoords(2);
                TVectorD mcCoords(2);
                TVectorD trackCoords(2);

                rcCoords[0] = h->getX();
                rcCoords[1] = h->getY();
              
        
                for (auto mch : static_cast<FwdHit*>(h)->_mcTrack->mHits)
                {
                  if(mch->getZ() > 200.0) continue;
                  if( mch->getSector() == disk) 
                  {
                    double arr[4] = {mch->getX(),mch->getY(),mch->getZ(),1.0};
                    TMatrixD mc4D(4,1,arr);
                    mc4D = mInverseM[sid] * mc4D;
                   
                    mcCoords[0] = mc4D(0,0);
                    mcCoords[1] = mc4D(1,0);
                    mcZ = mc4D(2,0);
                  }
                }      
              
                cout << "mcZ = " << mcZ << "    rcZ = " << rcZ << "       mcZ-rcZ = " << mcZ-rcZ << endl;

                for( int tp = 0; tp < fitTrack.getNumPointsWithMeasurement(); tp++)
                {
                  genfit::TrackPoint* point_meas_temp = fitTrack.getPointWithMeasurement(tp);
                  genfit::PlanarMeasurement* measPlanar = dynamic_cast<genfit::PlanarMeasurement*>(point_meas_temp->getRawMeasurement(0));
                  int stateid;
                  if (measPlanar) stateid = measPlanar->getPlaneId();
                  if (stateid == sid + 1000)
                  {
                    genfit::KalmanFitterInfo* fi = dynamic_cast<genfit::KalmanFitterInfo*>(point_meas_temp->getFitterInfo(fitTrack.getCardinalRep()));
                    genfit::ReferenceStateOnPlane* reference = new genfit::ReferenceStateOnPlane(*fi->getReferenceState());
                    TVectorD state = reference->getState();
                    genfit::AbsMeasurement* raw_meas = point_meas_temp->getRawMeasurement(0);
                    std::unique_ptr<const genfit::AbsHMatrix> HitHMatrix(raw_meas->constructHMatrix(fitTrack.getCardinalRep()));
                    TVectorD planeState(HitHMatrix->Hv(state));
                
                    trackCoords[0] = planeState(0);
                    trackCoords[1] = planeState(1);
                  }
                }

                double mcR = TMath::Sqrt(mcCoords[0]*mcCoords[0]+mcCoords[1]*mcCoords[1]);
                double mcP = TMath::ATan2(mcCoords[1],mcCoords[0]);
                double rcR = TMath::Sqrt(rcCoords[0]*rcCoords[0]+rcCoords[1]*rcCoords[1]);
                double rcP = TMath::ATan2(rcCoords[1],rcCoords[0]);
                double trackR = TMath::Sqrt(trackCoords[0]*trackCoords[0]+trackCoords[1]*trackCoords[1]);
                double trackP = TMath::ATan2(trackCoords[1],trackCoords[0]);

                double mcPhiStrip;
                int    mcPhiStripInt;
                double mcPhiStripPos;

                double rcPhiStrip;
                int    rcPhiStripInt;
                double rcPhiStripPos;

                double trackPhiStrip;
                int    trackPhiStripInt;
                double trackPhiStripPos;

                double stripWidth = TMath::Pi() / 6.0 / 128. ;
                if(sid % 3 == 0) // inner sensor
                {
                  if(mcP >= 0.0) 
                  {  
                    mcPhiStrip = mcP / stripWidth ; 
                    rcPhiStrip = rcP / stripWidth ; 
                    trackPhiStrip = trackP / stripWidth ; 
                  }
                  if(mcP < 0.0) 
                  {
                    mcPhiStrip = (mcP+2*TMath::Pi()) / stripWidth ; 
                    rcPhiStrip = (rcP+2*TMath::Pi()) / stripWidth ; 
                    trackPhiStrip = (trackP+2*TMath::Pi()) / stripWidth ; 
                  }
                  mcPhiStripInt = int(mcPhiStrip) ; 
                  mcPhiStripPos = (mcPhiStrip - double(mcPhiStripInt) - 0.5) * stripWidth;
                  cout << "mcP = " << mcP << "  stripWidth = " << stripWidth << "   mcPhiStrip = " << mcPhiStrip << "   mcPhiStripInt " << mcPhiStripInt << "   mcPhiStripPos = " << mcPhiStripPos << endl;

                  rcPhiStripInt = int(rcPhiStrip) ; 
                  rcPhiStripPos = (rcPhiStrip - double(rcPhiStripInt) - 0.5) * stripWidth;
                  cout << "rcP = " << rcP << "  stripWidth = " << stripWidth << "   rcPhiStrip = " << rcPhiStrip << "   rcPhiStripInt " << rcPhiStripInt << "   rcPhiStripPos = " << rcPhiStripPos << endl;

                  trackPhiStripInt = int(trackPhiStrip) ; 
                  trackPhiStripPos = (trackPhiStrip - double(rcPhiStripInt) - 0.5) * stripWidth;
                  cout << "trackP = " << trackP << "  stripWidth = " << stripWidth << "   trackPhiStrip = " << trackPhiStrip << "   trackPhiStripInt " << trackPhiStripInt << "   trackPhiStripPos = " << trackPhiStripPos << endl;
                  
                }

                int is = -1;
                for(int i = 0; i < 4; i++)
                {
                  if(rcR >= 5.0 + i * 2.875 && rcR < 5.0 + (i + 1) * 2.875 && sid % 3 == 0 ) is = i;
                  if(rcR >= 5.0 + (i + 4) * 2.875 && rcR < 5.0 + (i + 5) * 2.875 && sid % 3 != 0 ) is = i + 4;
                }
                if(is < 0) continue;

                std::string n;

                n = Form("TruthPhi_l%d_s%d",disk,is);                                mHist[n]->Fill(mcP); 
                n = Form("TruthPhiStrip_l%d_s%d",disk,is);                           mHist[n]->Fill(mcPhiStripPos); 
                n = Form("TruthR_l%d_s%d",disk,is);                                  mHist[n]->Fill(mcR);
                n = Form("MeasuredPhi_l%d_s%d",disk,is);                             mHist[n]->Fill(rcP);
                n = Form("MeasuredPhiStrip_l%d_s%d",disk,is);                        mHist[n]->Fill(rcPhiStripPos);
                n = Form("MeasuredR_l%d_s%d",disk,is);                               mHist[n]->Fill(rcR);
                n = Form("PredictedPhi_l%d_s%d",disk,is);                            mHist[n]->Fill(trackP);
                n = Form("PredictedPhiStrip_l%d_s%d",disk,is);                       mHist[n]->Fill(trackPhiStripPos);
                n = Form("PredictedR_l%d_s%d",disk,is);                              mHist[n]->Fill(trackR);

                n = Form("TruthMinusMeasuredPhi_l%d_s%d",disk,is);                 mHist[n]->Fill(mcP-rcP                );
                n = Form("TruthMinusMeasuredR_l%d_s%d",disk,is);                   mHist[n]->Fill(mcR-rcR                );
                n = Form("TruthMinusMeasuredY_l%d_s%d",disk,is);                   mHist[n]->Fill(mcCoords[1]-rcCoords[1]);
                n = Form("TruthMinusMeasuredX_l%d_s%d",disk,is);                   mHist[n]->Fill(mcCoords[0]-rcCoords[0]);

                n = Form("PredictedMinusTruthPhi_l%d_s%d",disk,is);                mHist[n]->Fill(trackP-mcP                );
                n = Form("PredictedMinusTruthR_l%d_s%d",disk,is);                  mHist[n]->Fill(trackR-mcR                );
                n = Form("PredictedMinusTruthY_l%d_s%d",disk,is);                  mHist[n]->Fill(trackCoords[1]-mcCoords[1]);
                n = Form("PredictedMinusTruthX_l%d_s%d",disk,is);                  mHist[n]->Fill(trackCoords[0]-mcCoords[0]);

                n = Form("PredictedMinusMeasuredPhi_l%d_s%d",disk,is);             mHist[n]->Fill(trackP-rcP                );
                n = Form("PredictedMinusMeasuredR_l%d_s%d",disk,is);               mHist[n]->Fill(trackR-rcR                );
                n = Form("PredictedMinusMeasuredY_l%d_s%d",disk,is);               mHist[n]->Fill(trackCoords[1]-rcCoords[1]);
                n = Form("PredictedMinusMeasuredX_l%d_s%d",disk,is);               mHist[n]->Fill(trackCoords[0]-rcCoords[0]);

                n = Form("PredictedMinusMeasuredVsTruthPhi_l%d_s%d",disk,is);      mHist[n]->Fill(mcP,         trackP-rcP                );
                n = Form("PredictedMinusMeasuredVsTruthR_l%d_s%d",disk,is);        mHist[n]->Fill(mcR,         trackR-rcR                );
                //n = Form("PredictedMinusMeasuredVsTruthY_l%d_s%d",disk,is);        mHist[n]->Fill(mcCoords[1], trackCoords[1]-rcCoords[1]);
                //n = Form("PredictedMinusMeasuredVsTruthX_l%d_s%d",disk,is);        mHist[n]->Fill(mcCoords[0], trackCoords[0]-rcCoords[0]);

                n = Form("PredictedMinusMeasuredVsTruthPhiStrip_l%d_s%d",disk,is);      mHist[n]->Fill(mcPhiStripPos,         trackP-rcP                );

                n = Form("PredictedMinusMeasuredVsTruthPhi2D_l%d_s%d",disk,is);  mHist2D[n]->Fill(mcP,         trackP-rcP                );
                n = Form("PredictedMinusMeasuredVsTruthR2D_l%d_s%d",disk,is);    mHist2D[n]->Fill(mcR,         trackR-rcR                );

               
                n = Form("PredictedMinusMeasuredVsTruthPhi2DStrip_l%d_s%d",disk,is);      mHist2D[n]->Fill(mcPhiStripPos,         trackP-rcP                );
                //n = Form("PredictedMinusMeasuredVsTruthY2D_l%d_s%d",disk,is);    mHist2D[n]->Fill(mcCoords[1], trackCoords[1]-rcCoords[1]);
                //n = Form("PredictedMinusMeasuredVsTruthX2D_l%d_s%d",disk,is);    mHist2D[n]->Fill(mcCoords[0], trackCoords[0]-rcCoords[0]);



                //n = Form("TruthPhi_l%d_s%d_n%d",disk,is,nhits);                                   mHist[n]->Fill(mcP); 
                //n = Form("TruthR_l%d_s%d_n%d",disk,is,nhits);                                     mHist[n]->Fill(mcR);
                //n = Form("MeasuredPhi_l%d_s%d_n%d",disk,is,nhits);                                mHist[n]->Fill(rcP);
                //n = Form("MeasuredR_l%d_s%d_n%d",disk,is,nhits);                                  mHist[n]->Fill(rcR);
                //n = Form("PredictedPhi_l%d_s%d_n%d",disk,is,nhits);                               mHist[n]->Fill(trackP);
                //n = Form("PredictedR_l%d_s%d_n%d",disk,is,nhits);                                 mHist[n]->Fill(trackR);

                //n = Form("TruthMinusMeasuredPhi_l%d_s%d_n%d",disk,is,nhits);                 mHist[n]->Fill(mcP-rcP                );
                //n = Form("TruthMinusMeasuredR_l%d_s%d_n%d",disk,is,nhits);                   mHist[n]->Fill(mcR-rcR                );
                //n = Form("TruthMinusMeasuredY_l%d_s%d_n%d",disk,is,nhits);                   mHist[n]->Fill(mcCoords[1]-rcCoords[1]);
                //n = Form("TruthMinusMeasuredX_l%d_s%d_n%d",disk,is,nhits);                   mHist[n]->Fill(mcCoords[0]-rcCoords[0]);

                //n = Form("PredictedMinusTruthPhi_l%d_s%d_n%d",disk,is,nhits);                mHist[n]->Fill(trackP-mcP                );
                //n = Form("PredictedMinusTruthR_l%d_s%d_n%d",disk,is,nhits);                  mHist[n]->Fill(trackR-mcR                );
                //n = Form("PredictedMinusTruthY_l%d_s%d_n%d",disk,is,nhits);                  mHist[n]->Fill(trackCoords[1]-mcCoords[1]);
                //n = Form("PredictedMinusTruthX_l%d_s%d_n%d",disk,is,nhits);                  mHist[n]->Fill(trackCoords[0]-mcCoords[0]);

                //n = Form("PredictedMinusMeasuredPhi_l%d_s%d_n%d",disk,is,nhits);             mHist[n]->Fill(trackP-rcP                );
                //n = Form("PredictedMinusMeasuredR_l%d_s%d_n%d",disk,is,nhits);               mHist[n]->Fill(trackR-rcR                );
                //n = Form("PredictedMinusMeasuredY_l%d_s%d_n%d",disk,is,nhits);               mHist[n]->Fill(trackCoords[1]-rcCoords[1]);
                //n = Form("PredictedMinusMeasuredX_l%d_s%d_n%d",disk,is,nhits);               mHist[n]->Fill(trackCoords[0]-rcCoords[0]);

                //n = Form("PredictedMinusMeasuredVsTruthPhi_l%d_s%d_n%d",disk,is,nhits);      mHist[n]->Fill(mcP,         trackP-rcP                );
                //n = Form("PredictedMinusMeasuredVsTruthR_l%d_s%d_n%d",disk,is,nhits);        mHist[n]->Fill(mcR,         trackR-rcR                );
                ////n = Form("PredictedMinusMeasuredVsTruthY_l%d_s%d",disk,is);        mHist[n]->Fill(mcCoords[1], trackCoords[1]-rcCoords[1]);
                ////n = Form("PredictedMinusMeasuredVsTruthX_l%d_s%d",disk,is);        mHist[n]->Fill(mcCoords[0], trackCoords[0]-rcCoords[0]);

                //n = Form("PredictedMinusMeasuredVsTruthPhiStrip_l%d_s%d_n%d",disk,is,nhits);      mHist[n]->Fill(mcPhiStripPos,         trackP-rcP                );

                //n = Form("PredictedMinusMeasuredVsTruthPhi2D_l%d_s%d_n%d",disk,is,nhits);  mHist2D[n]->Fill(mcP,         trackP-rcP                );
                //n = Form("PredictedMinusMeasuredVsTruthR2D_l%d_s%d_n%d",disk,is,nhits);    mHist2D[n]->Fill(mcR,         trackR-rcR                );

                //n = Form("PredictedMinusMeasuredVsTruthPhi2DStrip_l%d_s%d_n%d",disk,is,nhits);      mHist2D[n]->Fill(mcPhiStripPos,         trackP-rcP                );
                //n = Form("PredictedMinusMeasuredVsTruthY2D_l%d_s%d",disk,is);    mHist2D[n]->Fill(mcCoords[1], trackCoords[1]-rcCoords[1]);
                //n = Form("PredictedMinusMeasuredVsTruthX2D_l%d_s%d",disk,is);    mHist2D[n]->Fill(mcCoords[0], trackCoords[0]-rcCoords[0]);
              } 
            }

            return p;
        }
        return pOrig;
    } // refit with Si hits



    /* Generic method for fitting space points with GenFit
     *
     * 
     */
    TVector3 fitSpacePoints( vector<genfit::SpacepointMeasurement*> spoints, TVector3 &seedPos, TVector3 &seedMom ){
        
        // setup track reps
        auto trackRepPos = new genfit::RKTrackRep(mPdgPositron);
        auto trackRepNeg = new genfit::RKTrackRep(mPdgElectron);

        // setup track for fit with positive and negative reps
        auto mFitTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
        mFitTrack->addTrackRep(trackRepNeg);

        genfit::Track &fitTrack = *mFitTrack;

        // try adding the points to track and fitting
        try {
            for ( size_t i = 0; i < spoints.size(); i++ ){
                fitTrack.insertPoint(new genfit::TrackPoint(spoints[i], &fitTrack));
            }
            // do the fit against the two possible fits
            mFitter->processTrackWithRep(&fitTrack, trackRepPos);
            mFitter->processTrackWithRep(&fitTrack, trackRepNeg);

        } catch (genfit::Exception &e) {
            LOG_ERROR << "GenFit failed to fit track with: " << e.what() << endm;
        }

        try {
            fitTrack.checkConsistency();

            fitTrack.determineCardinalRep();
            auto cardinalRep = fitTrack.getCardinalRep();

            TVector3 p = cardinalRep->getMom(fitTrack.getFittedState(1, cardinalRep));
            // sucess, return momentum
            return p;
        } catch (genfit::Exception &e) {
            LOG_ERROR << "GenFit failed to fit track with: " << e.what() << endm;
        }
        return TVector3(0, 0, 0);
    }

    /* Fit a track 
     *
     * 
     */
    TVector3 fitTrack(Seed_t trackCand, double *Vertex = 0, TVector3 *McSeedMom = 0) {
        long long itStart = FwdTrackerUtils::nowNanoSecond();
        if (mGenHistograms) this->mHist["FitStatus"]->Fill("Total", 1);

        // The PV information, if we want to use it
        TVectorD pv(3);

        StarRandom rand = StarRandom::Instance();
        if (mSmearMcVertex/*0 == Vertex*/) { // randomized from simulation
            pv[0] = mVertexPos[0] + rand.gauss(mVertexSigmaXY);
            pv[1] = mVertexPos[1] + rand.gauss(mVertexSigmaXY);
            pv[2] = mVertexPos[2] + rand.gauss(mVertexSigmaZ);
        } else {
            pv[0] = Vertex[0];
            pv[1] = Vertex[1];
            pv[2] = Vertex[2];
        }
        ///cout << "vertex = " << pv[0] << "," << pv[1] << "," << pv[2] << endl;
        //pv[0] = 0.0;
        //pv[1] = 0.0;
        //pv[2] = 0.0;
        //cout << "vertex = " << pv[0] << "," << pv[1] << "," << pv[2] << endl;

        // get the seed info from our hits
        TVector3 seedMom, seedPos;
        float curv = seedState(trackCand, seedPos, seedMom);

        if (McSeedMom != nullptr) {
            seedMom = *McSeedMom;
        }

        // If we use the PV, use that as the start pos for the track
        if (mIncludeVertexInFit) {
            LOG_INFO << "Primary Vertex in fit (seed pos) @ " << TString::Format( "(%f, %f, %f)", pv[0], pv[1], pv[2] ).Data()  << endm;
            seedPos.SetXYZ(pv[0], pv[1], pv[2]);
        }

        if (mFitTrack){
            delete mFitTrack;
        }

        // create the track representations
        auto trackRepPos = new genfit::RKTrackRep(mPdgPositron);
        auto trackRepNeg = new genfit::RKTrackRep(mPdgElectron);

        // Create the track
        mFitTrack = new genfit::Track(trackRepPos, seedPos, seedMom);
        mFitTrack->addTrackRep(trackRepNeg);


        LOG_INFO
            << "seedPos : (" << seedPos.X() << ", " << seedPos.Y() << ", " << seedPos.Z() << " )"
            << ", seedMom : (" << seedMom.X() << ", " << seedMom.Y() << ", " << seedMom.Z() << " )"
            << ", seedMom : (" << seedMom.Pt() << ", " << seedMom.Eta() << ", " << seedMom.Phi() << " )"
            << endm;

        genfit::Track &fitTrack = *mFitTrack;

        size_t planeId(0);     // detector plane ID
        int hitId(0);       // hit ID

        // initialize the hit coords on plane
        TVectorD hitCoords(2);
        hitCoords[0] = 0;
        hitCoords[1] = 0;

        /******************************************************************************************************************
        * Include the Primary vertex if desired
        ******************************************************************************************************************/
        if (mIncludeVertexInFit) {

            TMatrixDSym hitCov3(3);
            hitCov3(0, 0) = mVertexSigmaXY * mVertexSigmaXY;
            hitCov3(1, 1) = mVertexSigmaXY * mVertexSigmaXY;
            hitCov3(2, 2) = mVertexSigmaZ * mVertexSigmaZ;
            hitCov3.Print();

            //LOG_INFO << "create vertex spacepoint measurement" << endm;
            genfit::SpacepointMeasurement *measurement = new genfit::SpacepointMeasurement(pv, hitCov3, 9999, ++hitId, nullptr);
            //LOG_INFO << "insert point" << endm;
            fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
        }
        //LOG_INFO << "mIncludeVertexInFit" << endm;

        /******************************************************************************************************************
		 * loop over the hits, add them to the track
		 ******************************************************************************************************************/
        for(int ip = 0; ip < 4; ip++) {
            for (auto h : trackCand) {
                planeId = h->getSector();
                if(planeId != ip) continue;
                cout << "FTT plane " << planeId << "added" << endl;
               
                hitCoords[0] = h->getX();
                hitCoords[1] = h->getY();
                
                genfit::PlanarMeasurement *measurement = new genfit::PlanarMeasurement(hitCoords, CovMatPlane(h), h->getSector(), ++hitId, nullptr);


                if (mFTTPlanes.size() <= planeId) {
                    LOG_WARN << "invalid VolumId -> out of bounds DetPlane, vid = " << planeId << endm;
                    return TVector3(0, 0, 0);
                }

                auto plane = mFTTPlanes[planeId];
                measurement->setPlane(plane, planeId);
                fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));

                if (abs(h->getZ() - plane->getO().Z()) > 0.05) {
                    LOG_WARN << "Z Mismatch h->z = " << h->getZ() << ", plane->z = "<< plane->getO().Z() <<", diff = " << abs(h->getZ() - plane->getO().Z()) << endm;
                }
            }
        } // loop on trackCand
        //LOG_INFO << "after loop on trackCand" << endm;

        /******************************************************************************************************************
		 * Do the fit
		 ******************************************************************************************************************/
        try {
            // do the fit
            mFitter->processTrackWithRep(&fitTrack, trackRepPos);
            mFitter->processTrackWithRep(&fitTrack, trackRepNeg);
            //LOG_INFO << "after doing the fit" << endm;
        } catch (genfit::Exception &e) {
            if (mGenHistograms) mHist["FitStatus"]->Fill("Exception", 1);
        }

        TVector3 p(0, 0, 0);

        /******************************************************************************************************************
		 * Now check the fit
		 ******************************************************************************************************************/
        try {
            //check
            fitTrack.checkConsistency();

            // find track rep with smallest chi2
            fitTrack.determineCardinalRep();
            auto cardinalRep = fitTrack.getCardinalRep();
            auto cardinalStatus = fitTrack.getFitStatus(cardinalRep);
            mFitStatus = *cardinalStatus; // save the status of last fit

            // Delete any previous track rep
            if (mTrackRep)
                delete mTrackRep;

            // Clone the cardinal rep for persistency
            mTrackRep = cardinalRep->clone(); // save the result of the fit
            if (fitTrack.getFitStatus(cardinalRep)->isFitConverged()) {
                LOG_INFO << "Track Fit converged" << endm;
            }
            if (fitTrack.getFitStatus(cardinalRep)->isFitConverged() && mGenHistograms ) {
                this->mHist["FitStatus"]->Fill("GoodCardinal", 1);
            }

            if (fitTrack.getFitStatus(trackRepPos)->isFitConverged() == false &&
                fitTrack.getFitStatus(trackRepNeg)->isFitConverged() == false) {
            
                LOG_WARN << "Fit Failed" << endm;

                p.SetXYZ(0, 0, 0);
                long long duration = (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e-6; // milliseconds
                if (mGenHistograms) {
                    this->mHist["FitStatus"]->Fill("Fail", 1);
                    this->mHist["FailedFitDuration"]->Fill(duration);
                }
                return p;
            } // neither track rep converged

            p = cardinalRep->getMom(fitTrack.getFittedState(1, cardinalRep));
            mQ = cardinalRep->getCharge(fitTrack.getFittedState(1, cardinalRep));
            mP = p;

            LOG_INFO << "track fit p = " << TString::Format( "(%f, %f, %f), q=%f", p.X(), p.Y(), p.Z(), mQ ).Data() << endm;

        } catch (genfit::Exception &e) {
            LOG_WARN << "Exception on track fit: " << e.what() << endm;
            p.SetXYZ(0, 0, 0);

            long long duration = (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e-6; // milliseconds
            if (mGenHistograms) {
                this->mHist["FitStatus"]->Fill("Exception", 1);
                this->mHist["FailedFitDuration"]->Fill(duration);
            }

            return p;
        } // try/catch 

        long long duration = (FwdTrackerUtils::nowNanoSecond() - itStart) * 1e-6; // milliseconds
        if (mGenHistograms) {
            this->mHist["FitStatus"]->Fill("Pass", 1);
            this->mHist["delta_fit_seed_pT"]->Fill(p.Pt() - seedMom.Pt());
            this->mHist["delta_fit_seed_eta"]->Fill(p.Eta() - seedMom.Eta());
            this->mHist["delta_fit_seed_phi"]->Fill(p.Phi() - seedMom.Phi());
            this->mHist["FitDuration"]->Fill(duration);
        }
        return p;
    }

    int getCharge() {
        return (int)mQ;
    }

    

    // Store the planes for FTT and FST
    vector<genfit::SharedPlanePtr> mFTTPlanes;
    vector<genfit::SharedPlanePtr> mFSTPlanes;
    vector<genfit::SharedPlanePtr> mFSTSensorPlanes;
    vector<genfit::SharedPlanePtr> mFSTPlanesInner;
    vector<genfit::SharedPlanePtr> mFSTPlanesOuter;

    void SetIncludeVertex( bool vert ) { mIncludeVertexInFit = vert; }

  protected:
    std::unique_ptr<genfit::AbsBField> mBField;

    FwdTrackerConfig mConfig; // main config object

    // optional histograms, off by default
    std::map<std::string, TH1 *> mHist;
    std::map<std::string, TH2 *> mHist2D;
    bool mGenHistograms = false;
    bool mMisaligned = false;

    int mModuleMap[3][12] = {{1,6,0,11,5,10,4,9,3,8,2,7},
                             {6,0,11,5,10,4,9,3,8,2,7,1},
                             {1,6,0,11,5,10,4,9,3,8,2,7}};

    int mSensorMap[3] = {2,0,1};

    TMatrixD mInverseM[108];

    // Main GenFit fitter instance
    std::unique_ptr<genfit::AbsKalmanFitter> mFitter = nullptr;
    std::unique_ptr<StFwdGbl> mGblFitter = nullptr;

    // PDG codes for the default plc type for fits
    const int mPdgPiPlus = 211;
    const int mPdgPiMinus = -211;
    const int mPdgPositron = 11;
    const int mPdgElectron = -11;


    // det z locations loaded from geom or config
    vector<double> mFSTZLocations, mFTTZLocations;

    // parameter ALIASED from mConfig wrt PV vertex
    double mVertexSigmaXY = 1;
    double mVertexSigmaZ = 30;
    vector<double> mVertexPos;
    bool mIncludeVertexInFit = false;
    bool mSmearMcVertex = false;
    bool mRefitGBL = false;

    // GenFit state
    genfit::FitStatus mFitStatus;
    genfit::AbsTrackRep *mTrackRep;
    genfit::Track *mFitTrack;

    // Fit results
    TVector3 mP;
    double mQ;
};

#endif
