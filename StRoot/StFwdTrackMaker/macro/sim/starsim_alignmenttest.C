// macro to instantiate the Geant3 from within
// STAR  C++  framework and get the starsim prompt
// To use it do
// root4star starsim.C
// This macro is capable of loading the misalignment for the FST geometry using dev2022m geomtag
// To ensure the misalignment tables load, the SDT needs to be set later than timestamp on the tables
// 	located under ./StarDb/Geometry/fst/ 
class St_geant_Maker;
St_geant_Maker *geant_maker = 0;

class StarGenEvent;
StarGenEvent   *event       = 0;

class StarPrimaryMaker;
StarPrimaryMaker *_primary  = 0;

class StarKinematics;
StarKinematics *kinematics  = 0;

int     _npart = 10;    // floor number of tracks per event
TString _part  = "mu+"; // particle to simulate
float   _ptmn  = 0.200; // min pT to simulate [GeV]
float   _ptmx  = 5.000; // max pT to simulate [GeV]
float   _etamn = 2.5;   // min eta to simulate
float   _etamx = 4.0;   // max eta to simulate
float   _phimn = 1.01229;   // min phi to simulate
float   _phimx = 1.60570;   // max phi to simulate

TString _geometry = "dev2022m";
TString DBV; 
TString SDT = "sdt20220421";

//______________________________________________________________________________________

void geometry( TString tag, Bool_t agml=true )
{
  TString cmd = "DETP GEOM "; cmd += tag;
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> LoadGeometry(cmd);
  //  if ( agml ) command("gexec $STAR_LIB/libxgeometry.so");
}

//______________________________________________________________________________________

void command( TString cmd )
{
  if ( !geant_maker ) geant_maker = (St_geant_Maker *)chain->GetMaker("geant");
  geant_maker -> Do( cmd );
}

//______________________________________________________________________________________

void trig( int n=1 )
{


  for ( int i=0; i<n; i++ ) {

    // Clear the chain from the previous event
    chain->Clear();

    _primary->SetVertex( 0.0, 0.0, 0.0 );
    _primary->SetSigma(0.0, 0.0, 0.0);

    kinematics->Kine( _npart, _part, _ptmn, _ptmx, _etamn, _etamx, _phimn, _phimx );

    // Generate the event
    chain->Make();

  }
}
//______________________________________________________________________________________

void Kinematics()
{
  gSystem->Load( "libKinematics.so");
  kinematics = new StarKinematics();
    
  _primary->AddGenerator(kinematics);
}
//______________________________________________________________________________________

void starsim ( int rngSeed=0,      
               int nevents=-1,    
               const char* outfile = "sim" )
{ 
 
  //
  //________________________________________________________
  //
  // Setup the big full chain
  //
  //________________________________________________________
   gSystem->Load( "libStarRoot.so" );
   gROOT->SetMacroPath(".:/star-sw/StRoot/macros/:./StRoot/macros:./StRoot/macros/graphics:./StRoot/macros/analysis:./StRoot/macros/test:./StRoot/macros/examples:./StRoot/macros/html:./StRoot/macros/qa:./StRoot/macros/calib:./StRoot/macros/mudst:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/graphics:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/analysis:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/test:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/examples:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/html:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/qa:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/calib:/afs/rhic.bnl.gov/star/packages/DEV/StRoot/macros/mudst:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/macros:/afs/rhic.bnl.gov/star/ROOT/36/5.34.38/.sl73_x8664_gcc485/rootdeb/tutorials");
  //

  gROOT->ProcessLine(".L bfc.C");
  {
    TString simple = _geometry; simple += " ";
    simple += SDT; simple += " ";
    simple += DBV; simple += " ";
    simple += " geant gstar FieldOff usexgeom misalign agml ";
    bfc(0, simple );
  }
  //
  //________________________________________________________
  //
  // Load in supporting libraries
  //________________________________________________________
  //
  gSystem->Load( "libVMC.so");

  gSystem->Load( "StarGeneratorUtil.so"  );
  gSystem->Load( "StarGeneratorEvent.so" );
  gSystem->Load( "StarGeneratorBase.so"  );
  gSystem->Load( "StarGeneratorDecay.so" ); 
  gSystem->Load( "libMathMore.so"        );
  // Force xgeometry.so from local repository to load   
  gSystem->Load( "xgeometry.so"          );


  //________________________________________________________
  //
  // Setup RNG seed and map all ROOT TRandom here
  //________________________________________________________
  // 
 
  StarRandom::seed( rngSeed ); // but will reset based on run number and event number
  StarRandom::capture();
  
  //
  // Create the primary event generator and insert it
  // before the geant maker
  //
  //  StarPrimaryMaker *
  _primary = new StarPrimaryMaker();
  {
    _primary -> SetFileName( "kinematics.starsim.root");
    chain -> AddBefore( "geant", _primary );
  }

  Kinematics();

  //
  // Initialize primary event generator and all sub makers
  //
  _primary -> Init();

  //
  // Setup geometry and set starsim to use agusread for input
  //
  command("gkine -4 0");

  TString fzdname = Form("gfile o %s.fzd",outfile);
  TString rooname = Form("%s.genevents.root",outfile);
  command( fzdname );
  _primary -> SetFileName( rooname );

  //
  // Trigger on nevents
  //
  trig( nevents );

  command("call agexit");  // Make sure that STARSIM exits properly

}

//______________________________________________________________________________________

void starsim_alignmenttest ( const char* outfile,
                             const int nevents,
                             const int np,
                             const char* part,
                             const float ptmn,
                             const float ptmx,
	                     const float etamn,
                             const float etamx,
	                     const float phimn,
                             const float phimx,
                             const int seed,
                             const char* geom,
                             const char* dbv = 0)
{
  _npart = np;
  _part  = part;
  _ptmn  = ptmn;
  _ptmx  = ptmx;
  _etamn = etamn;
  _etamx = etamx;
  _phimn = phimn;
  _phimx = phimx;
  _geometry = geom;
  if ( dbv ) DBV = dbv;
  starsim( seed, nevents, outfile );
}
//______________________________________________________________________________________

