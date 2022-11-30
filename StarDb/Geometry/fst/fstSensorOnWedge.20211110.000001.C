#include "TMath.h"

TDataSet *CreateTable() {
    if (!TClass::GetClass("St_Survey")) return 0;
Survey_st row;
St_Survey *tableSet = new St_Survey("fstSensorOnWedge",108);
//
  for(int i = 0; i < 108; i++){
    if(i != 36){
      memset(&row,0,tableSet->GetRowSize());
          row.Id   = i+1000;
          row.r00  = 1.0;
          row.r01  = 0.0;
          row.r02  = 0.0;
          row.r10  = 0.0;
          row.r11  = 1.0;
          row.r12  = 0.0;
          row.r20  = 0.0;
          row.r21  = 0.0;
          row.r22  = 1.0;
          row.t0   = 0.0;
          row.t1   = 0.0;
          row.t2   = 0.0;
          memcpy(&row.comment,"Identity Matrix\x00",1);
      tableSet->AddAt(&row);
    }
    else{
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
      memset(&row,0,tableSet->GetRowSize());
          row.Id   = i+1000;
          row.r00  = 1.0;//cost;
          row.r01  = 0.0;//-sint;
          row.r02  = 0.0;
          row.r10  = 0.0;//sint;
          row.r11  = 1.0;//cost;
          row.r12  = 0.0;
          row.r20  = 0.0;
          row.r21  = 0.0;
          row.r22  = 1.0;
          row.t0   = 0.0;//dx;
          row.t1   = 0.0;//dy;
          row.t2   = 0.0;
          memcpy(&row.comment,"Identity Matrix\x00",1);
      tableSet->AddAt(&row);
  
    }
  }

  return (TDataSet *)tableSet; 
}
