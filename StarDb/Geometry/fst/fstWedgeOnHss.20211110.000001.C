TDataSet *CreateTable() {
    if (!TClass::GetClass("St_Survey")) return 0;
Survey_st row;
St_Survey *tableSet = new St_Survey("fstWedgeOnHss",36);
//
for(int i = 0; i < 36; i++){
  //if(i != 12){
    memset(&row,0,tableSet->GetRowSize());
        row.Id   = i+1;
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
  //}
  //else{
  //  double cos02 = TMath::Cos(0.02);    
  //  double sin02 = TMath::Sin(0.02);  
  //  //double cos75 = TMath::Cos(1.309);    
  //  //double sin75 = TMath::Sin(1.309);  
  //  //double dx = -sin75*0.1; // dv = 500 um   
  //  //double dy =  cos75*0.1; // dv = 500 um   
  //  memset(&row,0,tableSet->GetRowSize());
  //      row.Id   = i+1;
  //      row.r00  = cos01;
  //      row.r01  = -sin02;
  //      row.r02  = 0.0;
  //      row.r10  = sin02;
  //      row.r11  = cos02;
  //      row.r12  = 0.0;
  //      row.r20  = 0.0;
  //      row.r21  = 0.0;
  //      row.r22  = 1.0;
  //      row.t0   = 0.0;
  //      row.t1   = 0.0;
  //      row.t2   = 0.0;
  //      memcpy(&row.comment,"Identity Matrix\x00",1);
  //  tableSet->AddAt(&row);
  //}
}
return (TDataSet *)tableSet; 
}
