TRandom3 *randa = new TRandom3(615);

double aValues(){

   double numbers = randa->Uniform();
   return -240 + 2*240*numbers;  
} 


double bValues(){

    double numbers = randa->Uniform();
    return -240 + 2*240*numbers;  
}


double phiValues(){

   return randa->Uniform()*2*TMath::Pi();  
}


double thetaValue(){ 
   double x = 0.0;
   double y = 1000.0;
   while (y > TMath::Power(x,2)) {
      x = randa->Uniform(0,1);
      y = randa->Uniform(0,1);
        }
    double th = TMath::ACos(x); 
    return th;
}

// Length distributions

double l1(double a, double b, double Theta, double Phi){

   double L = 123.0; //cm
   double R = 102.0; //cm
   double Rs = 800.; // cm 
   
   // start point of the line
   double x0 = Rs*TMath::Sin(Theta)*TMath::Cos(Phi) + a*TMath::Cos(Phi)*TMath::Cos(Theta) - b*TMath::Sin(Phi);
   double y0 = Rs*TMath::Sin(Theta)*TMath::Sin(Phi)+a*TMath::Sin(Phi)*TMath::Cos(Theta)+b*TMath::Cos(Phi);
   double z0 = Rs*TMath::Cos(Theta)-a*TMath::Sin(Theta);

   // interceptions upper plane  (Point P1) 
   double xt = x0+(L-z0)*TMath::Tan(Theta)*TMath::Cos(Phi); 
   double yt = y0+(L-z0)*TMath::Tan(Theta)*TMath::Sin(Phi);
   double zt = L;
    
   // interceptions lower plane (Point P2)
   double xb = x0-(z0)*TMath::Tan(Theta)*TMath::Cos(Phi);
   double yb = y0-(z0)*TMath::Tan(Theta)*TMath::Sin(Phi);
   double zb = 0.0;

   // interceptions with the face cylinder  
   double A = TMath::Power(TMath::Sin(Theta),2);
   double B = -2*TMath::Sin(Theta)*(x0*TMath::Cos(Phi)+y0*TMath::Sin(Phi));
   double C = x0*x0+y0*y0-R*R ;
   // only interceptions
   
   double Discriminante = B*B - 4*A*C;
    
   if (Discriminante > 0){
       double T1 = 1/(2*A)*(-B+TMath::Sqrt(B*B-4*A*C));
           double T2 = 1/(2*A)*(-B-TMath::Sqrt(B*B-4*A*C));
       // two interceptions 
           double t[2] = {T1, T2};
       // first solution (Point P3)
           double  x1 = x0-t[0]*TMath::Sin(Theta)*TMath::Cos(Phi);
           double  y1 = y0-t[0]*TMath::Sin(Theta)*TMath::Sin(Phi);
           double  z1 = z0-t[0]*TMath::Cos(Theta);
        // second solution (Point P4)
           double x2 = x0-t[1]*TMath::Sin(Theta)*TMath::Cos(Phi);
           double y2 = y0-t[1]*TMath::Sin(Theta)*TMath::Sin(Phi);
           double z2 = z0-t[1]*TMath::Cos(Theta);
          // first case
            if ( (xt*xt+yt*yt < R*R) && (xb*xb+yb*yb) < R*R){
                 double Puntos[6] = {xt, yt , L, xb , yb, 0.0};
                 return TMath::Sqrt(TMath::Power((Puntos[3]-Puntos[0]),2)+TMath::Power((Puntos[4]-Puntos[1]),2)+TMath::Power((Puntos[5]-Puntos[2]),2) );
            }
            // second case 
            if ((xt*xt+yt*yt) < R*R){
               if ((z1 > 0.0)  &&  (z1 < L)){
                   double Puntos1[6] = {xt, yt, L, x1, y1, z1};
                   return TMath::Sqrt(TMath::Power((Puntos1[3]-Puntos1[0]),2)+TMath::Power((Puntos1[4]-Puntos1[1]),2)+TMath::Power((Puntos1[5]-Puntos1[2]),2) );
               }
                else if ((z2 > 0.0) && (z2 < L)){
                   double Puntos2[6] = {xt, yt, L, x2, y2, z2};
                   return  TMath::Sqrt( TMath::Power((Puntos2[3]-Puntos2[0]),2)+TMath::Power((Puntos2[4]-Puntos2[1]),2)+TMath::Power((Puntos2[5]-Puntos2[2]),2) );
                }
            }
             // thirth case 
            if (((z1 > 0.0)  &&  (z1 < L)) && ((z2 > 0.0)  &&  (z2 < L))){
                double Puntos3[6] = {x1, y1, z1, x2, y2, z2};
                return TMath::Sqrt(TMath::Power((Puntos3[3]-Puntos3[0]),2)+TMath::Power((Puntos3[4]-Puntos3[1]),2)+TMath::Power((Puntos3[5]-Puntos3[2]),2) ); 
               }
               //  fourt case  
            if ((xb*xb+yb*yb) < R*R){
                 if ((z1 > 0.0)  &&  (z1 < L)){
                    double Puntos4[6] = {xb, yb, 0.0 , x1, y1, z1};
                    return TMath::Sqrt(TMath::Power((Puntos4[3]-Puntos4[0]),2)+TMath::Power((Puntos4[4]-Puntos4[1]),2)+TMath::Power((Puntos4[5]-Puntos4[2]),2) );
               }
                 else if ((z2 > 0.0) && (z2 < L)){
                   double Puntos5[6] = {xb, yb, 0.0, x2, y2, z2};
                   return TMath::Sqrt( TMath::Power((Puntos5[3]-Puntos5[0]),2)+TMath::Power((Puntos5[4]-Puntos5[1]),2)+TMath::Power((Puntos5[5]-Puntos5[2]),2) );
                  }
            }
            
      }
    
      return 0.0;
    
}

double len(){

 double a, b, theta, phi;
 double lenf = -20.0;
  
      while(lenf <= 0){
   //branch values of a 
     a = aValues();
   //branch values of b 
     b = bValues();
   //branch values of theta 
     theta = thetaValue();
   //branch values of phi
     phi = phiValues();
   //branch lengths 
     lenf = l1(a, b, theta, phi);
      }
      
return lenf;

}

void lengthvalues(){

TFile *lengthDvalues = new TFile("lengthDvalues.root", "RECREATE");
TTree *data = new TTree("data", "lengthDvalues");  

double l;
int NE = 100000;

data->Branch("NE",&NE, "NE/I");
data->Branch("l", &l, "l/D");

for (int i = 0; i<NE;i++){

l = len();

data->Fill();

}	

data->Write();
lengthDvalues->Close();
}


