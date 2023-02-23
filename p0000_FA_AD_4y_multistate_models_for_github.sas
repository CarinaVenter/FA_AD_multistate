*** case-cohort design   ***;
*** multistate modelling ***;
/*********************************************************************************************
Date: 11/03/2022
Written By: Kaci Pickett (adapted from code by Jaakko Nevalainen)
Program: 
Purpose: 
**********************************************************************************************/


/*********************************
Define libnames for this program
**********************************/
options nofmterr; 
libname ModData "&ProjModDat.";



data midi0;
	set ModData.Allergy_timing;
run;

/*proc contents data = midi0; run;*/


***************************************************************************;
****here is where we'll need to censor at 4 years;
data midi1;
  set midi0;
  tHH = 3.99; *censoring horizon;
  tHF = food_yrs; 
  tHA = ad_yrs;
  if tHH>0 and (tHA>0 or tHA=.) and (tHF>0 or tHF=.);																																															
run;

data t;
  set midi1;
run;

**set your knot positions;
%let cp1 = 0.5;
%let cp2 = 1;

* we need a large set of time variables *;
*H = No disease, A = Atopic Dermatitis (AD), F = Food Allergy(FA);
*HA = healthy to AD, FA = FA to AD, etc;
data timevars;
  set midi1;
  tH1 = min(&cp1,tHH);
  tH2 = min(max(0,tHH-&cp1),&cp2-&cp1);
  tH3 =  max(0,tHH-&cp2);

  if tHF ne . then do; 
    tF1 = min(&cp1,tHF);
    tF2 = min(max(0,tHF-&cp1),&cp2-&cp1);
    tF3 = max(0,tHA-&cp2);
  end;

  if tHA ne . then do; 
    tA1 = min(&cp1,tHA);
    tA2 = min(max(0,tHA-&cp1),&cp2-&cp1);
    tA3 =max(0,tHA-&cp2);
  end;

  if tHF ne . and tHA ne . and tHF<tHA then do; 
    tFA1 = max(0,min(&cp1,tHA)-tHF);
    tFA2 = max(0,min(&cp2,tHA)-max(&cp1,tHF));
    tFA3 = tHA-tHF-tFA1-tFA2;
  end;

  if tHF ne . and tHA ne . and tHF>tHA then do; 
    tAF1 = max(0,min(&cp1,tHF)-tHA);
    tAF2 = max(0,min(&cp2,tHF)-max(&cp1,tHA));
    tAF3 = tHF-tHA-tAF1-tAF2;
  end;

  if tHF ne . and tHA = . then do; 
    tFH1 = max(0,min(&cp1,tHH)-tHF);
    tFH2 = max(0,min(&cp2,tHH)-max(&cp1,tHF));
    tFH3 = tHH-tHF-tFH1-tFH2;
  end;

  if tHF = . and tHA ne . then do; 
    tAH1 = max(0,min(&cp1,tHH)-tHA);
    tAH2 = max(0,min(&cp2,tHH)-max(&cp1,tHA));
    tAH3 = tHH-tHA-tAH1-tAH2;
  end;
run;


*** *** piecewise exponential model *** ***;

proc nlmixed data = timevars;
***starting parameters for baseline hazards/transitions... d = Food allergy, a = atopic dermatitis, theta = hazard ratio(HR) for AD ->FA , tau = HR for FA ->AD;
  parms d1 = 0.0003 d2 = 0.00048 d3 = 0.0006 a1 = 0.003 a2 = 0.005 a3 = 0.0046 theta = 0.6 tau = -0.29;
  d1a = exp(theta)*d1;
  d2a = exp(theta)*d2; 
  d3a = exp(theta)*d3;

  a1d = exp(tau)*a1;
  a2d = exp(tau)*a2;
  a3d = exp(tau)*a3;

  HF = ((0 le tHF < &cp1)*d1 + (&cp1 le tHF < &cp2)*d2 + (&cp2 le tHF)*d3 );
  SF = exp(-d1*tF1-d2*tF2-d3*tF3);
  hA = ((0 le tHA < &cp1)*a1 + (&cp1 le tHA < &cp2)*a2 + (&cp2 le tHA)*a3 );
  SA = exp(-a1*tA1-a2*tA2-a3*tA3);
  SH = exp(-d1*tH1-d2*tH2-d3*tH3-a1*tH1-a2*tH2-a3*tH3); 

  *** total sampling probability ***;
  pS = 1; **0.1*SH+1*(1-SH);

  pHH = 1; * no transition from healthy state*;
  if food_yrs = . and ad_yrs = . then
    pHH = SH; 

  pHF = 1; * food allergy transition only *;
  if food_yrs ne . and ad_yrs = . then 
    pHF = HF*SF*exp(-a1*tF1-a2*tF2-a3*tF3)*exp(-a1d*tFH1-a2d*tFH2-a3d*tFH3);

  pHA = 1; * ad transition only *;
  if ad_yrs ne . and food_yrs = . then 
    pHA = hA*SA*exp(-d1*tA1-d2*tA2-d3*tA3)*exp(-d1a*tAH1-d2a*tAH2-d3a*tAH3);

  pHFA = 1; * first food allergy, then ad*;
  if ad_yrs ne . and food_yrs ne . and ad_yrs > food_yrs then do;
    HFA = ((0 le tHA < &cp1)*a1d + (&cp1 le tHA < &cp2)*a2d + (&cp2 le tHA)*a3d);
    pHFA = HF*SF*HFA*exp(-a1*tF1-a2*tF2-a3*tF3)*exp(-a1d*tFA1-a2d*tFA2-a3d*tFA3);
  end;

  pHAF = 1; * first ad, then food allergy *;
  if ad_yrs ne . and food_yrs ne . and ad_yrs < food_yrs then do;
    hAF = ((0 le tHF < &cp1)*d1a + (&cp1 le tHF < &cp2)*d2a + (&cp2 le tHF)*d3a );
    pHAF = hA*SA*hAF*exp(-d1*tA1-d2*tA2-d3*tA3)*exp(-d1a*tAF1-d2a*tAF2-d3a*tAF3);  
  end;

  ll = log(pHH*pHF*pHA*pHFA*pHAF/pS);
  model tHH~general(ll);
run;


