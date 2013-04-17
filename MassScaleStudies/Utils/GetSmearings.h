#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>


float GetSmearings(const float& scEta, const float& r9, const int& year, const int& isEB, const float& syst = 0.)
{ 
  float smear = 1.;
  float smearErr = 0.;
  
  //2012 
  if(year == 2012)
  {
    //EB LR9
    if( isEB == 1 && fabs(scEta) < 1. && r9 < 0.94) smear    = 0.0107;
    if( isEB == 1 && fabs(scEta) < 1. && r9 < 0.94) smearErr = 0.0025;
    if( isEB == 1 && fabs(scEta) > 1. && r9 < 0.94) smear    = 0.0194;
    if( isEB == 1 && fabs(scEta) > 1. && r9 < 0.94) smearErr = 0.0060;
    //EB HR9
    if( isEB == 1 && fabs(scEta) < 1. && r9 > 0.94) smear    = 0.0111;
    if( isEB == 1 && fabs(scEta) < 1. && r9 > 0.94) smearErr = 0.0023;
    if( isEB == 1 && fabs(scEta) > 1. && r9 > 0.94) smear    = 0.0155;
    if( isEB == 1 && fabs(scEta) > 1. && r9 > 0.94) smearErr = 0.0072;
    //EE LR9
    if( isEB == 0 && fabs(scEta) < 2. && r9 < 0.94) smear    = 0.0276;
    if( isEB == 0 && fabs(scEta) < 2. && r9 < 0.94) smearErr = 0.0033;
    if( isEB == 0 && fabs(scEta) > 2. && r9 < 0.94) smear    = 0.0371;
    if( isEB == 0 && fabs(scEta) > 2. && r9 < 0.94) smearErr = 0.0054;
    //EE HR9
    if( isEB == 0 && fabs(scEta) < 2. && r9 > 0.94) smear    = 0.0295;
    if( isEB == 0 && fabs(scEta) < 2. && r9 > 0.94) smearErr = 0.0093;
    if( isEB == 0 && fabs(scEta) > 2. && r9 > 0.94) smear    = 0.0370;
    if( isEB == 0 && fabs(scEta) > 2. && r9 > 0.94) smearErr = 0.0036;
  }

  //2011  // from AN2011-426
  if(year == 2011)
  {
    if( isEB == 1 && fabs(scEta) < 1. && r9 > 0.94) smear = 0.0089;
    if( isEB == 1 && fabs(scEta) < 1. && r9 < 0.94) smear = 0.0109;

    if( isEB == 1 && fabs(scEta) > 1. && r9 > 0.94) smear = 0.0156;
    if( isEB == 1 && fabs(scEta) > 1. && r9 < 0.94) smear = 0.0203;

    if( isEB == 0 && fabs(scEta) < 2. && r9 > 0.94) smear = 0.0303;
    if( isEB == 0 && fabs(scEta) < 2. && r9 < 0.94) smear = 0.0326;

    if( isEB == 0 && fabs(scEta) > 2. && r9 > 0.94) smear = 0.0318;
    if( isEB == 0 && fabs(scEta) > 2. && r9 < 0.94) smear = 0.0331;
  }
  
  return smear + syst*smearErr;
};
