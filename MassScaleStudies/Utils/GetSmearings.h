#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>


float GetSmearings(const float& scEta, const float& r9, const int& year, const int& isEB){ 

  float smear = 1.;

  //2012 
  if(year == 2012){
    //EB LR9
    if( isEB == 1 && fabs(scEta) < 1. && r9 < 0.94) smear = 0.0075;
    if( isEB == 1 && fabs(scEta) > 1. && r9 < 0.94) smear = 0.0075;
    //EB HR9
    if( isEB == 1 && fabs(scEta) < 1. && r9 > 0.94) smear = 0.0075;
    if( isEB == 1 && fabs(scEta) > 1. && r9 > 0.94) smear = 0.0075;
    //EE LR9
    if( isEB == 0 && fabs(scEta) < 2. && r9 < 0.94) smear = 0.0075;
    if( isEB == 0 && fabs(scEta) > 2. && r9 < 0.94) smear = 0.0075;
    //EE HR9
    if( isEB == 0 && fabs(scEta) < 2. && r9 > 0.94) smear = 0.0075;
    if( isEB == 0 && fabs(scEta) > 2. && r9 > 0.94) smear = 0.0075;
  }

  //2011  // from AN2011-426
  if(year == 2011){

    if( isEB == 1 && fabs(scEta) < 1. && r9 > 0.94) smear = 0.0089;
    if( isEB == 1 && fabs(scEta) < 1. && r9 < 0.94) smear = 0.0109;

    if( isEB == 1 && fabs(scEta) > 1. && r9 > 0.94) smear = 0.0156;
    if( isEB == 1 && fabs(scEta) > 1. && r9 < 0.94) smear = 0.0203;

    if( isEB == 0 && fabs(scEta) < 2. && r9 > 0.94) smear = 0.0303;
    if( isEB == 0 && fabs(scEta) < 2. && r9 < 0.94) smear = 0.0326;

    if( isEB == 0 && fabs(scEta) > 2. && r9 > 0.94) smear = 0.0318;
    if( isEB == 0 && fabs(scEta) > 2. && r9 < 0.94) smear = 0.0331;
  }
  return smear;
};

