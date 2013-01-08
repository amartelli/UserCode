#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>


float GetShervingCorrections(const float& scEta, const float& r9, const int& runId){ 

  float corr = 1.;

  //2012 EB
  if( fabs(scEta) < 1. && r9 < 0.94) {
    if(runId >    190645 && runId <  190781) corr = 1.0020;
    if(runId >    190782 && runId <  191042) corr = 1.0079;
    if(runId >    191043 && runId <  193555) corr = 0.9991;
    if(runId >    193556 && runId <  194150) corr = 0.9976;
    if(runId >    194151 && runId <  194532) corr = 0.9981;
    if(runId >    194533 && runId <  195113) corr = 0.9983;
    if(runId >    195114 && runId <  195915) corr = 0.9984;
    if(runId >    195916 && runId <  198115) corr = 0.9975;
    if(runId >    198116 && runId <  199803) corr = 1.0011;
    if(runId >    199804 && runId <  200048) corr = 1.0023;
    if(runId >    200049 && runId <  200151) corr = 1.0036;
    if(runId >    200152 && runId <  200490) corr = 1.0013;
    if(runId >    200491 && runId <  200531) corr = 1.0035;
    if(runId >    200532 && runId <  201656) corr = 1.0017;
    if(runId >    201657 && runId <  202305) corr = 1.0025;
    if(runId >    202305 && runId <  203002) corr = 1.0038;
  }
  if( fabs(scEta) < 1. && r9 > 0.94) {
    if(runId > 190645 && runId <  190781) corr = 0.9964;
    if(runId > 190782 && runId <  191042) corr = 1.0024;
    if(runId > 191043 && runId <  193555) corr = 0.9935;
    if(runId > 193556 && runId <  194150) corr = 0.9920;
    if(runId > 194151 && runId <  194532) corr = 0.9925;
    if(runId > 194533 && runId <  195113) corr = 0.9927;
    if(runId > 195114 && runId <  195915) corr = 0.9929;
    if(runId > 195916 && runId <  198115) corr = 0.9919;
    if(runId > 198116 && runId <  199803) corr = 0.9955;
    if(runId > 199804 && runId <  200048) corr = 0.9967;
    if(runId > 200049 && runId <  200151) corr = 0.9980;
    if(runId > 200152 && runId <  200490) corr = 0.9958;
    if(runId > 200491 && runId <  200531) corr = 0.9979;
    if(runId > 200532 && runId <  201656) corr = 0.9961;
    if(runId > 201657 && runId <  202305) corr = 0.9969;
    if(runId > 202305 && runId <  203002) corr = 0.9982;
  }

  //2012 EB 4
  if(fabs(scEta) > 1. && fabs(scEta) < 1.4442 && r9 < 0.94){
    if(runId >190645 && runId <  190781) corr = 1.0028;
    if(runId >190782 && runId <  191042) corr = 1.0058;
    if(runId >191043 && runId <  193555) corr = 0.9997;
    if(runId >193556 && runId <  194150) corr = 0.9951;
    if(runId >194151 && runId <  194532) corr = 0.9963;
    if(runId >194533 && runId <  195113) corr = 0.9981;
    if(runId >195114 && runId <  195915) corr = 0.9974;
    if(runId >195916 && runId <  198115) corr = 0.9964;
    if(runId >198116 && runId <  199803) corr = 0.9995;
    if(runId >199804 && runId <  200048) corr = 1.0006;
    if(runId >200049 && runId <  200151) corr = 1.0012;
    if(runId >200152 && runId <  200490) corr = 1.0004;
    if(runId >200491 && runId <  200531) corr = 1.0012;
    if(runId >200532 && runId <  201656) corr = 0.9996;
    if(runId >201657 && runId <  202305) corr = 1.0002;
    if(runId >202305 && runId <  203002) corr = 1.0008;
  }
  if(fabs(scEta) > 1. && fabs(scEta) < 1.4442 && r9 > 0.94){
    if(runId >190645 && runId <  190781) corr = 0.9893;
    if(runId >190782 && runId <  191042) corr = 0.9923;
    if(runId >191043 && runId <  193555) corr = 0.9861;
    if(runId >193556 && runId <  194150) corr = 0.9814;
    if(runId >194151 && runId <  194532) corr = 0.9826;
    if(runId >194533 && runId <  195113) corr = 0.9844;
    if(runId >195114 && runId <  195915) corr = 0.9838;
    if(runId >195916 && runId <  198115) corr = 0.9827;
    if(runId >198116 && runId <  199803) corr = 0.9859;
    if(runId >199804 && runId <  200048) corr = 0.9870;
    if(runId >200049 && runId <  200151) corr = 0.9877;
    if(runId >200152 && runId <  200490) corr = 0.9868;
    if(runId >200491 && runId <  200531) corr = 0.9876;
    if(runId >200532 && runId <  201656) corr = 0.9860;
    if(runId >201657 && runId <  202305) corr = 0.9866;
    if(runId >202305 && runId <  203002) corr = 0.9872;
  }

  //2012 EE
  if(fabs(scEta) > 1.566 && fabs(scEta) < 2. && r9 < 0.94){
    if(runId >190645 && runId <  190781) corr = 0.9937; 
    if(runId >190782 && runId <  191042) corr = 0.9977 ;
    if(runId >191043 && runId <  193555) corr = 0.9960 ;
    if(runId >193556 && runId <  194150) corr = 0.9962 ;
    if(runId >194151 && runId <  194532) corr = 0.9980 ;
    if(runId >194533 && runId <  195113) corr = 0.9999 ;
    if(runId >195114 && runId <  195915) corr = 1.0007 ;
    if(runId >195916 && runId <  198115) corr = 1.0017 ;
    if(runId >198116 && runId <  199803) corr = 0.9959 ;
    if(runId >199804 && runId <  200048) corr = 0.9959 ;
    if(runId >200049 && runId <  200151) corr = 0.9976 ;
    if(runId >200152 && runId <  200490) corr = 0.9988 ;
    if(runId >200491 && runId <  200531) corr = 0.9981 ;
    if(runId >200532 && runId <  201656) corr = 0.9970 ;
    if(runId >201657 && runId <  202305) corr = 0.9980 ;
    if(runId >202305 && runId <  203002) corr = 1.0000 ;
  }
  if(fabs(scEta) > 1.566 && fabs(scEta) < 2. && r9 > 0.94){
    if(runId >190645 && runId <  190781) corr = 0.9871; 
    if(runId >190782 && runId <  191042) corr = 0.9911; 
    if(runId >191043 && runId <  193555) corr = 0.9894; 
    if(runId >193556 && runId <  194150) corr = 0.9896; 
    if(runId >194151 && runId <  194532) corr = 0.9914; 
    if(runId >194533 && runId <  195113) corr = 0.9934; 
    if(runId >195114 && runId <  195915) corr = 0.9942; 
    if(runId >195916 && runId <  198115) corr = 0.9952; 
    if(runId >198116 && runId <  199803) corr = 0.9893; 
    if(runId >199804 && runId <  200048) corr = 0.9893; 
    if(runId >200049 && runId <  200151) corr = 0.9910; 
    if(runId >200152 && runId <  200490) corr = 0.9922; 
    if(runId >200491 && runId <  200531) corr = 0.9915; 
    if(runId >200532 && runId <  201656) corr = 0.9904; 
    if(runId >201657 && runId <  202305) corr = 0.9914; 
    if(runId >202305 && runId <  203002) corr = 0.9934; 
  }

  //2012 EE far
  if(fabs(scEta) > 2. && fabs(scEta) < 2.5 && r9 < 0.94){
    if(runId > 190645 && runId <  190781) corr = 0.9958; 
    if(runId > 190782 && runId <  191042) corr = 1.0005; 
    if(runId > 191043 && runId <  193555) corr = 0.9982; 
    if(runId > 193556 && runId <  194150) corr = 0.9990; 
    if(runId > 194151 && runId <  194532) corr = 0.9993; 
    if(runId > 194533 && runId <  195113) corr = 0.9996; 
    if(runId > 195114 && runId <  195915) corr = 0.9997; 
    if(runId > 195916 && runId <  198115) corr = 0.9987; 
    if(runId > 198116 && runId <  199803) corr = 1.0041; 
    if(runId > 199804 && runId <  200048) corr = 1.0055; 
    if(runId > 200049 && runId <  200151) corr = 1.0097; 
    if(runId > 200152 && runId <  200490) corr = 1.0065; 
    if(runId > 200491 && runId <  200531) corr = 1.0096; 
    if(runId > 200532 && runId <  201656) corr = 1.0063; 
    if(runId > 201657 && runId <  202305) corr = 1.0116; 
    if(runId > 202305 && runId <  203002) corr = 1.0135; 
  }
  if(fabs(scEta) > 2. && fabs(scEta) < 2.5 && r9 > 0.94){
    if(runId > 190645 && runId <  190781) corr = 0.9839;
    if(runId > 190782 && runId <  191042) corr = 0.9886;
    if(runId > 191043 && runId <  193555) corr = 0.9864;
    if(runId > 193556 && runId <  194150) corr = 0.9872;
    if(runId > 194151 && runId <  194532) corr = 0.9874;
    if(runId > 194533 && runId <  195113) corr = 0.9878;
    if(runId > 195114 && runId <  195915) corr = 0.9878;
    if(runId > 195916 && runId <  198115) corr = 0.9869;
    if(runId > 198116 && runId <  199803) corr = 0.9923;
    if(runId > 199804 && runId <  200048) corr = 0.9937;
    if(runId > 200049 && runId <  200151) corr = 0.9980;
    if(runId > 200152 && runId <  200490) corr = 0.9948;
    if(runId > 200491 && runId <  200531) corr = 0.9979;
    if(runId > 200532 && runId <  201656) corr = 0.9945;
    if(runId > 201657 && runId <  202305) corr = 0.9999;
    if(runId > 202305 && runId <  203002) corr = 1.0018;
  }


  //2011 EB
  if(fabs(scEta) < 1. && r9 < 0.94 ){
    if(runId > 160431 && runId < 167913) corr = 1. - ( -0.0029);
    if(runId > 170000 && runId < 172619) corr = 1. - ( -0.0025);
    if(runId > 172620 && runId < 173692) corr = 1. - ( -0.0025);
    if(runId > 175860 && runId < 177139) corr = 1. - ( -0.0023);
    if(runId > 177140 && runId < 178421) corr = 1. - ( -0.0025);
    if(runId > 178424 && runId < 180252) corr = 1. - ( -0.0031);
  }
  if(fabs(scEta) < 1. && r9 > 0.94 ){
    if(runId > 160431 && runId < 167913) corr = 1.00 + ( -0.0095);
    if(runId > 170000 && runId < 172619) corr = 1.00 + ( -0.0091);
    if(runId > 172620 && runId < 173692) corr = 1.00 + ( -0.0091);
    if(runId > 175860 && runId < 177139) corr = 1.00 + ( -0.0089);
    if(runId > 177140 && runId < 178421) corr = 1.00 + ( -0.009); 
    if(runId > 178424 && runId < 180252) corr = 1.00 + ( -0.0097);
  }

  //2011 EB 4
  if(fabs(scEta) > 1. && fabs(scEta) < 1.4442 && r9 < 0.94 ){
    if(runId > 160431 && runId < 167913) corr = 1. - ( 0.0094); 
    if(runId > 170000 && runId < 172619) corr = 1. - ( 0.0112); 
    if(runId > 172620 && runId < 173692) corr = 1. - ( 0.0096); 
    if(runId > 175860 && runId < 177139) corr = 1. - ( 0.0109); 
    if(runId > 177140 && runId < 178421) corr = 1. - ( 0.0105); 
    if(runId > 178424 && runId < 180252) corr = 1. - ( 0.0095); 
  }
  if(fabs(scEta) > 1. && fabs(scEta) < 1.4442 && r9 > 0.94 ){
    if(runId > 160431 && runId < 167913) corr = 1.00 + ( -0.0024);
    if(runId > 170000 && runId < 172619) corr = 1.00 + ( -0.0006);
    if(runId > 172620 && runId < 173692) corr = 1.00 + ( -0.0023);
    if(runId > 175860 && runId < 177139) corr = 1.00 + ( -0.001); 
    if(runId > 177140 && runId < 178421) corr = 1.00 + ( -0.0013);
    if(runId > 178424 && runId < 180252) corr = 1.00 + ( -0.0024);
  }

  //2011 EE
  if(fabs(scEta) > 1.566 && fabs(scEta) < 2. && r9 < 0.94 ){
    if(runId > 160431 && runId < 167913) corr = 1. - ( 0.0044); 
    if(runId > 170000 && runId < 172619) corr = 1. - ( 0.0012); 
    if(runId > 172620 && runId < 173692) corr = 1. - ( 0.0013); 
    if(runId > 175860 && runId < 177139) corr = 1. - ( -0.0027);
    if(runId > 177140 && runId < 178421) corr = 1. - ( -0.0028);
    if(runId > 178424 && runId < 180252) corr = 1. - ( -0.006); 
  }
 if(fabs(scEta) > 1.566 && fabs(scEta) < 2. && r9 > 0.94 ){
   if(runId > 160431 && runId < 167913) corr = 1.00 + ( -0.0006);
   if(runId > 170000 && runId < 172619) corr = 1.00 + ( -0.0038);
   if(runId > 172620 && runId < 173692) corr = 1.00 + ( -0.0037);
   if(runId > 175860 && runId < 177139) corr = 1.00 + ( -0.0078);
   if(runId > 177140 && runId < 178421) corr = 1.00 + ( -0.0079);
   if(runId > 178424 && runId < 180252) corr = 1.00 + ( -0.0111);
 }

  //2011 EE far
 if(fabs(scEta) > 2. && fabs(scEta) < 2.5 && r9 < 0.94){
   if(runId > 160431 && runId < 167913) corr = 1. - ( 0.0079); 
   if(runId > 170000 && runId < 172619) corr = 1. - ( 0.0072); 
   if(runId > 172620 && runId < 173692) corr = 1. - ( 0.0057); 
   if(runId > 175860 && runId < 177139) corr = 1. - ( 0.0077); 
   if(runId > 177140 && runId < 178421) corr = 1. - ( 0.0085); 
   if(runId > 178424 && runId < 180252) corr = 1. - ( 0.0086); 
 }
 if(fabs(scEta) > 2. && fabs(scEta) < 2.5 && r9 > 0.94){
   if(runId > 160431 && runId < 167913) corr = 1.00 + ( -0.0032); 
   if(runId > 170000 && runId < 172619) corr = 1.00 + ( -0.0038); 
   if(runId > 172620 && runId < 173692) corr = 1.00 + ( -0.0053); 
   if(runId > 175860 && runId < 177139) corr = 1.00 + ( -0.0033); 
   if(runId > 177140 && runId < 178421) corr = 1.00 + ( -0.0025); 
   if(runId > 178424 && runId < 180252) corr = 1.00 + ( -0.0024); 
 }

  return corr;

};

