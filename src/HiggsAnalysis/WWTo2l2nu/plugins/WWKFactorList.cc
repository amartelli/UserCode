#include "HiggsAnalysis/WWTo2l2nu/plugins/WWKFactorList.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <stdexcept>

using namespace std;


const unsigned WWKfactorList::lineSize_ = 10000;


WWKfactorList::WWKfactorList(const char* name, unsigned nbinspt,
                               double minpt, double maxpt, double value)
  : TH1D(name, name, nbinspt, minpt, maxpt), alternativeK_(1.), alternativeNNLOK_(1.)
{

  GetXaxis()->SetTitle("#pt");
  GetYaxis()->SetTitle("Kfactor");
  if(value>0)
  {
    for(int ieta=1; ieta<=GetNbinsX(); ieta++)
    {
      SetBinContent(ieta, value);
    }
  }

}

WWKfactorList::WWKfactorList(const char* name, const char* mapfile)
  : alternativeK_(1.), alternativeNNLOK_(1.)
{

  SetTitle(mapfile);
  GetXaxis()->SetTitle("#eta");
  GetYaxis()->SetTitle("E");
  if( ! ReadMapFile(mapfile) )
  {
    string err = "WWKfactorList::WWKfactorList : cannot read file ";
    err += mapfile;
    throw invalid_argument(err);
  }
}


bool WWKfactorList::WriteMapFile(const char* mapfile)
{
  // open the file
  ofstream outf(mapfile);
  if( !outf.good() )
  {
    cout<<"WWKfactorList::Write : cannot open file "<<mapfile<<endl;
    return false;
  }

  outf<<(*this)<<endl;
  if(!outf.good() )
  {
    cerr<<"WWKfactorList::Write : corrupted file "<<mapfile<<endl;
    return false;
  }
  else
  {
    mapFile_ = mapfile;
    return true;
  }
}



bool WWKfactorList::ReadMapFile(const char* mapfile)
{
  // open the file
  ifstream inf(mapfile);
  if( !inf.good() )
  {
    return false;
  }
  // first data describes the map: histo bin, max, min
  int nbinspt=0;
  double minpt=0;
  double maxpt=0;
  std::string dummy;
  inf>>nbinspt;inf>>dummy;
  inf>>minpt;inf>>dummy;
  inf>>maxpt;inf>>dummy;
  inf>>this->alternativeK_;inf>>dummy;
  inf>>this->alternativeNNLOK_;inf>>dummy;inf>>dummy;
  SetBins(nbinspt, minpt, maxpt);

  char s[lineSize_];
  // get position in stream
  int pos=inf.tellg(); int j=0;
  // parse map data
  do
  {
    inf.seekg(pos);
    inf.getline(s,lineSize_);
    pos = inf.tellg();
    //    cout<<"LINE"<< s<<endl;
    if(string(s).empty())
    {
      continue; // remove empty lines
    }
    istringstream lin(s);
    double dataw;

    if (lin.good())
    {
      lin>>dataw;   lin>>dataw;
      //  cout<<"LINE "<<dataw <<endl;
      SetBinContent(j, dataw);
      j++;
    }

  }
  while(inf.good());

  if(inf.eof())
  {
    mapFile_ = mapfile;
    return true;
  }
  else return false;
  mapFile_ = mapfile;
  return true;
}




ostream& operator<<(ostream& outf, const WWKfactorList& rm)
{

  if(!outf.good() ) return outf;

  // first data describes the map
  outf<<rm.GetNbinsX()<<endl;
  outf<<rm.GetXaxis()->GetXmin()<<endl;
  outf<<rm.GetXaxis()->GetXmax()<<endl;

  for(int ieta=0; ieta<=rm.GetNbinsX(); ieta++){
    outf<<ieta<<" "<<rm.GetBinContent(ieta)<<"\t";
  }
  outf<<endl;
  return outf;
}

