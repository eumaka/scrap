
#ifndef __Atestmod_H__
#define __Atestmod_H__

#include <fun4all/SubsysReco.h>
#include <string>
#include <vector>



//Forward declerations
class PHCompositeNode;
class TFile; 
class TTree;
class RawTowerContainer;
class RawTowerGeomContainer;
class TH2F;
class TProfile;
class SvtxTrackMap;

//class TowerInfoContainerv1;

//TowerInfoContainerv1* towers


//Brief: basic ntuple and histogram creation for sim evaluation
class Atestmod: public SubsysReco
{
 public: 
  //Default constructor
    Atestmod(const std::string &name="Atestmod");

  //Initialization, called for initialization
  int Init(PHCompositeNode *);

  int InitRun(PHCompositeNode *); 

  //Process Event, called for each event
  int process_event(PHCompositeNode *);

  //End, write and close files
  int EndRun(PHCompositeNode *);
  int End(PHCompositeNode *);

  //Change output filename
  void set_filename(const char* file)
  { if(file) _outfile_name = file; }

 private:
  //output filename
  std::string _outfile_name;
    
  //Event counter
  int _event;
  
  TTree* _event_tree;

  std::vector<int> _s;
  std::vector<int> _ns;
  std::vector<int> _t;
  std::vector<float> _x;
  std::vector<float> _y;
  std::vector<float> _z;
  std::vector<float> _e;
  
  
  std::vector<int> _sc;
  std::vector<int> _nsc;
  std::vector<int> _tc;
  std::vector<float> _xc;
  std::vector<float> _yc;
  std::vector<float> _zc;
  std::vector<float> _calibe;
  

  //User modules
  void fill_tree(PHCompositeNode*);

  //Get all the nodes
  int GetNodes(PHCompositeNode *);
 // unsigned int getchannel(unsigned int tower_key);
 

};

#endif //* __Atestmod_H__ *//
