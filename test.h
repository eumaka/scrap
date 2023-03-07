// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANE_EPFINDERRECO_H
#define EVENTPLANE_EPFINDERRECO_H

#include <fun4all/SubsysReco.h>

#include <string>

//Forward declarations
class PHCompositeNode;
class EpFinder;
class EpInfo;
class RawTowerContainer;
class RawTowerGeomContainer;
class PHG4HitContainer;
class TH1D;

class EpFinderReco : public SubsysReco
{
 public:
  EpFinderReco(const std::string &name = "EpFinderReco");
  ~EpFinderReco() override;

  int Init(PHCompositeNode *) override;

  int process_event(PHCompositeNode *) override;

  int End(PHCompositeNode *) override;

  int ResetEvent(PHCompositeNode * /*topNode*/) override;

  void set_algo_node(const std::string &algonode) { _algonode = algonode; }

  void set_ep_mode(int do_ep)
  {
    _do_ep = do_ep;
  }

  void Detector(const std::string &d)
  {
    detector = d;
  }

 private:
 
  TH1D *hcent = nullptr;
  
  void GetEventPlanes(PHCompositeNode *);
  int GetNodes(PHCompositeNode *);
  int CreateNodes(PHCompositeNode *);

  std::string _algonode = "EVENT_PLANE";
  int _do_ep = 0;

  std::string detector = "EPD";

  EpFinder *EpFinder_det[2] = nullptr;
  
  EpInfo *_EpInfo_det[2] = nullptr;


  std::string CaliTowerNodeName;
  std::string TowerGeomNodeName;
  std::string EPNodeName;
};

#endif  //* EVENTPLANE_EPFINDERRECO_H *//
