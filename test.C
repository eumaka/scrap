#include "EpFinderReco.h"

#include "EpFinder.h"
#include "EpInfo.h"    // for EpInfo
#include "EpInfov1.h"  // for EpInfo

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <epd/EpdGeom.h>
#include <centrality/CentralityInfo.h>
#include <centrality/CentralityInfov1.h>

#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfov1.h>


#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>  // for TH2F
#include <TSystem.h>
#include <TVector3.h>  // for TVector3

#include <algorithm>  // for max
#include <cmath>      // for M_PI
#include <cstdlib>    // for NULL, exit, getenv
#include <iostream>
#include <map>  // for _Rb_tree_const_iterator
#include <utility>
#include <vector>  // for vector

TH2F *mapCalib = nullptr;

using namespace std;

EpFinderReco::EpFinderReco(const std::string &name)
  : SubsysReco(name)
  , detector("EPD")
{
}

EpFinderReco::~EpFinderReco()
{
  for(int i = 0; i < 2; i++)
  {
    delete EpFinder_det[i];
  }
}

int EpFinderReco::Init(PHCompositeNode *topNode)
{
	
  for(int i = 0; i < 2; i++)
  {
    EpFinder_det[i] = new EpFinder(1, 3);
  }
  
    return CreateNodes(topNode);
}

int EpFinderReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *AlgoNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", _algonode));
  if (!AlgoNode)
  {
    AlgoNode = new PHCompositeNode(_algonode);
    dstNode->addNode(AlgoNode);
  }

   for(int i = 0; i < 2; i++)
   {
     EpInfo *EpInfo_det = new EpInfov1();
     PHIODataNode<PHObject> *EpInfo_det_node = new PHIODataNode<PHObject>(EpInfo_det,Form("EpInfo_det%i",i), "PHObject");
     AlgoNode->addNode(EpInfo_det_node);
   }

    return Fun4AllReturnCodes::EVENT_OK;
}

int EpFinderReco::process_event(PHCompositeNode *topNode)
{

  GetNodes(topNode);
  GetEventPlanes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int EpFinderReco::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int EpFinderReco::ResetEvent(PHCompositeNode * /*topNode*/)
{
	
  for(int i = 0; i < 2; i++)
  {
    EpFinder_det[i]->ResetEvent();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//

void EpFinderReco::GetEventPlanes(PHCompositeNode *topNode)
{
  
      CentralityInfov1 *cent = findNode::getClass<CentralityInfov1>(topNode, "CentralityInfo");
      if (!cent)
      {
        std::cout << " ERROR -- can't find CentralityInfo node" << std::endl;
        return;
      }

      int _b = cent->get_centile(CentralityInfo::PROP::epd_NS)/10;
	     
     TowerInfoContainerv1 *_epd_towerinfos_calib = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_EPD");
     if (!_epd_towerinfos_calib)
      {
          std::cout << "Could not locate SEPD CALIB tower info node " << std::endl;
          exit(1);
      }

     EpdGeom *epdtilegeom = findNode::getClass<EpdGeom>(topNode,"TOWERGEOM_EPD");
     if (!epdtilegeom)
      {
          std::cout << "Could not locate SEPD geometry node " << std::endl;
          exit(1);
      }


      std::vector<EpHit> tnehits;
      tnehits.clear();

      std::vector<EpHit> tsehits;
      tsehits.clear();
      
      float tile_phi = 0.; float tile_z = 0.; float tile_e = 0.;
      int eta_bin = -1;
      unsigned int ntowers = _epd_towerinfos_calib->size();
      for (unsigned int ch = 0; ch < ntowers;  ch++)
      {
         TowerInfo *raw_tower = _epd_towerinfos_calib->get_tower_at_channel(ch);
         unsigned int thiskey =_epd_towerinfos_calib->encode_epd(ch);
	 tile_phi = epdtilegeom->phi(thiskey);
	 tile_z = epdtilegeom->z(thiskey);  
	 eta_bin = raw_tower->getTowerEtaBin(thiskey);
	      
	 //truncate tile energies
	 if ((raw_tower->get_energy()) < 0.2) continue; 
	 std::cout<<eta_bin<<"\t"<<_b<<std::endl;
//          tile_e = (raw_tower->get_energy() < mapCalib->GetBinContent(_b + 1, i + 1)) ? raw_tower->get_energy() : mapCalib->GetBinContent(_b + 1, i + 1);
	      
	 if(tile_z > 0)
	 {
           EpHit newHit;
           newHit.nMip = tile_e;
           newHit.phi = tile_phi;
           tnehits.push_back(newHit);
	 }
	 else if(tile_z < 0)
	  {
           EpHit newHit;
           newHit.nMip = tile_e;
           newHit.phi = tile_phi;
           tsehits.push_back(newHit);
	 }
	      
      }//end loop over sepd tower info
	
 
     EpFinder_det[0]->Results(tsehits, 0, _EpInfo_det[0]);
     EpFinder_det[1]->Results(tnehits, 0, _EpInfo_det[1]);
	
     tsehits.clear(); 
     tnehits.clear();
 
  return;
}

int EpFinderReco::GetNodes(PHCompositeNode *topNode)
{
 
  for(int i=0; i<2; i++)
  {
     _EpInfo_det[i] = findNode::getClass<EpInfo>(topNode,Form("EpInfo_det%i",i));
     if (!_EpInfo_det[i]) 
     {
	  std::cout << PHWHERE << ": Could not find node:"<< "_EpInfo_det" <<i<< std::endl;
          return Fun4AllReturnCodes::ABORTEVENT;
     }
  }
 
  return Fun4AllReturnCodes::EVENT_OK;
	
}




