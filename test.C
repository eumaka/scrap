#include "EpFinderReco.h"

#include "EpFinder.h"
#include "EpInfo.h"    // for EpInfo
#include "EpInfov1.h"  // for EpInfo

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <centrality/CentralityInfo.h>
#include <centrality/CentralityInfov1.h>

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
  , detector("CEMC")
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
  

  if (detector == "EPD")
  {
      double centrange[11] = { 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100.};
      hcent = new TH1D("hcent", "cent hist", 10, centrange);
     
      std::string sepdmapname;
      const char *Calibroot = getenv("CALIBRATIONROOT");
      if (Calibroot)
      {
        sepdmapname = Calibroot;
      }
      else
      {
        std::cout << "no CALIBRATIONROOT environment variable" << std::endl;
        gSystem->Exit(1);
      }

      sepdmapname += "/EPD/Calibmap/sEPDCalibMap.root";

      TFile *file = new TFile(sepdmapname.c_str());
      mapCalib = (TH2F *) file->Get("epdcalib");
      if (!mapCalib)
      {
        std::cout << "ERROR: mapCalib is NULL" << std::endl;
        gSystem->Exit(1);
      }
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
  
  else if (_do_ep == 3)
  {
    std::vector<EpHit> nehits;
    nehits.clear();

    std::vector<EpHit> sehits;
    sehits.clear();

    int _side = -1;
    int phi_bin = -1;
    int eta_bin = -1;
    int num_phi = 24;
    int num_eta = 16;
    float meanPhi = 0.;
    float nmips = 0.;

    double _epd_e[2][16][24] = {{{0.0}}};
    double _epd_e_calib[2][16][24] = {{{0.0}}};

    PHG4HitContainer::ConstRange e_range = e_hit_container->getHits();
    for (PHG4HitContainer::ConstIterator e_itr = e_range.first; e_itr != e_range.second; e_itr++)
    {
      PHG4Hit *ehit = e_itr->second;
      if (!ehit)
      {
        continue;
      }
   
      if (ehit->get_light_yield() < 0.0) continue;
      if ((e_itr->second->get_t(0) > -50) && (e_itr->second->get_t(1) < 50))
      {
        TVector3 ehitPos(ehit->get_avg_x(), ehit->get_avg_y(), ehit->get_avg_z()); 
        if (ehit->get_z(0) > 0)
        {
          _side = 1;
        }
        else
        {
          _side = 0;
        }

        if (_side == 1)
        {
          //do tile segmentation

          if ((ehitPos.Eta() >= 2.01) && (ehitPos.Eta() < 4.29))
          {
            phi_bin = GetPhiBin(ehitPos.Phi(), num_phi);
            eta_bin = heta->FindBin(ehitPos.Eta()) - 1;
          }

          if ((ehitPos.Eta() >= 4.29) && (ehitPos.Eta() < 4.96))
          {
            phi_bin = GetPhiBin(ehitPos.Phi(), 0.5 * num_phi);
            eta_bin = heta->FindBin(ehitPos.Eta()) - 1;
          }
        }

        // south wheel
        else
        {
          if ((ehitPos.Eta() <= -2.01) && (ehitPos.Eta() > -4.29))
          {
            phi_bin = GetPhiBin(ehitPos.Phi(), num_phi);
            eta_bin = heta->FindBin(-1 * ehitPos.Eta()) - 1;
          }

          if ((ehitPos.Eta() <= -4.29) && (ehitPos.Eta() > -4.96))
          {
            phi_bin = GetPhiBin(ehitPos.Phi(), 0.5 * num_phi);
            eta_bin = heta->FindBin(-1 * ehitPos.Eta()) - 1;
          }
        }

        if ((eta_bin >= 0) && (phi_bin >= 0)) _epd_e[_side][eta_bin][phi_bin] += ehit->get_light_yield();

        if (_do_sepd_calib)
        {
          nmips = ehit->get_light_yield() / _sepdmpv;
          if ((eta_bin >= 0) && (phi_bin >= 0)) _epd_e_calib[_side][eta_bin][phi_bin] += nmips;
        }
      }
    }

    for (int i = 0; i < num_eta; i++)
    {
      for (int j = 0; j < num_phi; j++)
      {
        if (_epd_e[1][i][j] > 0.0)
        {
          meanPhi = GetMeanPhi(j, num_phi);
          if (i == 15) meanPhi = GetMeanPhi(j, 0.5 * num_phi);

          EpHit newHit;
          newHit.nMip = _epd_e[1][i][j];
          newHit.phi = meanPhi;
          nehits.push_back(newHit);
        }
      }
    }

    EpFinder_1->Results(nehits, 0, _EPD_EpInfoN);
    nehits.clear();

    for (int i = 0; i < num_eta; i++)
    {
      for (int j = 0; j < num_phi; j++)
      {
        if (_epd_e[0][i][j] > 0.0)
        {
          meanPhi = GetMeanPhi(j, num_phi);
          if (i == 15) meanPhi = GetMeanPhi(j, 0.5 * num_phi);

          EpHit newHit;
          newHit.nMip = _epd_e[0][i][j];
          newHit.phi = meanPhi;
          sehits.push_back(newHit);
        }
      }
    }

    EpFinder_2->Results(sehits, 0, _EPD_EpInfoS);
    sehits.clear();

    if (_do_sepd_calib)
    {
      if (Verbosity() >= 1) std::cout << "using centrality dependent Nmip truncation" << std::endl;

      CentralityInfov1 *cent = findNode::getClass<CentralityInfov1>(topNode, "CentralityInfo");
      if (!cent)
      {
        std::cout << " ERROR -- can't find CentralityInfo node" << std::endl;
        return;
      }

      int _b = hcent->FindBin(cent->get_centile(CentralityInfo::PROP::epd_NS)) - 1;

      std::vector<EpHit> tnehits;
      tnehits.clear();

      std::vector<EpHit> tsehits;
      tsehits.clear();

      for (int i = 0; i < num_eta; i++)
      {
        for (int j = 0; j < num_phi; j++)
        {
          if (_epd_e_calib[1][i][j] > 0.0)
          {
            meanPhi = GetMeanPhi(j, num_phi);
            if (i == 15) meanPhi = GetMeanPhi(j, 0.5 * num_phi);

            if (_epd_e_calib[1][i][j] < 0.2) continue;

            float calib_e = (_epd_e_calib[1][i][j] < mapCalib->GetBinContent(_b + 1, i + 1)) ? _epd_e_calib[1][i][j] : mapCalib->GetBinContent(_b + 1, i + 1);

            EpHit newHit;
            newHit.nMip = calib_e;
            newHit.phi = meanPhi;
            tnehits.push_back(newHit);
          }
        }
      }

      EpFinder_3->Results(tnehits, 0, _EPD_EpInfoN_calib);
      tnehits.clear();

      for (int i = 0; i < num_eta; i++)
      {
        for (int j = 0; j < num_phi; j++)
        {
          if (_epd_e_calib[0][i][j] > 0.0)
          {
            meanPhi = GetMeanPhi(j, num_phi);
            if (i == 15) meanPhi = GetMeanPhi(j, 0.5 * num_phi);

            if (_epd_e_calib[0][i][j] < 0.2) continue;

            float calib_e = (_epd_e_calib[0][i][j] < mapCalib->GetBinContent(_b + 1, i + 1)) ? _epd_e_calib[0][i][j] : mapCalib->GetBinContent(_b + 1, i + 1);

            EpHit newHit;
            newHit.nMip = calib_e;
            newHit.phi = meanPhi;
            tsehits.push_back(newHit);
          }
        }
      }

      EpFinder_4->Results(tsehits, 0, _EPD_EpInfoS_calib);
      tsehits.clear();
    }
  }

  return;
}

int EpFinderReco::GetNodes(PHCompositeNode *topNode)
{
  if (_do_ep == 0)
  {
    if (detector == "CEMC" || detector == "HCALOUT" || detector == "HCALIN")
    {
      //EpInfo nodes
      EPNodeName = "EPINFO_" + detector;
      _CALO_EpInfo = findNode::getClass<EpInfo>(topNode, EPNodeName);
      if (!_CALO_EpInfo)
      {
        std::cout << PHWHERE << ": Could not find node: " << EPNodeName << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }

      // Detector towers nodes
      CaliTowerNodeName = "TOWER_CALIB_" + detector;
      _calib_towers = findNode::getClass<RawTowerContainer>(topNode, CaliTowerNodeName);

      if (!_calib_towers)
      {
        std::cout << PHWHERE << ": Could not find node: " << CaliTowerNodeName << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }

      // Detector geometry nodes
      TowerGeomNodeName = "TOWERGEOM_" + detector;
      rawtowergeom = findNode::getClass<RawTowerGeomContainer>(topNode, TowerGeomNodeName);
      if (!rawtowergeom)
      {
        std::cout << PHWHERE << ": Could not find node: " << TowerGeomNodeName << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
    }
    else
    {
      std::cout << PHWHERE
                << " Detector choice does not match, use a single calorimeter for this mode, exiting"
                << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
  }

  else if (_do_ep == 1)
  {
    _CEMCHCAL_EpInfo = findNode::getClass<EpInfo>(topNode, "EPINFO_CEMCHCAL");
    if (!_CEMCHCAL_EpInfo)
    {
      std::cout << PHWHERE << " EPINFO_CEMCHCAL node not found on node tree" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    cemctowers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
    if (!cemctowers)
    {
      std::cout << PHWHERE << ": Could not find node TOWER_CALIB_CEMC " << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    hcalotowers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");
    if (!hcalotowers)
    {
      std::cout << PHWHERE << ": Could not find node TOWER_CALIB_HCALOUT" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    hcalitowers = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
    if (!hcalitowers)
    {
      std::cout << PHWHERE << ": Could not find node TOWER_CALIB_HCALIN " << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    cemctowergeom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if (!cemctowergeom)
    {
      std::cout << PHWHERE << ": Could not find node TOWERGEOM_CEMC" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    hcalotowergeom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    if (!hcalotowergeom)
    {
      std::cout << PHWHERE << ": Could not find node TOWERGEOM_HCALOUT" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    hcalitowergeom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    if (!hcalitowergeom)
    {
      std::cout << PHWHERE << ": Could not find node TOWERGEOM_HCALIN " << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  else if (_do_ep == 2)
  {
    if (detector != "BBC")
    {
      std::cout << PHWHERE
                << " Detector choice does not match, use BBC for this mode, exiting"
                << std::endl;
      gSystem->Exit(1);
      exit(1);
    }

    else
    {
      //EpInfo nodes
      _BBC_EpInfoN = findNode::getClass<EpInfo>(topNode, "EPINFO_BBC_North");
      if (!_BBC_EpInfoN)
      {
        std::cout << PHWHERE << ": Could not find node: EPINFO_BBC_North" << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }

      _BBC_EpInfoS = findNode::getClass<EpInfo>(topNode, "EPINFO_BBC_South");
      if (!_BBC_EpInfoN)
      {
        std::cout << PHWHERE << ": Could not find node: EPINFO_BBC_South" << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      // G4 Hit nodes
      b_hit_container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_" + detector);
      if (!b_hit_container)
      {
        std::cout << PHWHERE << ": Could not find node: G4HIT_BBC " << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
    }
  }

 
  for(int i=0; i<2; i++)
  {

    _CEMC_EpInfo_EtaSlice[i] = findNode::getClass<EpInfo>(topNode,Form("EpInfo_det%i",i));
    if (!_CEMC_EpInfo_EtaSlice[i]) {
      cout << PHWHERE << " _CEMC_EpInfo_EtaSlice"<<i<<" node not found on node tree"
           << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }


 
        _EPD_EpInfoN_calib = findNode::getClass<EpInfo>(topNode, "EPINFO_EPD_North_calib");
        if (!_EPD_EpInfoN_calib)
        {
          std::cout << PHWHERE << ": Could not find node: EPINFO_EPD_North_calib" << std::endl;
          return Fun4AllReturnCodes::ABORTEVENT;
        }

        _EPD_EpInfoS_calib = findNode::getClass<EpInfo>(topNode, "EPINFO_EPD_South_calib");
        if (!_EPD_EpInfoS_calib)
        {
          std::cout << PHWHERE << ": Could not find node: EPINFO_EPD_South_calib" << std::endl;
          return Fun4AllReturnCodes::ABORTEVENT;
        }
      }

  }
  return Fun4AllReturnCodes::EVENT_OK;
}




