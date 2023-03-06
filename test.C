#include "Atestmod.h"

#include <epd/EpdGeom.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>


#include <centrality/CentralityInfo.h>
#include <centrality/CentralityInfov1.h>

#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfov1.h>


#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Particle.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>


#include <phgeom/PHGeomUtility.h>

#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom_Spacalv1.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom_Spacalv3.h>
#include <g4detectors/PHG4CellDefs.h>

#include <TTree.h>
#include <TH2D.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <g4jets/Jet.h>
#include <g4jets/JetMap.h>
#include <g4jets/JetMapv1.h>
#include <g4main/PHG4Utils.h>
#include <epd/EPDDefs.h>

#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <boost/format.hpp>
#include <boost/math/special_functions/sign.hpp>

std::ofstream outfile1("channel.txt", std::ios::app);
std::ofstream outfile2("Calibtowers_channel.txt", std::ios::app);



using namespace std;

//----------------------------------------------------------------------------//
//-- Constructor:
//--  simple initialization
//----------------------------------------------------------------------------//

Atestmod::Atestmod(const string &name) :
  SubsysReco(name)
{
	//initialize
	_event = 0;
	_outfile_name = ".root";

}

//----------------------------------------------------------------------------//
//-- Init():
//--   Intialize all histograms, trees, and ntuples
//----------------------------------------------------------------------------//
int Atestmod::Init(PHCompositeNode *topNode) {

    std::cout << PHWHERE << " Opening file " << _outfile_name << std::endl;
    PHTFileServer::get().open(_outfile_name, "RECREATE");
    PHTFileServer::get().cd(_outfile_name);
    
    _event_tree = new TTree("event", "EPDQA => event info");
    _event_tree->Branch("event", &_event, "_event/I");
    _event_tree->Branch("cent", &_cent, "_cent/F");
    _event_tree->Branch("x", &_x);
    _event_tree->Branch("y", &_y);
    _event_tree->Branch("z", &_z);
    _event_tree->Branch("s", &_s);
    _event_tree->Branch("t", &_t);
    _event_tree->Branch("ns", &_ns);
    _event_tree->Branch("e", &_e);
    
    _event_tree->Branch("xc", &_xc);
    _event_tree->Branch("yc", &_yc);
    _event_tree->Branch("zc", &_zc);
    _event_tree->Branch("sc", &_sc);
    _event_tree->Branch("tc", &_tc);
    _event_tree->Branch("nsc", &_nsc);
    _event_tree->Branch("calibe", &_calibe);
 
    return Fun4AllReturnCodes::EVENT_OK;
}

int Atestmod::InitRun(PHCompositeNode *topNode)
{
	return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- process_event():
//--   Call user instructions for every event.
//--   This function contains the analysis structure.
//----------------------------------------------------------------------------//

int Atestmod::process_event(PHCompositeNode *topNode)
{
	_event++;

// 	GetNodes(topNode);

	fill_tree(topNode);

	return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- End():
//--   End method, wrap everything up
//----------------------------------------------------------------------------//

int Atestmod::EndRun(PHCompositeNode *topNode)
{

    return Fun4AllReturnCodes::EVENT_OK;

}

int Atestmod::End(PHCompositeNode *topNode)
{

	PHTFileServer::get().cd(_outfile_name);
	PHTFileServer::get().write(_outfile_name);

	return Fun4AllReturnCodes::EVENT_OK;
}

//----------------------------------------------------------------------------//
//-- fill_tree():
//--   Fill the various trees...
//----------------------------------------------------------------------------//

void Atestmod::fill_tree(PHCompositeNode *topNode)
{

    cout << _event << endl;
    
    bool dooutfile = "true";

    
    string towerinfonodename = "TOWERINFO_SIM_EPD";
    TowerInfoContainerv1 *_epd_towerinfos = findNode::getClass<TowerInfoContainerv1>(topNode, towerinfonodename.c_str());
    
    if (!_epd_towerinfos)
      {
          std::cout << "Could not locate SEPD tower info node " << std::endl;
          exit(1);
      }

     string towerinfonodename_calib = "TOWERINFO_CALIB_EPD";
    TowerInfoContainerv1 *_epd_towerinfos_calib = findNode::getClass<TowerInfoContainerv1>(topNode, towerinfonodename_calib.c_str());
    
    if (!_epd_towerinfos_calib)
      {
          std::cout << "Could not locate SEPD CALIB tower info node " << std::endl;
          exit(1);
      }
    
    string geominfonodename = "TOWERGEOM_EPD";
    EpdGeom *epdtilegeom = findNode::getClass<EpdGeom>(topNode, geominfonodename.c_str());
    
    if (!epdtilegeom)
      {
          std::cout << "Could not locate SEPD geometry node " << std::endl;
          exit(1);
      }
	
	
    CentralityInfov1 *cent = findNode::getClass<CentralityInfov1>(topNode, "CentralityInfo");
    if (!cent)
    {
           std::cout << " ERROR -- can't find CentralityInfo node" << std::endl;
           exit(1);
    }
    _cent = cent->get_centile(CentralityInfo::PROP::bimp);
    
    float tile_phi = 0.; float tile_r = 0.; float tile_z = 0.;
    std::tuple<unsigned int, unsigned int, unsigned int > simindex;
    int _arm = -1; int _sector = -1; int _tile = -1; float x = 0.;
    float y = 0.;
 
    
    unsigned int ntowers = _epd_towerinfos->size();
    for (unsigned int ch = 0; ch < ntowers;  ch++)
    {
           TowerInfo *raw_tower = _epd_towerinfos->get_tower_at_channel(ch);
           unsigned int thiskey =_epd_towerinfos->encode_epd(ch);
          
           tile_phi = epdtilegeom->phi(thiskey);
           tile_r = epdtilegeom->r(thiskey);
           tile_z = epdtilegeom->z(thiskey);
           simindex = epdtilegeom->id_to_side_sector_tile(thiskey); 
           _arm = get<0>(simindex);
           _sector = get<1>(simindex);
           _tile = get<2>(simindex);
	    
	    x = tile_r * cos(tile_phi);
            y = tile_r * sin(tile_phi);
           
           _x.push_back(x);
           _y.push_back(y);
           _z.push_back(tile_z);
           _e.push_back(raw_tower->get_energy()); 
           _s.push_back(_sector);
           _ns.push_back(_arm);
           _t.push_back(_tile);
           
           if(dooutfile)
           {
              if((raw_tower->get_energy() > 0.) && (_event < 10e3)) 
              outfile1 << _tile << "\t" << _sector << "\t" << _arm << "\t" << ch << "\t" << tile_phi << "\t"<< tile_r << "\t" << tile_z << "\t" << raw_tower->get_energy() << std::endl;
           }
      }



    float tile_phic = 0.; float tile_rc = 0.; float tile_zc = 0.;
    std::tuple<unsigned int, unsigned int, unsigned int > simindexc;
    int _armc = -1; int _sectorc = -1; int _tilec = -1; float xc = 0.;
    float yc =0.;
    
    unsigned int ntowers2 = _epd_towerinfos_calib->size();
    for (unsigned int ch2 = 0; ch2 < ntowers2;  ch2++)
    {
           TowerInfo *raw_tower2 = _epd_towerinfos_calib->get_tower_at_channel(ch2);
           unsigned int thiskey2 =_epd_towerinfos->encode_epd(ch2);
           
           tile_phic = epdtilegeom->phi(thiskey2);
           tile_rc = epdtilegeom->r(thiskey2);
           tile_zc = epdtilegeom->z(thiskey2);
           simindexc = epdtilegeom->id_to_side_sector_tile(thiskey2); 
           _armc = get<0>(simindexc);
           _sectorc = get<1>(simindexc);
           _tilec = get<2>(simindexc);
	    
	    xc = tile_rc * cos(tile_phic);
            yc = tile_rc * sin(tile_phic);
           
           _xc.push_back(x);
           _yc.push_back(y);
           _zc.push_back(tile_zc);
           _calibe.push_back(raw_tower2->get_energy()); 
           _sc.push_back(_sectorc);
           _nsc.push_back(_armc);
           _tc.push_back(_tilec);
           
            if(dooutfile)
            {
              if((raw_tower2->get_energy() > 0.) && (_event < 10e3)) 
              outfile2 << ch2 << "\t" << raw_tower2->get_energy() << std::endl;
            }
    }


_event_tree->Fill();

return;

}

//----------------------------------------------------------------------------//
//-- GetNodes():
//--   Get all the all the required nodes off the node tree
//----------------------------------------------------------------------------//
/*
int Atestmod::GetNodes(PHCompositeNode * topNode) {
    
    
//    string towerinfonodename = "TOWERINFO_SIM_EPD";
//    towers = findNode::getClass<TowerInfoContainerv1>(topNode, towerinfonodename.c_str());
//
//    if (!towers)
//     {
//        std::cout << PHWHERE << ": Could not find node " << towerinfonodename.c_str() << std::endl;
//        return Fun4AllReturnCodes::ABORTEVENT;
//     }

       
     return Fun4AllReturnCodes::EVENT_OK;
    
    
}

unsigned int Atestmod::getchannel(unsigned int tower_key)
{
  int channels_per_sector = -1;
  int supersector = -1;
  int nchannelsperpacket = -1;
  int maxphibin = -1;
  int maxetabin = -1;
  int etabinoffset[4] = {0};
  int phibinoffset[4] = {0};
  int etabinmap[4] = {0};

 
  channels_per_sector = 31;
  supersector = channels_per_sector * 12;
  unsigned int ns_sector = tower_key >> 20U;
  unsigned int rbin = (tower_key - (ns_sector << 20U)) >> 10U;
  unsigned int phibin = tower_key - (ns_sector << 20U) - (rbin << 10U);
  int epdchnlmap[16][2] = {{0, 0}, {1, 2}, {3, 4}, {5, 6}, {7, 8}, {9, 10}, {11, 12}, {13, 14}, {15, 16}, {17, 18}, {19, 20}, {21, 22}, {23, 24}, {25, 26}, {27, 28}, {29, 30}};
  int sector = phibin / 2;
  int channel = 0;
  if (rbin > 0)
  {
      channel = epdchnlmap[rbin][phibin - 2 * sector];
  }
  else
  {
      channel = 0;
  }
 
  unsigned int index = 0;
  index = ns_sector * supersector + sector * channels_per_sector + channel;
  return index;
}

*/
