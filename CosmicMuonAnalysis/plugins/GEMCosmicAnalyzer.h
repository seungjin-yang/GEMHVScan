#ifndef DQMOffline_Muon_GEMCosmicAnalyzer_h
#define DQMOffline_Muon_GEMCosmicAnalyzer_h


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Framework/interface/ESHandle.h"


#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"


class GEMCosmicAnalyzer : public edm::EDAnalyzer {
public:
  explicit GEMCosmicAnalyzer(const edm::ParameterSet&);
  ~GEMCosmicAnalyzer(){};

private:
  // NOTE DT: wheel, station, sector, superlayer, layer, wire
  // NOTE CSCDetId(int iendcap, int istation, int iring, int ichamber, int ilayer)
  // NOTE RPCDetId(int region, int ring, int station, int sector, int layer, int subsector, int roll);
  // NOTE GEMDetId(int region, int ring, int station, int layer, int chamber, int roll)
  struct DetIdCandidate {
    int detector;
    int subdet;
    int wheel;
    int station;
    int sector;
    int superlayer;
    int layer;
    int wire;
    int region; // iendcap
    int ring;
    int chamber;
    int subsector;
    int roll;

    DetIdCandidate() : detector(0), subdet(0), wheel(-10), station(-10),
                       sector(-10), superlayer(-10), layer(-10), wire(-10),
                       region(-10), ring(-10), chamber(-10), subsector(-10),
                       roll(-10) {}

  };



  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();

  // NOTE
  void resetBranchEvent();
  void resetBranchMuon();
  void resetBranchChamber();

  const GEMEtaPartition*  findEtaPartition(const GEMChamber*, const GlobalPoint&);
  const GEMRecHit* findClosetHit(const float, const GEMRecHitCollection::range&);
  DetIdCandidate getDetIdCandidate(unsigned int);

  // NOTE
  MuonServiceProxy *muon_service_;
  edm::EDGetTokenT<GEMRecHitCollection> rechit_token_;
  edm::EDGetTokenT<GEMSegmentCollection> segment_token_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muon_token_;

  edm::Service<TFileService> file_service_;
  TTree* tree_;
  TH1F* h_count_;
  TH2F* h_rechit_neg_; // GE-1/1
  TH2F* h_rechit_pos_; // GE+1/1
  TH1F* h_rechit_;

  // event
  unsigned int b_run_;
  unsigned int b_lumi_block_;
  unsigned long long b_event_;

  // muon
  int b_index_;
  bool b_is_incoming_;
  float b_pt_;
  float b_eta_;
  float b_phi_;
  float b_d0_;
  float b_d0_error_;
  float b_norm_chi2_;
  int b_ndof_;
  int b_num_valid_hits_;
  int b_rechits_size_;
  int b_num_seg_csc_st1_;
  int b_num_seg_csc_st2_;
  int b_num_seg_csc_st3_;
  int b_num_seg_csc_st4_;
  int b_has_seg_csc_;
  float b_start_x_err_;

  int b_in_detector_;
  int b_in_subdet_;
  int b_in_wheel_;
  int b_in_station_;
  int b_in_sector_;
  int b_in_superlayer_;
  int b_in_layer_;
  int b_in_wire_;
  int b_in_region_;
  int b_in_ring_;
  int b_in_chamber_;
  int b_in_subsector_;
  int b_in_roll_;

  int b_out_detector_;
  int b_out_subdet_;
  int b_out_wheel_;
  int b_out_station_;
  int b_out_sector_;
  int b_out_superlayer_;
  int b_out_layer_;
  int b_out_wire_;
  int b_out_region_;
  int b_out_ring_;
  int b_out_chamber_;
  int b_out_subsector_;
  int b_out_roll_;

  // destination
  int b_region_;
  int b_station_;
  int b_layer_;
  int b_chamber_;
  int b_ieta_;
  float b_path_length_;
  float b_dest_x_err_;
  int b_num_same_ieta_hit_;
  // 
  int b_bx_;
  int b_strip_;
  float b_first_cluster_strip_;
  int b_cluster_size_;
  float b_residual_x_;
  float b_residual_y_;
  float b_residual_phi_;
  float b_pull_x_;
  float b_pull_y_;
};


#endif  // DQMOffline_Muon_GEMCosmicAnalyzer_h
