#include "GEMHVScan/CosmicMuonAnalysis/plugins/GEMCosmicAnalyzer.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

GEMCosmicAnalyzer::GEMCosmicAnalyzer(const edm::ParameterSet& pset) {
  rechit_token_ = consumes<GEMRecHitCollection>(pset.getParameter<edm::InputTag>("recHitTag"));
  segment_token_ = consumes<GEMSegmentCollection>(pset.getParameter<edm::InputTag>("segmentTag"));

  muons_token_ = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonsTag"));
  muonsFromCosmics_token_ = consumes<edm::View<reco::Muon> >(pset.getParameter<edm::InputTag>("muonsFromCosmicsTag"));

  auto muon_service_parameter = pset.getParameter<edm::ParameterSet>("ServiceParameters");
  muon_service_ = new MuonServiceProxy(muon_service_parameter, consumesCollector());

  tree_muons_ = file_service_->make<TTree>("muons", "muons");
  tree_muonsFromCosmics_ = file_service_->make<TTree>("muonsFromCosmics", "muonsFromCosmics");

  makeBranch(tree_muons_);
  makeBranch(tree_muonsFromCosmics_);
}


void GEMCosmicAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup) {
  analyzeMuon(event, setup, tree_muons_, muons_token_);
  analyzeMuon(event, setup, tree_muonsFromCosmics_, muonsFromCosmics_token_);
}

void GEMCosmicAnalyzer::analyzeMuon(const edm::Event& event,
                                    const edm::EventSetup& setup,
                                    TTree* tree,
                                    const edm::EDGetTokenT<edm::View<reco::Muon> >& muon_token) {

  edm::Handle<GEMRecHitCollection> rechit_collection;
  event.getByToken(rechit_token_, rechit_collection);
  if (not rechit_collection.isValid()) {
    edm::LogError("GEMCosmicAnalyzer") << "GEMRecHitCollection is invalid" << std::endl;
    return;
  }

  edm::Handle<GEMSegmentCollection> segment_collection;
  event.getByToken(segment_token_, segment_collection);
  if (not segment_collection.isValid()) {
    edm::LogError("GEMCosmicAnalyzer") << "GEMSegmentCollection is invalid" << std::endl;
    return;
  }

  edm::Handle<edm::View<reco::Muon> > muon_view;
  event.getByToken(muon_token, muon_view);
  if (not muon_view.isValid()) {
    edm::LogError("GEMCosmicAnalyzer") << "View<Muon> is invalid" << std::endl;
    return;
  }

  edm::ESHandle<GEMGeometry> gem;
  setup.get<MuonGeometryRecord>().get(gem);
  if (not gem.isValid()) {
    edm::LogError("GEMCosmicAnalyzer") << "GEMGeometry is invalid" << std::endl;
    return;
  }

  edm::ESHandle<TransientTrackBuilder> transient_track_builder;
  setup.get<TransientTrackRecord>().get("TransientTrackBuilder", transient_track_builder);
  if (not transient_track_builder.isValid()) {
    edm::LogError("GEMCosmicAnalyzer") << "TransientTrackRecord is invalid" << std::endl;
    return;
  }

  muon_service_->update(setup);

  edm::ESHandle<Propagator>&& propagator_any = muon_service_->propagator("SteppingHelixPropagatorAny");
  if (not propagator_any.isValid()) {
    edm::LogError("GEMCosmicAnalyzer") << "SteppingHelixPropagatorAny is invalid" << std::endl;
    return;
  }

  edm::ESHandle<Propagator>&& propagator_along = muon_service_->propagator("SteppingHelixPropagatorAlong");
  if (not propagator_along.isValid()) {
    edm::LogError("GEMCosmicAnalyzer") << "SteppingHelixPropagatorAlong is invalid" << std::endl;
    return;
  }

  edm::ESHandle<Propagator>&& propagator_opposite = muon_service_->propagator("SteppingHelixPropagatorOpposite");
  if (not propagator_opposite.isValid()) {
    edm::LogError("GEMCosmicAnalyzer") << "SteppingHelixPropagatorOpposite is invalid" << std::endl;
    return;
  }

  // NOTE
  resetBranchEvent();

  b_event_ = event.id().event();
  b_lumi_block_ = event.luminosityBlock();
  b_run_ = event.run();

  for (size_t idx = 0; idx < muon_view->size(); ++idx) {
    const edm::RefToBase<reco::Muon>&& muon_ref = muon_view->refAt(idx);
    const reco::Muon* muon = muon_ref.get();

    const reco::Track* track = nullptr;
    if (muon->outerTrack().isNonnull()) {
      track = muon->outerTrack().get();
    }

    if (track == nullptr) {
      edm::LogError("GEMCosmicAnalyzer") << "failed to get muon track" << std::endl;
      continue;
    }

    const auto&& transient_track = transient_track_builder->build(track);
    if (not transient_track.isValid()) {
      edm::LogError("GEMCosmicAnalyzer") << "failed to build TransientTrack" << std::endl;
      continue;
    }

    // NOTE
    float x2_in = track->innerPosition().mag2();
    float x2_out = track->outerPosition().mag2();
    float p2_in = track->innerMomentum().mag2();
    float p2_out = track->outerMomentum().mag2();
    unsigned int raw_inner_det_id = track->innerDetId();
    unsigned int raw_outer_det_id = track->outerDetId();

    bool is_insideout = x2_in > x2_out;

    if (is_insideout) {
      std::swap(x2_in, x2_out);
      std::swap(p2_in, p2_out);
      std::swap(raw_inner_det_id, raw_outer_det_id);
    }

    bool is_incoming = p2_out > p2_in;

    // choose true innermost measurement state
    const auto&& start_state = is_insideout ? transient_track.outermostMeasurementState() : transient_track.innermostMeasurementState();
    auto& propagator = is_incoming ? propagator_along : propagator_opposite;
    // propagation_direction = PropagationDirection::anyDirection;

    // NOTE
    const auto&& impact_parameter = transient_track.stateAtBeamLine().transverseImpactParameter();
    const auto&& inner_id = getDetIdCandidate(raw_inner_det_id);
    const auto&& outer_id = getDetIdCandidate(raw_outer_det_id);

    resetBranchMuon();

    b_index_ = idx;
    b_is_incoming_ = is_incoming;
    b_pt_ = muon->pt();
    b_eta_ = muon->eta();
    b_phi_ = muon->phi();
    b_d0_ = impact_parameter.value();
    b_d0_error_ = impact_parameter.error();
    b_norm_chi2_ = track->normalizedChi2();
    b_ndof_ = track->ndof();
    b_num_valid_hits_ = track->numberOfValidHits();
    b_rechits_size_ = track->recHitsSize();
    b_num_seg_csc_st1_ = muon->numberOfSegments(1, MuonSubdetId::CSC);
    b_num_seg_csc_st2_ = muon->numberOfSegments(2, MuonSubdetId::CSC);
    b_num_seg_csc_st3_ = muon->numberOfSegments(3, MuonSubdetId::CSC);
    b_num_seg_csc_st4_ = muon->numberOfSegments(4, MuonSubdetId::CSC);
    b_has_seg_csc_ = (b_num_seg_csc_st1_ + b_num_seg_csc_st2_ + b_num_seg_csc_st3_ + b_num_seg_csc_st4_) > 0;

    b_in_detector_ = inner_id.detector;
    b_in_subdet_ = inner_id.subdet;
    b_in_wheel_ = inner_id.wheel;
    b_in_station_ = inner_id.station;
    b_in_sector_ = inner_id.sector;
    b_in_superlayer_ = inner_id.superlayer;
    b_in_layer_ = inner_id.layer;
    b_in_wire_ = inner_id.wire;
    b_in_region_ = inner_id.region;
    b_in_ring_ = inner_id.ring;
    b_in_chamber_ = inner_id.chamber;
    b_in_subsector_ = inner_id.subsector;
    b_in_roll_ = inner_id.roll;

    b_out_detector_ = outer_id.detector;
    b_out_subdet_ = outer_id.subdet;
    b_out_wheel_ = outer_id.wheel;
    b_out_station_ = outer_id.station;
    b_out_sector_ = outer_id.sector;
    b_out_superlayer_ = outer_id.superlayer;
    b_out_layer_ = outer_id.layer;
    b_out_wire_ = outer_id.wire;
    b_out_region_ = outer_id.region;
    b_out_ring_ = outer_id.ring;
    b_out_chamber_ = outer_id.chamber;
    b_out_subsector_ = outer_id.subsector;
    b_out_roll_ = outer_id.roll;

    b_start_x_err_ = std::sqrt(start_state.localError().positionError().xx());

    for (const GEMRegion* region : gem->regions()) {
      bool is_opposite_region = (muon->eta() * region->region() < 0);
      if (is_incoming xor is_opposite_region)
        continue;

      for (const GEMStation* station : region->stations()) {
        for (const GEMSuperChamber* super_chamber : station->superChambers()) {
          for (const GEMChamber* chamber : super_chamber->chambers()) {

            const BoundPlane& bound_plane = chamber->surface();
            const auto& [dest_state, path_length] = propagator->propagateWithPath(start_state, bound_plane);
            if (not dest_state.isValid()) {
              edm::LogInfo("GEMCosmicAnalyzer") << "failed to propagate ," << chamber->id() << std::endl;
              continue;
            }

            const GlobalPoint&& dest_global_pos = dest_state.globalPosition();
            const GEMEtaPartition* eta_partition = findEtaPartition(chamber, dest_global_pos);
            if (eta_partition == nullptr) {
              edm::LogInfo("GEMCosmicAnalyzer") << "failed to find GEMEtaPartition" << std::endl; 
              continue;
            }

            resetBranchChamber();

            // NOTE
            const GEMDetId&& gem_id = eta_partition->id();
            const auto rechit_range = rechit_collection->get(gem_id);
            const LocalPoint&& dest_local_pos = eta_partition->toLocal(dest_global_pos);
            const LocalError&& dest_local_err = dest_state.localError().positionError();

            b_region_ = gem_id.region();
            b_station_ = gem_id.station();
            b_layer_ = gem_id.layer();
            b_chamber_ = gem_id.chamber();
            b_ieta_ = gem_id.roll();

            b_path_length_ = path_length;
            b_num_same_ieta_hit_ = std::distance(rechit_range.first, rechit_range.second);

            if (auto hit = findClosetHit(dest_local_pos.x(), rechit_range)) {
              const LocalPoint&& hit_local_pos = hit->localPosition();
              const GlobalPoint&& hit_global_pos = eta_partition->toGlobal(hit_local_pos);
              const LocalError&& hit_local_err = hit->localPositionError();

              b_bx_ = hit->BunchX();
              b_strip_ = eta_partition->strip(hit_local_pos);
              b_first_cluster_strip_ = hit->firstClusterStrip();
              b_cluster_size_ = hit->clusterSize();

              b_residual_x_ = dest_local_pos.x() - hit_local_pos.x();
              b_residual_y_ = dest_local_pos.y() - hit_local_pos.y();
              b_residual_phi_ = reco::deltaPhi(dest_global_pos.barePhi(), hit_global_pos.barePhi());
              // TODO RdPhi

              b_pull_x_ = b_residual_x_ / std::sqrt(dest_local_err.xx() + hit_local_err.xx());
              b_pull_y_ = b_residual_y_ / std::sqrt(dest_local_err.yy() + hit_local_err.yy());

            }

            tree->Fill();

          } // GEMChamber
        } // GEMSuperChamber
      } // GEMStation
    } // GEMRegion
  } // reco::Muon

  // NOTE
}


const GEMEtaPartition* GEMCosmicAnalyzer::findEtaPartition(
    const GEMChamber* chamber, const GlobalPoint& global_point) {

  for (const GEMEtaPartition* eta_partition : chamber->etaPartitions()) {
    const LocalPoint&& local_point = eta_partition->toLocal(global_point);
    const LocalPoint local_point_2d(local_point.x(), local_point.y(), 0.0f);
    if (eta_partition->surface().bounds().inside(local_point_2d)) 
      return eta_partition;
  }

  return nullptr;
}


const GEMRecHit* GEMCosmicAnalyzer::findClosetHit(
    const float track_local_x, const GEMRecHitCollection::range& range) {

  float min_residual_x = std::numeric_limits<float>::infinity();;
  const GEMRecHit* closest_hit = nullptr;

  for (auto hit = range.first; hit != range.second; ++hit) {
    float residual_x = std::fabs(track_local_x - hit->localPosition().x());
    if (residual_x < min_residual_x) {
      min_residual_x = residual_x;
      closest_hit = &(*hit);
    }
  }

  return closest_hit;
}


GEMCosmicAnalyzer::DetIdCandidate GEMCosmicAnalyzer::getDetIdCandidate(unsigned int raw_det_id) {
  DetId det_id{raw_det_id};
  std::string label = "unknown";

  DetIdCandidate id;
  id.detector = static_cast<int>(det_id.det());

  if (det_id.det() == DetId::Detector::Muon) {
    id.subdet = static_cast<int>(det_id.subdetId());

    switch (det_id.subdetId()) {
      case MuonSubdetId::DT: {
        const DTWireId dt_id{det_id};
        id.wheel = dt_id.wheel();
        id.station = dt_id.station();
        id.sector = dt_id.sector();
        id.superlayer = dt_id.superLayer();
        id.layer = dt_id.layer();
        id.wire = dt_id.wire();
        break;
      }

      case MuonSubdetId::CSC: {
        const CSCDetId csc_id{det_id};
        id.region = static_cast<int>(csc_id.zendcap());
        id.station = csc_id.station();
        id.ring = csc_id.ring();
        id.chamber = csc_id.chamber();
        id.layer = csc_id.layer();
        break;
      }

      case MuonSubdetId::RPC: {
        const RPCDetId rpc_id{det_id};
        id.region = rpc_id.region();
        id.ring = rpc_id.ring();
        id.station = rpc_id.station();
        id.sector = rpc_id.sector();
        id.layer = rpc_id.layer();
        id.subsector = rpc_id.subsector();
        id.roll = rpc_id.roll();
        break;
      }

      case MuonSubdetId::GEM: {
        const GEMDetId gem_id{det_id};
        id.region = gem_id.region();
        id.ring = gem_id.ring();
        id.station = gem_id.station();
        id.layer = gem_id.layer();
        id.chamber = gem_id.chamber();
        id.roll = gem_id.roll();
        break;
      }

      default: {
        break;
      }
    } // switch (det_id.subdetId())
  }

  return id;
}


void GEMCosmicAnalyzer::makeBranch(TTree* tree) {
  // event
  tree->Branch("event", &b_event_);
  tree->Branch("lumi_block", &b_lumi_block_);
  tree->Branch("run", &b_run_);

  // muon
  tree->Branch("index", &b_index_);
  tree->Branch("is_incoming", &b_is_incoming_);
  tree->Branch("pt", &b_pt_);
  tree->Branch("eta", &b_eta_);
  tree->Branch("phi", &b_phi_);
  tree->Branch("d0", &b_d0_);
  tree->Branch("d0_error", &b_d0_error_);
  tree->Branch("norm_chi2", &b_norm_chi2_);
  tree->Branch("ndof", &b_ndof_);
  tree->Branch("num_valid_hits", &b_num_valid_hits_);
  tree->Branch("rechits_size", &b_rechits_size_);

  tree->Branch("num_seg_csc_st1", &b_num_seg_csc_st1_);
  tree->Branch("num_seg_csc_st2", &b_num_seg_csc_st2_);
  tree->Branch("num_seg_csc_st3", &b_num_seg_csc_st3_);
  tree->Branch("num_seg_csc_st4", &b_num_seg_csc_st4_);
  tree->Branch("has_seg_csc", &b_has_seg_csc_);
  tree->Branch("start_x_err", &b_start_x_err_);

  tree->Branch("in_detector", &b_in_detector_);
  tree->Branch("in_subdet", &b_in_subdet_);
  tree->Branch("in_wheel", &b_in_wheel_);
  tree->Branch("in_station", &b_in_station_);
  tree->Branch("in_sector", &b_in_sector_);
  tree->Branch("in_superlayer", &b_in_superlayer_);
  tree->Branch("in_layer", &b_in_layer_);
  tree->Branch("in_wire", &b_in_wire_);
  tree->Branch("in_region", &b_in_region_);
  tree->Branch("in_ring", &b_in_ring_);
  tree->Branch("in_chamber", &b_in_chamber_);
  tree->Branch("in_subsector", &b_in_subsector_);
  tree->Branch("in_roll", &b_in_roll_);

  tree->Branch("out_detector", &b_out_detector_);
  tree->Branch("out_subdet", &b_out_subdet_);
  tree->Branch("out_wheel", &b_out_wheel_);
  tree->Branch("out_station", &b_out_station_);
  tree->Branch("out_sector", &b_out_sector_);
  tree->Branch("out_superlayer", &b_out_superlayer_);
  tree->Branch("out_layer", &b_out_layer_);
  tree->Branch("out_wire", &b_out_wire_);
  tree->Branch("out_region", &b_out_region_);
  tree->Branch("out_ring", &b_out_ring_);
  tree->Branch("out_chamber", &b_out_chamber_);
  tree->Branch("out_subsector", &b_out_subsector_);
  tree->Branch("out_roll", &b_out_roll_);

  // NOTE destination
  tree->Branch("region", &b_region_);
  tree->Branch("station", &b_station_);
  tree->Branch("layer", &b_layer_);
  tree->Branch("chamber", &b_chamber_);
  tree->Branch("ieta", &b_ieta_);

  tree->Branch("path_length", &b_path_length_);
  tree->Branch("dest_state_x_err", &b_dest_x_err_);
  tree->Branch("num_same_ieta_hit", &b_num_same_ieta_hit_);

  tree->Branch("bx", b_bx_);
  tree->Branch("strip", b_strip_);
  tree->Branch("first_cluster_strip", &b_first_cluster_strip_);
  tree->Branch("cluster_size", &b_cluster_size_);

  tree->Branch("residual_x", &b_residual_x_);
  tree->Branch("residual_y", &b_residual_y_);
  tree->Branch("residual_phi", &b_residual_phi_);
  tree->Branch("pull_x", &b_pull_x_);
  tree->Branch("pull_y", &b_pull_y_);
}


void GEMCosmicAnalyzer::resetBranchEvent() {
  b_run_ = 0;
  b_lumi_block_ = 0;
  b_event_ = 0;
}


void GEMCosmicAnalyzer::resetBranchMuon() {
  b_index_ = -1;
  b_is_incoming_ = false;
  b_pt_ = -100.0f;
  b_eta_ = -100.0f;
  b_phi_ = -100.0f;
  b_d0_ = -100.0f;
  b_d0_error_ = -100.0f;
  b_norm_chi2_ = -100.0f;

  b_ndof_ = -1;
  b_num_valid_hits_ = -1;
  b_rechits_size_ = -1;
  b_num_seg_csc_st1_ = -1;
  b_num_seg_csc_st2_ = -1;
  b_num_seg_csc_st3_ = -1;
  b_num_seg_csc_st4_ = -1;
  b_has_seg_csc_ = false;
  b_start_x_err_ = -100.0f;
}


void GEMCosmicAnalyzer::resetBranchChamber() {
  b_region_ = 0;
  b_station_ = 0;
  b_layer_ = 0;
  b_chamber_ = 0;
  b_ieta_ = 0;

  b_path_length_ = -100.0f;
  b_dest_x_err_ = -100.0f;
  b_num_same_ieta_hit_ = -1;
  
  b_bx_ = -100;
  b_strip_ = -1.0f;
  b_first_cluster_strip_ = -1;
  b_cluster_size_ = -1;

  b_residual_x_ = -100.0f;
  b_residual_y_ = -100.0f;
  b_residual_phi_ = -100.0f;
  b_pull_x_ = -100.0f;
  b_pull_y_ = -100.0f;
}


void GEMCosmicAnalyzer::beginJob(){}
void GEMCosmicAnalyzer::endJob(){}

DEFINE_FWK_MODULE(GEMCosmicAnalyzer);
