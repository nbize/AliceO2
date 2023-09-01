// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TrackAtVertexSpec.cxx
/// \brief Implementation of a data processor to extrapolate the tracks to the vertex
///
/// \author Philippe Pillot, Subatech

#include "GlobalTrackingWorkflow/ExtrapMuonTrackSpec.h"

#include <chrono>
#include <stdexcept>
#include <list>

#include <gsl/span>
#include <filesystem>

#include <TMath.h>
#include <TGeoManager.h>
#include <TGeoGlobalMagField.h>

#include "CommonUtils/NameConf.h"

#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Output.h"
#include "Framework/Task.h"
#include "Framework/CCDBParamSpec.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"

#include "MathUtils/Cartesian.h"
#include "Field/MagneticField.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "MCHBase/TrackBlock.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackExtrap.h"

namespace o2
{
namespace globaltracking
{

struct ExtrapMuonTrackStruct {
  o2::mch::TrackParamStruct paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  int mchTrackIdx = 0;
};

using namespace std;
using namespace o2::framework;

class ExtrapMuonTrackTask
{
 public:
 ExtrapMuonTrackTask(std::shared_ptr<globaltracking::DataRequest> dr, std::shared_ptr<o2::base::GRPGeomRequest> gr)
    : mDataRequest(dr), mGGCCDBRequest(gr) {}
  //~ExtrapMuonTrackTask() override = default;
  //_________________________________________________________________________________________________
  void init(framework::InitContext& ic)
  {
    /// Prepare the track extrapolation tools
    LOG(info) << "initializing track extrapolation to vertex";
    o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);

    auto stop = [this]() {
      LOG(info) << "track propagation to vertex duration = " << mElapsedTime.count() << " s";
    };
    ic.services().get<CallbackService>().set<CallbackService::Id::Stop>(stop);
  }

  //_________________________________________________________________________________________________
  void run(framework::ProcessingContext& pc)
  {
    /// propagate the MCH tracks to the vertex for each event in the TF and send the results

    // get the ROFs, tracks and vertices
    auto rofs = pc.inputs().get<gsl::span<o2::mch::ROFRecord>>("rofs");
    auto tracks = pc.inputs().get<gsl::span<o2::mch::TrackMCH>>("tracks");
    // auto vertices = pc.inputs().get<gsl::span<math_utils::Point3D<double>>>("vertices");
    math_utils::Point3D<double> vertices = {0., 0., 0.}; // For now define the vertex at 0
    // get the data
    o2::globaltracking::RecoContainer recoData;
    recoData.collectData(pc, *mDataRequest.get());
    o2::base::GRPGeomHelper::instance().checkUpdates(pc); 

    // if (vertices.size() != rofs.size()) { // to be added once we will get the vertices from ITS
    //   throw length_error("number of vertices different from number of events");
    // }

    // for each event, propagate the tracks to the vertex
    mTracksAtVtx.clear();
    int nTracksTot(0);
    for (const auto& rof : rofs) {
      auto tStart = std::chrono::high_resolution_clock::now();
      nTracksTot += extrapTracksToVertex(tracks.subspan(rof.getFirstIdx(), rof.getNEntries()), vertices);
      auto tEnd = std::chrono::high_resolution_clock::now();
      mElapsedTime += tEnd - tStart;
    }

    // create the output message
    auto msgOut = pc.outputs().make<char>(Output{"MCH", "TRACKSATVERTEX", 0, Lifetime::Timeframe},
                                          mTracksAtVtx.size() * sizeof(int) + nTracksTot * sizeof(ExtrapMuonTrackStruct));

    // write the tracks
    writeTracks(msgOut.data());
  }

 private:
  std::shared_ptr<o2::globaltracking::DataRequest> mDataRequest;
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  const o2::globaltracking::RecoContainer* mRecoCont = nullptr;
  gsl::span<const o2::dataformats::TrackMCHMID> mMCHMIDMatches;         ///< input MCH MID Matches
  //_________________________________________________________________________________________________
  int extrapTracksToVertex(gsl::span<const o2::mch::TrackMCH> tracks, const math_utils::Point3D<double>& vertex)
  {
    /// compute the tracks parameters at vertex, at DCA and at the end of the absorber
    /// return the number of tracks successfully propagated to the vertex

    auto& tracksAtVtx = mTracksAtVtx.emplace_back();
    int trackIdx(-1);

    for (const auto& track : tracks) {

      // create a new track at vertex pointing to the current track (index within the current event)
      auto& trackAtVtx = tracksAtVtx.emplace_back();
      trackAtVtx.mchTrackIdx = ++trackIdx;

      // extrapolate to vertex
      o2::mch::TrackParam trackParamAtVertex(track.getZ(), track.getParameters());
      if (!o2::mch::TrackExtrap::extrapToVertex(trackParamAtVertex, vertex.x(), vertex.y(), vertex.z(), 0., 0.)) {
        tracksAtVtx.pop_back();
        continue;
      }
      trackAtVtx.paramAtVertex.x = trackParamAtVertex.getNonBendingCoor();
      trackAtVtx.paramAtVertex.y = trackParamAtVertex.getBendingCoor();
      trackAtVtx.paramAtVertex.z = trackParamAtVertex.getZ();
      trackAtVtx.paramAtVertex.px = trackParamAtVertex.px();
      trackAtVtx.paramAtVertex.py = trackParamAtVertex.py();
      trackAtVtx.paramAtVertex.pz = trackParamAtVertex.pz();
      trackAtVtx.paramAtVertex.sign = trackParamAtVertex.getCharge();

      // extrapolate to DCA
      o2::mch::TrackParam trackParamAtDCA(track.getZ(), track.getParameters());
      if (!o2::mch::TrackExtrap::extrapToVertexWithoutBranson(trackParamAtDCA, vertex.z())) {
        tracksAtVtx.pop_back();
        continue;
      }
      double dcaX = trackParamAtDCA.getNonBendingCoor() - vertex.x();
      double dcaY = trackParamAtDCA.getBendingCoor() - vertex.y();
      trackAtVtx.dca = TMath::Sqrt(dcaX * dcaX + dcaY * dcaY);
      LOG(info) << "DCA calculation : ";
      LOG(info) << "dcaX = " << dcaX;
      LOG(info) << "dcaY = " << dcaY;
      LOG(info) << "dca = " << trackAtVtx.dca;

      // extrapolate to the end of the absorber
      o2::mch::TrackParam trackParamAtRAbs(track.getZ(), track.getParameters());
      if (!o2::mch::TrackExtrap::extrapToZ(trackParamAtRAbs, -505.)) {
        tracksAtVtx.pop_back();
        continue;
      }
      double xAbs = trackParamAtRAbs.getNonBendingCoor();
      double yAbs = trackParamAtRAbs.getBendingCoor();
      trackAtVtx.rAbs = TMath::Sqrt(xAbs * xAbs + yAbs * yAbs);
    }

    return tracksAtVtx.size();
  }

  bool processMCHMIDmatches(){
    // const auto& inp = *mRecoCont;

    // Load MCHMID matches
    // mMCHMIDMatches = inp.getMCHMIDMatches();
    // LOG(info) << "Loaded " << mMCHMIDMatches.size() << " MCHMID matches";
  }

  //_________________________________________________________________________________________________
  void writeTracks(char* bufferPtr) const
  {
    /// write the track informations for each event in the message payload

    for (const auto& tracksAtVtx : mTracksAtVtx) {
      // write the number of tracks
      int nTracks = tracksAtVtx.size();
      memcpy(bufferPtr, &nTracks, sizeof(int));
      bufferPtr += sizeof(int);

      // write the tracks
      if (nTracks > 0) {
        memcpy(bufferPtr, tracksAtVtx.data(), nTracks * sizeof(ExtrapMuonTrackStruct));
        bufferPtr += nTracks * sizeof(ExtrapMuonTrackStruct);
      }
    }
  }

  std::vector<std::vector<ExtrapMuonTrackStruct>> mTracksAtVtx{}; ///< list of tracks extrapolated to vertex for each event
  std::chrono::duration<double> mElapsedTime{};              ///< timer
};

void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    if (matcher == ConcreteDataMatcher("GLO", "GRPMAGFIELD", 0)) {
      o2::mch::TrackExtrap::setField();
    }
    return;
  }
}

//_________________________________________________________________________________________________
o2::framework::DataProcessorSpec getExtrapMuonTrackSpec(const char* specName)
{
  std::vector<InputSpec> inputs;
  inputs.emplace_back("rofs", "MCH", "TRACKROFS", 0, Lifetime::Timeframe);
  inputs.emplace_back("tracks", "MCH", "TRACKS", 0, Lifetime::Timeframe);
  inputs.emplace_back("clusters", "MCH", "TRACKCLUSTERS", 0, Lifetime::Timeframe);
  auto dataRequest = std::make_shared<o2::globaltracking::DataRequest>();
  //dataRequest->requestMCHMIDMatches(false); // Request MCHMID Matches. Labels are not used
  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              true,                              // GRPECS=true
                                                              true,                              // GRPLHCIF
                                                              true,                              // GRPMagField
                                                              false,                             // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              inputs,
                                                              true); // query only once all objects except mag.field
  return DataProcessorSpec{
    specName,
    inputs,
    Outputs{OutputSpec{{"tracksatvertex"}, "MCH", "TRACKSATVERTEX", 0, Lifetime::Timeframe}},
    AlgorithmSpec{adaptFromTask<ExtrapMuonTrackTask>(dataRequest,ggRequest)},
    Options{
      {"grp-file", VariantType::String, o2::base::NameConf::getGRPFileName(), {"Name of the grp file"}},
      // {"l3Current", VariantType::Float, -30000.0f, {"L3 current"}},
      // {"dipoleCurrent", VariantType::Float, -6000.0f, {"Dipole current"}}
      }};
}

} // namespace globaltracking
} // namespace o2
