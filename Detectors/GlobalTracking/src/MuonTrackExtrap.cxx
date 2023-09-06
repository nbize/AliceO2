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

#include "GlobalTracking/MuonTrackExtrap.h"

using namespace o2::globaltracking;

//_________________________________________________________
void MuonTrackExtrap::init()
{
  LOG(info) << "Initializing Muon track extrapolation";
}

//_________________________________________________________
void MuonTrackExtrap::run(const o2::globaltracking::RecoContainer& inp)
{
  mRecoCont = &inp;

  if (!testProcessMCHMID()){
    return;
  }
  // mStartIR = inp.startIR;
  
  // clear();

  // if (!processMCHMIDMatches()) {
  //   return;
  // }

  // get the ROFs, tracks and vertices
    // auto rofs = pc.inputs().get<gsl::span<ROFRecord>>("rofs");
    // auto tracks = pc.inputs().get<gsl::span<TrackMCH>>("tracks");
    // auto vertices = pc.inputs().get<gsl::span<math_utils::Point3D<double>>>("vertices");

    // if (vertices.size() != rofs.size()) {
    //   throw length_error("number of vertices different from number of events");
    // }

    // for each event, propagate the tracks to the vertex
    // mTracksAtVtx.clear();
    // int iVertex(-1);
    // int nTracksTot(0);
    // for (const auto& rof : rofs) {
    //   auto tStart = std::chrono::high_resolution_clock::now();
    //   nTracksTot += extrapTracksToVertex(tracks.subspan(rof.getFirstIdx(), rof.getNEntries()), vertices[++iVertex]);
    //   auto tEnd = std::chrono::high_resolution_clock::now();
    //   mElapsedTime += tEnd - tStart;
    // }

    // // create the output message
    // auto msgOut = pc.outputs().make<char>(Output{"MCH", "TRACKSATVERTEX", 0, Lifetime::Timeframe},
    //                                       mTracksAtVtx.size() * sizeof(int) + nTracksTot * sizeof(TrackAtVtxStruct));

    // // write the tracks
    // writeTracks(msgOut.data());
}

void MuonTrackExtrap::finalize()
{
  LOG(info) << "Finalizing...";
}

//_________________________________________________________
void MuonTrackExtrap::clear()
{
  // mMCHROFTimes.clear();
  // mMCHWork.clear();
  // mMFTROFTimes.clear();
  // mMFTWork.clear();
  // mMFTClusters.clear();
  // mMatchedTracks.clear();
  // mMatchLabels.clear();
  // mMFTTrackROFContMapping.clear();
  // mMatchingInfo.clear();
}

bool MuonTrackExtrap::testProcessMCHMID()
{
  const auto& inp = *mRecoCont;
  // Test load MCH-MID tracks
  mMCHTracks = inp.getMCHTracks();
  mMCHTrackROFRec = inp.getMCHTracksROFRecords();
  mMCHMIDMatches = inp.getMCHMIDMatches();

  LOG(info) << "Number of MCH tracks : " << mMCHTracks.size();
  LOG(info) << "Number of MCH ROFs : " << mMCHTrackROFRec.size();
  LOG(info) << "Number of MCH-MID matches : " << mMCHMIDMatches.size();

  return true;
}

//_________________________________________________________
// bool MuonTrackExtrap::prepareMCHData()
// {
//   const auto& inp = *mRecoCont;

//   // Load MCH tracks
//   mMCHTracks = inp.getMCHTracks();
//   mMCHTrackROFRec = inp.getMCHTracksROFRecords();

//   int nROFs = mMCHTrackROFRec.size();
//   LOG(info) << "Loaded " << mMCHTracks.size() << " MCH Tracks in " << nROFs << " ROFs";
//   if (mMCHTracks.empty()) {
//     return false;
//   }
//   mMCHID2Work.clear();
//   mMCHID2Work.resize(mMCHTracks.size(), -1);
//   static int BCDiffErrCount = 0;
//   constexpr int MAXBCDiffErrCount = 2;

//   for (int irof = 0; irof < nROFs; irof++) {
//     const auto& rofRec = mMCHTrackROFRec[irof];

//     int nBC = rofRec.getBCData().differenceInBC(mStartIR);
//     if (nBC < 0) {
//       if (BCDiffErrCount++ < MAXBCDiffErrCount) {
//         LOGP(alarm, "wrong bunches diff. {} for current IR {} wrt 1st TF orbit {} in MCH data", nBC, rofRec.getBCData().asString(), mStartIR.asString());
//       }
//     }
//     float tMin = nBC * o2::constants::lhc::LHCBunchSpacingMUS;
//     float tMax = (nBC + rofRec.getBCWidth()) * o2::constants::lhc::LHCBunchSpacingMUS;
//     auto mchTime = rofRec.getTimeMUS(mStartIR).first;

//     mMCHROFTimes.emplace_back(tMin, tMax); // MCH ROF min/max time
//     LOG(debug) << "MCH ROF # " << irof << " " << rofRec.getBCData() << " [tMin;tMax] = [" << tMin << ";" << tMax << "]";
//     int trlim = rofRec.getFirstIdx() + rofRec.getNEntries();
//     for (int it = rofRec.getFirstIdx(); it < trlim; it++) {
//       auto& trcOrig = mMCHTracks[it];
//       // working copy MCH track propagated to matching plane and converted to the forward track format
//       o2::mch::TrackParam tempParam(trcOrig.getZ(), trcOrig.getParameters(), trcOrig.getCovariances());
//       if (!o2::mch::TrackExtrap::extrapToVertexWithoutBranson(tempParam, 0.)) {
//         LOG(warning) << "MCH track propagation to matching plane failed!";
//         continue;
//       }
//       // auto convertedTrack = MCHtoFwd(tempParam);
//       // auto& thisMCHTrack = mMCHWork.emplace_back(TrackLocMCH{convertedTrack, {tMin, tMax}});
//       // thisMCHTrack.setMCHTrackID(it);
//       // thisMCHTrack.setTimeMUS(mchTime);
//     }
//   }
//   return true;
// }

// //_________________________________________________________
// bool MuonTrackExtrap::processMCHMIDMatches()
// {

//     const auto& inp = *mRecoCont;

//     // Load MCHMID matches
//     mMCHMIDMatches = inp.getMCHMIDMatches();

//     LOG(info) << "Loaded " << mMCHMIDMatches.size() << " MCHMID matches";

//     for (const auto& MIDMatch : mMCHMIDMatches) {
//       const auto& MCHId = MIDMatch.getMCHRef().getIndex();
//       const auto& MIDId = MIDMatch.getMIDRef().getIndex();
//       // auto& thisMuonTrack = mMCHWork[mMCHID2Work[MCHId]];
//       // LOG(debug) << " MCHId: " << MCHId << " --> mMCHID2Work[MCHId]:" << mMCHID2Work[MCHId];
//       // const auto& IR = MIDMatch.getIR();
//       // int nBC = IR.differenceInBC(mStartIR);
//       // float tMin = nBC * o2::constants::lhc::LHCBunchSpacingMUS;
//       // float tMax = (nBC + 1) * o2::constants::lhc::LHCBunchSpacingMUS;
//       // thisMuonTrack.setMIDTrackID(MIDId);
//       // thisMuonTrack.setTimeMUS(MIDMatch.getTimeMUS(mStartIR).first);
//       // thisMuonTrack.tBracket.set(tMin, tMax);
//       // thisMuonTrack.setMIDMatchingChi2(MIDMatch.getMatchChi2OverNDF());
    
//   }
//   return true;
// }

//_________________________________________________________________________________________________
int MuonTrackExtrap::extrapTracksToVertex(gsl::span<const o2::mch::TrackMCH> tracks, const math_utils::Point3D<double>& vertex)
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

//_________________________________________________________________________________________________
void MuonTrackExtrap::writeTracks(char* bufferPtr) const
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

MuonTrackExtrap::MuonTrackExtrap()
{
  // Do nothing
}