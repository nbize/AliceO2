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
  mStartIR = inp.startIR;

  clear();

  if (!prepareMCHTracks()){
    return;
  }
  if (!extrapMCHMIDTracks()){
    return;
  }
}

void MuonTrackExtrap::finalize()
{
  LOG(info) << "Finalizing...";
}

//_________________________________________________________
void MuonTrackExtrap::clear()
{
  mMCHROFTimes.clear();
  mMCHFinalTracks.clear();
  mTracksAtVtx.clear();
  mDCA.clear();
  mDCAx.clear();
  mDCAy.clear();
  mP.clear();
  mPt.clear();
}

bool MuonTrackExtrap::extrapMCHMIDTracks()
{
  const auto& inp = *mRecoCont;
  // Load MCH-MID tracks
  mMCHMIDMatches = inp.getMCHMIDMatches();

  math_utils::Point3D<double> vertex = {0., 0., 0.}; // For now define the vertex at 0
  auto& tracksAtVtx = mTracksAtVtx.emplace_back();
  auto& vecDCA = mDCA;
  auto& vecDCAx = mDCAx;
  auto& vecDCAy = mDCAy;
  auto& vecP = mP;
  auto& vecPt = mPt;
  int trackIdx(-1);

  LOG(info) << "Number of MCH-MID matches : " << mMCHMIDMatches.size();

  for (const auto& MIDMatch : mMCHMIDMatches) {
    const auto& MCHId = MIDMatch.getMCHRef().getIndex();
    const auto& MIDId = MIDMatch.getMIDRef().getIndex();

    auto& thisMuonTrack = mMCHFinalTracks[mMCHID2Work[MCHId]];
    LOG(debug) << " MCHId: " << MCHId << " --> mMCHID2Work[MCHId]:" << mMCHID2Work[MCHId];

    LOG(debug) << "Muon track Z = " << thisMuonTrack.getZ();

    // create a new track at vertex pointing to the current track (index within the current event)
    auto& trackAtVtx = tracksAtVtx.emplace_back();
    trackAtVtx.mchTrackIdx = ++trackIdx;

    // extrapolate to vertex
    o2::mch::TrackParam trackParamAtVertex(thisMuonTrack.getZ(), thisMuonTrack.getParameters());
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

    LOG(debug) << "Muon track after extrapolation, Z = " << trackAtVtx.paramAtVertex.z;

    // extrapolate to DCA
    o2::mch::TrackParam trackParamAtDCA(thisMuonTrack.getZ(), thisMuonTrack.getParameters());
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
    double p = trackParamAtVertex.p();
    double px = trackParamAtVertex.px();
    double py = trackParamAtVertex.py();
    double pt = TMath::Sqrt(px * px + py * py);
   
    vecDCA.emplace_back(trackAtVtx.dca);
    vecDCAx.emplace_back(dcaX);
    vecDCAy.emplace_back(dcaY);
    vecP.emplace_back(p);
    vecPt.emplace_back(pt);

    // extrapolate to the end of the absorber
    o2::mch::TrackParam trackParamAtRAbs(thisMuonTrack.getZ(), thisMuonTrack.getParameters());
    if (!o2::mch::TrackExtrap::extrapToZ(trackParamAtRAbs, -505.)) {
      tracksAtVtx.pop_back();
      continue;
    }
    double xAbs = trackParamAtRAbs.getNonBendingCoor();
    double yAbs = trackParamAtRAbs.getBendingCoor();
    trackAtVtx.rAbs = TMath::Sqrt(xAbs * xAbs + yAbs * yAbs);
  }

  return true;
}

//_________________________________________________________
bool MuonTrackExtrap::prepareMCHTracks()
{
  const auto& inp = *mRecoCont;

  // Load MCH tracks
  mMCHTracks = inp.getMCHTracks();
  mMCHTrackROFRec = inp.getMCHTracksROFRecords();

  int nROFs = mMCHTrackROFRec.size();
  LOG(info) << "Loaded " << mMCHTracks.size() << " MCH Tracks in " << nROFs << " ROFs";
  if (mMCHTracks.empty()) {
    return false;
  }
  mMCHFinalTracks.reserve(mMCHTracks.size());
  mMCHID2Work.clear();
  mMCHID2Work.resize(mMCHTracks.size(), -1);
  static int BCDiffErrCount = 0;
  constexpr int MAXBCDiffErrCount = 2;

  for (int irof = 0; irof < nROFs; irof++) {
    const auto& rofRec = mMCHTrackROFRec[irof];

    int nBC = rofRec.getBCData().differenceInBC(mStartIR);
    if (nBC < 0) {
      if (BCDiffErrCount++ < MAXBCDiffErrCount) {
        LOGP(alarm, "wrong bunches diff. {} for current IR {} wrt 1st TF orbit {} in MCH data", nBC, rofRec.getBCData().asString(), mStartIR.asString());
      }
    }
    float tMin = nBC * o2::constants::lhc::LHCBunchSpacingMUS;
    float tMax = (nBC + rofRec.getBCWidth()) * o2::constants::lhc::LHCBunchSpacingMUS;
    auto mchTime = rofRec.getTimeMUS(mStartIR).first;

    mMCHROFTimes.emplace_back(tMin, tMax); // MCH ROF min/max time
    LOG(debug) << "MCH ROF # " << irof << " " << rofRec.getBCData() << " [tMin;tMax] = [" << tMin << ";" << tMax << "]";
    int trlim = rofRec.getFirstIdx() + rofRec.getNEntries();
    for (int it = rofRec.getFirstIdx(); it < trlim; it++) {
      auto& trcOrig = mMCHTracks[it];
      int nWorkTracks = mMCHFinalTracks.size();
      mMCHID2Work[it] = nWorkTracks;
      auto& mchTrack = mMCHFinalTracks.emplace_back(trcOrig);
    }
  }
  return true;
}

MuonTrackExtrap::MuonTrackExtrap()
{
  // Do nothing
}