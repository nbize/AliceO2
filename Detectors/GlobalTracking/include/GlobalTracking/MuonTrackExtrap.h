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

/// \file MuonTrackExtrap.h
/// \brief Class to perform MCH-MID extrapolation from CTFs
/// \author nicolas.bize@cern.ch adapted from philippe.pillot@cern.ch and rafael.pezzi@cern.ch MCH and GlobalFwd workflows

#ifndef ALICEO2_GLOBTRACKING_MUONTRACKEXTRAP_
#define ALICEO2_GLOBTRACKING_MUONTRACKEXTRAP_

#include <Rtypes.h>
#include <array>
#include <vector>
#include <string>
#include <gsl/span>
#include <TStopwatch.h>
#include "CommonConstants/LHCConstants.h"
#include "CommonUtils/ConfigurationMacroHelper.h"
#include "CommonDataFormat/BunchFilling.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "MFTBase/Constants.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DetectorsBase/GeometryManager.h"
#include "TGeoManager.h"

namespace o2
{

namespace mch
{
class TrackMCH;
}

namespace globaltracking
{

struct ExtrapMuonTrackStruct {
  o2::mch::TrackParamStruct paramAtVertex{};
  double dca = 0.;
  double rAbs = 0.;
  int mchTrackIdx = 0;
};

class MuonTrackExtrap
{
 public:
  MuonTrackExtrap();
  ~MuonTrackExtrap() = default;

  void run(const o2::globaltracking::RecoContainer& inp);
  void init();
  void finalize();
  void clear();

  const std::vector<double>& getDCA() const { return mDCA; }
  const std::vector<double>& getDCAx() const { return mDCAx; }
  const std::vector<double>& getDCAy() const { return mDCAy; }
  const std::vector<double>& getP() const { return mP; }
  const std::vector<double>& getPt() const { return mPt; }
  // void setMCTruthOn(bool v) { mMCTruthON = v; }

 private:
  bool prepareMCHTracks();
  bool extrapMCHMIDTracks();

  std::vector<std::vector<ExtrapMuonTrackStruct>> mTracksAtVtx{}; ///< list of tracks extrapolated to vertex for each event
  std::vector<double> mDCA;
  std::vector<double> mDCAx;
  std::vector<double> mDCAy;
  std::vector<double> mP;
  std::vector<double> mPt;

  o2::InteractionRecord mStartIR{0, 0}; ///< IR corresponding to the start of the TF

  const o2::globaltracking::RecoContainer* mRecoCont = nullptr;
  gsl::span<const o2::mch::TrackMCH> mMCHTracks;                ///< input MCH tracks
  gsl::span<const o2::mch::ROFRecord> mMCHTrackROFRec;          ///< MCH tracks ROFRecords
  gsl::span<const o2::dataformats::TrackMCHMID> mMCHMIDMatches; ///< input MCH MID Matches

  std::vector<o2::math_utils::Bracket<float>> mMCHROFTimes; ///< min/max times of MCH ROFs in \mus
  std::vector<int> mMCHID2Work;                             ///< MCH track id to ensure correct indexing for matching
  std::vector<o2::mch::TrackMCH> mMCHFinalTracks;
};

} // namespace globaltracking
} // namespace o2

#endif
