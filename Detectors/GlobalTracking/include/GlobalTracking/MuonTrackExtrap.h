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

/// \file MatchGlobalFwd.h
/// \brief Class to perform MFT MCH (and MID) matching
/// \author rafael.pezzi@cern.ch

#ifndef ALICEO2_GLOBTRACKING_MUONEXTRAPTRACK_
#define ALICEO2_GLOBTRACKING_MUONEXTRAPTRACK_

#include <Rtypes.h>
#include <array>
#include <vector>
#include <string>
#include <gsl/span>
#include <TStopwatch.h>
#include "CommonConstants/LHCConstants.h"
#include "CommonUtils/ConfigurationMacroHelper.h"
#include "CommonDataFormat/BunchFilling.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "MFTTracking/IOUtils.h"
#include "MFTBase/Constants.h"
#include "DataFormatsMFT/TrackMFT.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "DataFormatsMCH/ROFRecord.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "ReconstructionDataFormats/GlobalFwdTrack.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/MatchInfoFwd.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "GlobalTracking/MatchGlobalFwdParam.h"
#include "DetectorsBase/GeometryManager.h"
#include "TGeoManager.h"
// #include <TGeoGlobalMagField.h>

// #include "MathUtils/Cartesian.h"
// #include "Field/MagneticField.h"
// #include "MCHBase/TrackBlock.h"

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
  // void setBz(float bz) { mBz = bz; }

  ///< set Bunch filling and init helpers for validation by BCs
  // void setBunchFilling(const o2::BunchFilling& bf);

  // void setMCTruthOn(bool v) { mMCTruthON = v; }

 private:
  // void fillBuiltinFunctions();

  bool prepareMCHData();
  bool processMCHMIDMatches();

  bool testProcessMCHMID();

  int extrapTracksToVertex(gsl::span<const o2::mch::TrackMCH> tracks, const math_utils::Point3D<double>& vertex);
  void writeTracks(char* bufferPtr) const;

  std::vector<std::vector<ExtrapMuonTrackStruct>> mTracksAtVtx{}; ///< list of tracks extrapolated to vertex for each event

  // float mBz = -5.f;                       ///< nominal Bz in kGauss
  // o2::InteractionRecord mStartIR{0, 0}; ///< IR corresponding to the start of the TF

  // o2::BunchFilling mBunchFilling;
  // std::array<int16_t, o2::constants::lhc::LHCMaxBunches> mClosestBunchAbove; // closest filled bunch from above
  // std::array<int16_t, o2::constants::lhc::LHCMaxBunches> mClosestBunchBelow; // closest filled bunch from below

  const o2::globaltracking::RecoContainer* mRecoCont = nullptr;
  gsl::span<const o2::mch::TrackMCH> mMCHTracks;                        ///< input MCH tracks
  gsl::span<const o2::mch::ROFRecord> mMCHTrackROFRec;                  ///< MCH tracks ROFRecords
  gsl::span<const o2::dataformats::TrackMCHMID> mMCHMIDMatches;         ///< input MCH MID Matches

  // std::vector<o2::math_utils::Bracket<float>> mMCHROFTimes;                          ///< min/max times of MCH ROFs in \mus
  // std::vector<int> mMCHID2Work;                                ///< MCH track id to ensure correct indexing for matching

  // bool mMCTruthON = false;      ///< Flag availability of MC truth
  // int mSaveMode = 0;            ///< Output mode [0 = SaveBestMatch; 1 = SaveAllMatches; 2 = SaveTrainingData]
  // TGeoManager* mGeoManager;
};

} // namespace globaltracking
} // namespace o2

#endif
