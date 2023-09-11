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

/// \file ExtrapMuonTrackSpec.cxx.cxx
/// \brief Implementation of a data processor to extrapolate the tracks to the vertex
///
/// \author Nicolas Bizé, adapted from Philippe Pillot, Subatech and Rafael Pezzi work

#include <chrono>
#include <stdexcept>
#include <list>
#include <filesystem>

#include <TMath.h>

#include "CommonUtils/NameConf.h"

#include "Framework/ConfigParamRegistry.h"
#include "Framework/ControlService.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/Lifetime.h"
#include "Framework/Output.h"
#include "Framework/Task.h"
#include "Framework/CCDBParamSpec.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DataFormatsParameters/GRPObject.h"

#include "GlobalTracking/MuonTrackExtrap.h"
#include "GlobalTrackingWorkflow/ExtrapMuonTrackSpec.h"

// using namespace std;
using namespace o2::framework;
using GTrackID = o2::dataformats::GlobalTrackID;

namespace o2
{
namespace globaltracking
{

class ExtrapMuonTrackDPL : public Task
{
 public:
  ExtrapMuonTrackDPL(std::shared_ptr<DataRequest> dr, std::shared_ptr<o2::base::GRPGeomRequest> gr)
    : mDataRequest(dr), mGGCCDBRequest(gr) {}
  ~ExtrapMuonTrackDPL() override = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext& pc) final;
  void endOfStream(EndOfStreamContext& ec) final;
  void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj) final;

 private:
  void updateTimeDependentParams(ProcessingContext& pc);

  std::shared_ptr<DataRequest> mDataRequest;
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  o2::globaltracking::MuonTrackExtrap mExtrap; // Extrapolation engine

  std::chrono::duration<double> mElapsedTime{};
  TStopwatch mTimer; ///< timer
};

//_________________________________________________________________________________________________
void ExtrapMuonTrackDPL::init(InitContext& ic)
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
void ExtrapMuonTrackDPL::run(ProcessingContext& pc)
{
  mTimer.Start(false);
  RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest.get());
  updateTimeDependentParams(pc); // Make sure this is called after recoData.collectData, which may load some conditions

  mExtrap.run(recoData);

  // make outputs
  for (int i = 0; i < mExtrap.getDCA().size(); i++) {
    pc.outputs().snapshot(Output{"GLO", "DCA", 0, Lifetime::Timeframe}, mExtrap.getDCA()[i]);
    pc.outputs().snapshot(Output{"GLO", "DCAx", 0, Lifetime::Timeframe}, mExtrap.getDCAx()[i]);
    pc.outputs().snapshot(Output{"GLO", "DCAy", 0, Lifetime::Timeframe}, mExtrap.getDCAy()[i]);
    pc.outputs().snapshot(Output{"GLO", "p", 0, Lifetime::Timeframe}, mExtrap.getP()[i]);
    pc.outputs().snapshot(Output{"GLO", "pt", 0, Lifetime::Timeframe}, mExtrap.getPt()[i]);
  }

  mTimer.Stop();
}

void ExtrapMuonTrackDPL::endOfStream(EndOfStreamContext& ec)
{
  LOGF(info, "Extrapolation total timing: Cpu: %.3e Real: %.3e s in %d slots",
       mTimer.CpuTime(), mTimer.RealTime(), mTimer.Counter() - 1);
}

void ExtrapMuonTrackDPL::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    if (matcher == ConcreteDataMatcher("GLO", "GRPMAGFIELD", 0)) {
      o2::mch::TrackExtrap::setField();
    }
    return;
  }
}

void ExtrapMuonTrackDPL::updateTimeDependentParams(ProcessingContext& pc)
{
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  static bool initOnceDone = false;
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    mExtrap.init();
  }
  // we may have other params which need to be queried regularly
}

//_________________________________________________________________________________________________
o2::framework::DataProcessorSpec getExtrapMuonTrackSpec(const char* specName)
{
  // std::vector<InputSpec> inputs;
  std::vector<OutputSpec> outputs;

  auto dataRequest = std::make_shared<DataRequest>();
  o2::dataformats::GlobalTrackID::mask_t src = o2::dataformats::GlobalTrackID::getSourcesMask("MCH-MID");

  dataRequest->requestTracks(src, false);
  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              false,                             // GRPECS=true
                                                              false,                             // GRPLHCIF
                                                              true,                              // GRPMagField
                                                              false,                             // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              dataRequest->inputs,
                                                              true); // query only once all objects except mag.field

  outputs.emplace_back("GLO", "DCA", 0, Lifetime::Timeframe);
  outputs.emplace_back("GLO", "DCAx", 0, Lifetime::Timeframe);
  outputs.emplace_back("GLO", "DCAy", 0, Lifetime::Timeframe);
  outputs.emplace_back("GLO", "p", 0, Lifetime::Timeframe);
  outputs.emplace_back("GLO", "pt", 0, Lifetime::Timeframe);

  return DataProcessorSpec{
    specName,
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<ExtrapMuonTrackDPL>(dataRequest, ggRequest)},
    Options{
      // {"grp-file", VariantType::String, o2::base::NameConf::getGRPFileName(), {"Name of the grp file"}},
      // {"l3Current", VariantType::Float, -30000.0f, {"L3 current"}},
      // {"dipoleCurrent", VariantType::Float, -6000.0f, {"Dipole current"}}
    }};
}

} // namespace globaltracking
} // namespace o2
