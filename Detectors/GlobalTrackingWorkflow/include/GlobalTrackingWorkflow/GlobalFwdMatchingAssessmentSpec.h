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

#ifndef GLOFWD_ASSESSMENT_DEVICE_H
#define GLOFWD_ASSESSMENT_DEVICE_H

/// @file   GlobalFwdMatchingAssessmentSpec.h

#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "Framework/CCDBParamSpec.h"
#include "GlobalTracking/MatchGlobalFwdAssessment.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsParameters/GRPObject.h"
#include "TStopwatch.h"

using namespace o2::framework;

namespace o2
{
namespace globaltracking
{
class GlobalFwdAssessmentSpec : public Task
{
 public:
  GlobalFwdAssessmentSpec(std::shared_ptr<DataRequest> dr, std::shared_ptr<o2::base::GRPGeomRequest> gr, bool useMC, bool processGen, bool midFilterDisabled, bool finalizeAnalysis = false) : mDataRequest(dr),
                                                                                                                                                                                               mGGCCDBRequest(gr),
                                                                                                                                                                                               mUseMC(useMC),
                                                                                                                                                                                               mMIDFilterDisabled(midFilterDisabled),
                                                                                                                                                                                               mProcessGen(processGen),
                                                                                                                                                                                               mFinalizeAnalysis(finalizeAnalysis){};
  //~GlobalFwdAssessmentSpec() override = default;
  void init(o2::framework::InitContext& ic) final;
  void run(o2::framework::ProcessingContext& pc) final;
  void endOfStream(o2::framework::EndOfStreamContext& ec) final;
  void finaliseCCDB(ConcreteDataMatcher& matcher, void* obj) final;

 private:
  void sendOutput(DataAllocator& output);
  std::shared_ptr<DataRequest> mDataRequest;
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  std::unique_ptr<o2::globaltracking::GloFwdAssessment> mGloFwdAssessment;
  bool mUseMC = true;
  bool mProcessGen = false;
  bool mFinalizeAnalysis = false;
  bool mMIDFilterDisabled = false;
  enum TimerIDs { SWTot,
                  SWQCAsync,
                  SWTrackables,
                  SWGenerated,
                  SWRecoAndTrue,
                  SWAnalysis,
                  NStopWatches };
  static constexpr std::string_view TimerName[] = {"Total",
                                                   "ProcessAsync",
                                                   "ProcessTrackables",
                                                   "ProcessGenerated",
                                                   "ProcessRecoAndTrue",
                                                   "Analysis"};
  TStopwatch mTimer[NStopWatches];
};

DataProcessorSpec getGlobaFwdAssessmentSpec(bool useMC, bool processGen, bool midFilterDisabled, bool finalizeAnalysis = false);

} // namespace globaltracking
} // namespace o2

#endif