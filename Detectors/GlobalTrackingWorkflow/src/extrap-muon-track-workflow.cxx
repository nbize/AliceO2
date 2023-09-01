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

/// \file tracks-to-tracks-at-vertex-workflow.cxx
/// \brief Implementation of a DPL device to run the algorithm that extrapolates tracks to the vertex
///
/// \author Philippe Pillot, Subatech

#include "CommonUtils/ConfigurableParam.h"
#include "GlobalTrackingWorkflow/ExtrapMuonTrackSpec.h"

using namespace o2::framework;

// we need to add workflow options before including Framework/runDataProcessing
void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  // option allowing to set parameters
  std::vector<o2::framework::ConfigParamSpec> options{
    {"configKeyValues", VariantType::String, "", {"Semicolon separated key=value strings"}}};

  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

WorkflowSpec defineDataProcessing(const ConfigContext& configcontext)
{
  // auto dataRequest = std::make_shared<globaltracking::DataRequest>();

  // auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
  //                                                             true,                              // GRPECS=true
  //                                                             true,                              // GRPLHCIF
  //                                                             true,                              // GRPMagField
  //                                                             false,                             // askMatLUT
  //                                                             o2::base::GRPGeomRequest::Aligned, // geometry
  //                                                             inputs,
  //                                                             true); // query only once all objects except mag.field
  o2::conf::ConfigurableParam::updateFromString(configcontext.options().get<std::string>("configKeyValues"));
  return WorkflowSpec{o2::mch::getExtrapMuonTrackSpec("mch-extrap-muon-track")};
}
