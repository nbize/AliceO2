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

/// @file   GlobalFwdTrackWriterSpec.cxx

#include <vector>
#include "GlobalTrackingWorkflow/ExtrapMuonWriterSpec.h"
#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "GlobalTracking/MuonTrackExtrap.h"

using namespace o2::framework;

namespace o2
{
namespace globaltracking
{

template <typename T>
using BranchDefinition = MakeRootTreeWriterSpec::BranchDefinition<T>;

DataProcessorSpec getExtrapMuonWriterSpec()
{
  return MakeRootTreeWriterSpec("extrap-muon-writer",
                                "muonTrackExtrapTrueMeanVertex.root",
                                "TrackInfos",
                                BranchDefinition<double>{InputSpec{"dca", "GLO", "DCA", 0}, "dca"},
                                BranchDefinition<double>{InputSpec{"dcax", "GLO", "DCAx", 0}, "dcax"},
                                BranchDefinition<double>{InputSpec{"dcay", "GLO", "DCAy", 0}, "dcay"},
                                BranchDefinition<double>{InputSpec{"p", "GLO", "p", 0}, "p"},
                                BranchDefinition<double>{InputSpec{"pt", "GLO", "pt", 0}, "pt"},
                                BranchDefinition<double>{InputSpec{"ptOrig", "GLO", "ptOrig", 0}, "ptOrig"},
                                BranchDefinition<double>{InputSpec{"rabs", "GLO", "rabs", 0}, "rabs"},
                                BranchDefinition<double>{InputSpec{"x", "GLO", "x", 0}, "x"},
                                BranchDefinition<double>{InputSpec{"y", "GLO", "y", 0}, "y"},
                                BranchDefinition<double>{InputSpec{"z", "GLO", "z", 0}, "z"},
                                BranchDefinition<double>{InputSpec{"xAtDCA", "GLO", "xAtDCA", 0}, "xAtDCA"},
                                BranchDefinition<double>{InputSpec{"yAtDCA", "GLO", "yAtDCA", 0}, "yAtDCA"},
                                BranchDefinition<double>{InputSpec{"zAtDCA", "GLO", "zAtDCA", 0}, "zAtDCA"})();
}

} // namespace globaltracking
} // namespace o2
