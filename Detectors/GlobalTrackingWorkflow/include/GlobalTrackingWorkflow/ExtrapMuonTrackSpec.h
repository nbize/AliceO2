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

/// \file TrackAtVertexSpec.h
/// \brief Definition of a data processor to extrapolate the tracks to the vertex
///
/// \author Philippe Pillot, Subatech

#ifndef O2_MCHMID_TRACKEXTRAPSPEC_H_
#define O2_MCHMID_TRACKEXTRAPSPEC_H_

#include "Framework/DataProcessorSpec.h"
// #include "DetectorsBase/GRPGeomHelper.h"
// #include "DataFormatsGlobalTracking/RecoContainer.h"

namespace o2
{
namespace globaltracking
{
o2::framework::DataProcessorSpec getExtrapMuonTrackSpec(const char* specName = "mchmid-track-extrap");

} // end namespace globaltracking
} // end namespace o2

#endif // O2_MCHMIDTRACKEXTRAPSPEC_H_
