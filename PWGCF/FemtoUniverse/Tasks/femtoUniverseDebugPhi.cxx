// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoUniverseDebugPhi.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for Phis
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw, zuzanna.chochulska.stud@pw.edu.pl

#include <fairlogger/Logger.h>
#include <cstdint>
#include <iostream>
#include <vector>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "DataFormatsParameters/GRPObject.h"

#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseParticleHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseEventHisto.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUtils.h"

using namespace o2;
using namespace o2::analysis::femtoUniverse;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct femtoUniverseDebugPhi {
  SliceCache cache;

  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodePhi{"ConfPDGCodePartOne", 3122, "Phi - PDG code"};
    Configurable<uint32_t> ConfCutPhi{"ConfCutPhi", 338, "Phi - Selection bit from cutCulator"};
    ConfigurableAxis ConfPhiTempFitVarBins{"ConfPhiTempFitVarBins", {300, 0.95, 1.}, "Phi: binning of the TempFitVar in the pT vs. TempFitVar plot"};
    ConfigurableAxis ConfPhiTempFitVarpTBins{"ConfPhiTempFitVarpTBins", {20, 0.5, 4.05}, "Phi: pT binning of the pT vs. TempFitVar plot"};
  } ConfPhigroup;

  struct : o2::framework::ConfigurableGroup {
    Configurable<int> ConfPDGCodeChildPos{"ConfPDGCodeChildPos", 2212, "Positive Child - PDG code"};
    Configurable<int> ConfPDGCodeChildNeg{"ConfPDGCodeChildNeg", 211, "Negative Child- PDG code"};
    Configurable<uint32_t> ConfCutChildPos{"ConfCutChildPos", 150, "Positive Child of Phi - Selection bit from cutCulator"};
    Configurable<uint32_t> ConfCutChildNeg{"ConfCutChildNeg", 149, "Negative Child of Phi - Selection bit from cutCulator"};
    Configurable<float> ConfChildPosPidnSigmaMax{"ConfChildPosPidnSigmaMax", 3.f, "Positive Child of Phi - Selection bit from cutCulator"};
    Configurable<float> ConfChildNegPidnSigmaMax{"ConfChildNegPidnSigmaMax", 3.f, "Negative Child of Phi - Selection bit from cutCulator"};
    Configurable<int> ConfChildPosIndex{"ConfChildPosIndex", 1, "Positive Child of Phi - Index from cutCulator"};
    Configurable<int> ConfChildNegIndex{"ConfChildNegIndex", 0, "Negative Child of Phi - Index from cutCulator"};
    Configurable<std::vector<float>> ConfChildPIDnSigmaMax{"ConfChildPIDnSigmaMax", std::vector<float>{4.f, 3.f}, "Phi child sel: Max. PID nSigma TPC"};
    Configurable<int> ConfChildnSpecies{"ConfChildnSpecies", 2, "Number of particle spieces (for Phi children) with PID info"};
    ConfigurableAxis ConfChildTempFitVarBins{"ConfChildTempFitVarBins", {300, -0.15, 0.15}, "Phi child: binning of the TempFitVar in the pT vs. TempFitVar plot"};
    ConfigurableAxis ConfChildTempFitVarpTBins{"ConfChildTempFitVarpTBins", {20, 0.5, 4.05}, "Phi child: pT binning of the pT vs. TempFitVar plot"};
  } ConfChildGroup;

  using FemtoFullParticles = soa::Join<aod::FDParticles, aod::FDExtParticles>;
  Partition<FemtoFullParticles> partsOne = (aod::femtouniverseparticle::partType == uint8_t(aod::femtouniverseparticle::ParticleType::kPhi)); // && ((aod::femtouniverseparticle::cut & ConfPhigroup.ConfCutPhi) == ConfPhigroup.ConfCutPhi);
  Preslice<FemtoFullParticles> perCol = aod::femtouniverseparticle::fdCollisionId;

  /// Histogramming
  FemtoUniverseEventHisto eventHisto;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kPhiChild, 3> posChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kPhiChild, 4> negChildHistos;
  FemtoUniverseParticleHisto<aod::femtouniverseparticle::ParticleType::kPhi> PhiHistos;

  /// Histogram output
  HistogramRegistry EventRegistry{"Event", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry PhiRegistry{"FullPhiQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    eventHisto.init(&EventRegistry);
    posChildHistos.init(&PhiRegistry, ConfChildGroup.ConfChildTempFitVarpTBins, ConfChildGroup.ConfChildTempFitVarBins, false, ConfChildGroup.ConfPDGCodeChildPos.value, true);
    negChildHistos.init(&PhiRegistry, ConfChildGroup.ConfChildTempFitVarpTBins, ConfChildGroup.ConfChildTempFitVarBins, false, ConfChildGroup.ConfPDGCodeChildNeg, true);
    PhiHistos.init(&PhiRegistry, ConfPhigroup.ConfPhiTempFitVarpTBins, ConfPhigroup.ConfPhiTempFitVarBins, false, ConfPhigroup.ConfPDGCodePhi.value, true);
  }

  /// Porduce QA plots for Phi selection in FemtoUniverse framework
  void process(o2::aod::FDCollision const& col, FemtoFullParticles const& parts)
  {
    auto groupPartsOne = partsOne->sliceByCached(aod::femtouniverseparticle::fdCollisionId, col.globalIndex(), cache);
    eventHisto.fillQA(col);
    for (auto& part : groupPartsOne) {
      if (!part.has_children()) {
        continue;
      }
      const auto& posChild = parts.iteratorAt(part.index() - 2);
      const auto& negChild = parts.iteratorAt(part.index() - 1);
      if (posChild.globalIndex() != part.childrenIds()[0] || negChild.globalIndex() != part.childrenIds()[1]) {
        LOG(warn) << "Indices of Phi children do not match";
        continue;
      }
      // check cuts on Phi children
      if ((posChild.partType() == uint8_t(aod::femtouniverseparticle::ParticleType::kPhiChild) /*&& (posChild.cut() & ConfChildGroup.ConfCutChildPos) == ConfChildGroup.ConfCutChildPos*/) &&
          (negChild.partType() == uint8_t(aod::femtouniverseparticle::ParticleType::kPhiChild) /*&& (negChild.cut() & ConfChildGroup.ConfCutChildNeg) == ConfChildGroup.ConfCutChildNeg*/) //&&
          // isFullPIDSelected(posChild.pidcut(), posChild.p(), 999.f, ConfChildGroup.ConfChildPosIndex.value, ConfChildGroup.ConfChildnSpecies.value, ConfChildGroup.ConfChildPIDnSigmaMax.value, ConfChildGroup.ConfChildPosPidnSigmaMax.value, 1.f) &&
          // isFullPIDSelected(negChild.pidcut(), negChild.p(), 999.f, ConfChildGroup.ConfChildNegIndex.value, ConfChildGroup.ConfChildnSpecies.value, ConfChildGroup.ConfChildPIDnSigmaMax.value, ConfChildGroup.ConfChildNegPidnSigmaMax.value, 1.f)
      ) {
        PhiHistos.fillQA<false, true>(part);
        posChildHistos.fillQA<false, true>(posChild);
        negChildHistos.fillQA<false, true>(negChild);
      }
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<femtoUniverseDebugPhi>(cfgc),
  };
  return workflow;
}
