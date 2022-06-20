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

///
/// \file   TrackletsTask.cxx
/// \author My Name
///

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <sstream>
#include <string>
#include <tuple>
#include "QualityControl/QcInfoLogger.h"
#include "TRD/TrackletsTask.h"
#include <Framework/InputRecord.h>
#include <Framework/InputRecordWalker.h>
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/NoiseCalibration.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "CCDB/BasicCCDBManager.h"
//#include "TRDQC/StatusHelper.h"

namespace o2::quality_control_modules::trd
{

TrackletsTask::~TrackletsTask()
{
}

void TrackletsTask::drawLinesMCM(TH2F* histo, Double_t len = 40, Int_t m = 1)
{

  TLine* l;
  Int_t nPos[o2::trd::constants::NSTACK - 1] = { 16, 32, 44, 60 };

  for (Int_t iStack = 0; iStack < o2::trd::constants::NSTACK - 1; ++iStack) {
    l = new TLine(nPos[iStack] - 0.5, -0.5, nPos[iStack] - 0.5, len);
    l->SetLineStyle(2);
    histo->GetListOfFunctions()->Add(l);
  }

  for (Int_t iLayer = 0; iLayer < o2::trd::constants::NLAYER * m; ++iLayer) {
    l = new TLine(-0.5, iLayer * 8 - 0.5, 75.5, iLayer * 8 - 0.5);
    l = new TLine(0.5, iLayer * 8 - 0.5, 75.5, iLayer * 8 - 0.5);
    l->SetLineStyle(2);
    histo->GetListOfFunctions()->Add(l);
  }
}

void TrackletsTask::retrieveCCDBSettings()
{
  // std::string a = mCustomParameters["noisetimestamp"];
  // mTimestamp = a;//std::stol(a,nullptr,10);
  // long int ts = mTimestamp ? mTimestamp : o2::ccdb::getCurrentTimestamp();
  // TODO come back and all for different time stamps
  long int ts = o2::ccdb::getCurrentTimestamp();
  ILOG(Info, Support) << "Getting noisemap from ccdb - timestamp: " << ts << ENDM;
  auto& mgr = o2::ccdb::BasicCCDBManager::instance();
  mgr.setURL("http://alice-ccdb.cern.ch");
  mgr.setTimestamp(ts);
  mNoiseMap = mgr.get<o2::trd::NoiseStatusMCM>("/TRD/Calib/NoiseMapMCM");
  if (mNoiseMap == nullptr) {
    ILOG(Info, Support) << "mNoiseMap is null, no noisy mcm reduction" << ENDM;
  }

  // HalfChamberMask(mgr);
}

void TrackletsTask::buildHistograms()
{
  for (Int_t sm = 0; sm < o2::trd::constants::NSECTOR; ++sm) {
    std::string label = fmt::format("TrackletHCMCM_{0}", sm);
    std::string title = fmt::format("MCM in Tracklets data stream for sector {0}", sm);
    moHCMCM[sm].reset(new TH2F(label.c_str(), title.c_str(), 76, -0.5, 75.5, 8 * 5, -0.5, 8 * 5 - 0.5));
    moHCMCM[sm]->GetYaxis()->SetTitle("ROB in stack");
    moHCMCM[sm]->GetXaxis()->SetTitle("mcm in rob in layer");
    getObjectsManager()->startPublishing(moHCMCM[sm].get());
    getObjectsManager()->setDefaultDrawOptions(moHCMCM[sm]->GetName(), "COLZ");
    drawLinesMCM(moHCMCM[sm].get());
  }
  mTrackletSlope.reset(new TH1F("trackletslope", "uncalibrated Slope of tracklets", 1024, -6.0, 6.0)); // slope is 8 bits in the tracklet
  getObjectsManager()->startPublishing(mTrackletSlope.get());
  mTrackletSlopeRaw.reset(new TH1F("trackletsloperaw", "Raw Slope of tracklets", 256, 0, 256)); // slope is 8 bits in the tracklet
  getObjectsManager()->startPublishing(mTrackletSlopeRaw.get());
  mTrackletHCID.reset(new TH1F("tracklethcid", "Tracklet distribution over Halfchambers", 1080, 0, 1080));
  getObjectsManager()->startPublishing(mTrackletHCID.get());
  mTrackletPosition.reset(new TH1F("trackletpos", "Uncalibrated Position of Tracklets", 1400, -70, 70));
  getObjectsManager()->startPublishing(mTrackletPosition.get());
  mTrackletPositionRaw.reset(new TH1F("trackletposraw", "Raw Position of Tracklets", 2048, 0, 2048));
  getObjectsManager()->startPublishing(mTrackletPositionRaw.get());
  mTrackletsPerEvent.reset(new TH1F("trackletsperevent", "Number of Tracklets per event", 2500, 0, 250000));
  getObjectsManager()->startPublishing(mTrackletsPerEvent.get());

  for (Int_t sm = 0; sm < o2::trd::constants::NSECTOR; ++sm) {
    std::string label = fmt::format("TrackletHCMCMnoise_{0}", sm);
    std::string title = fmt::format("MCM in Tracklets data stream for sector {0} noise in", sm);
    moHCMCMn[sm].reset(new TH2F(label.c_str(), title.c_str(), 76, -0.5, 75.5, 8 * 5, -0.5, 8 * 5 - 0.5));
    moHCMCMn[sm]->GetYaxis()->SetTitle("ROB in stack");
    moHCMCMn[sm]->GetXaxis()->SetTitle("mcm in rob in layer");
    getObjectsManager()->startPublishing(moHCMCMn[sm].get());
    getObjectsManager()->setDefaultDrawOptions(moHCMCMn[sm]->GetName(), "COLZ");
    drawLinesMCM(moHCMCMn[sm].get());
  }
  mTrackletSlopen.reset(new TH1F("trackletslopenoise", "uncalibrated Slope of tracklets noise in", 1024, -6.0, 6.0)); // slope is 8 bits in the tracklet
  getObjectsManager()->startPublishing(mTrackletSlopen.get());
  mTrackletSlopeRawn.reset(new TH1F("trackletsloperawnoise", "Raw Slope of tracklets noise in", 256, 0, 256)); // slope is 8 bits in the tracklet
  getObjectsManager()->startPublishing(mTrackletSlopeRawn.get());
  mTrackletHCIDn.reset(new TH1F("tracklethcidnoise", "Tracklet distribution over Halfchambers noise in", 1080, 0, 1080));
  getObjectsManager()->startPublishing(mTrackletHCIDn.get());
  mTrackletPositionn.reset(new TH1F("trackletposnoise", "Uncalibrated Position of Tracklets noise in", 1400, -70, 70));
  getObjectsManager()->startPublishing(mTrackletPositionn.get());
  mTrackletPositionRawn.reset(new TH1F("trackletposrawnoise", "Raw Position of Tracklets noise in", 2048, 0, 2048));
  getObjectsManager()->startPublishing(mTrackletPositionRawn.get());
  mTrackletsPerEventn.reset(new TH1F("trackletspereventn", "Number of Tracklets per event noise in", 2500, 0, 250000));
  getObjectsManager()->startPublishing(mTrackletsPerEventn.get());

  // TCanvas* c[6];

  for (Int_t iLayer = 0; iLayer < 6; ++iLayer) {
    std::string labelc = fmt::format("TCperMCMinLayer{0}", iLayer);
    c[iLayer].reset(new TCanvas(labelc.c_str()));

    c[iLayer]->cd();

    std::string label = fmt::format("TCperMCMinLayer{0}", iLayer);
    std::string title = fmt::format("Tracklet count per MCM in layer {0}", iLayer);
    hLayers[iLayer].reset(new TH2F(label.c_str(), title.c_str(), 76, -0.5, 75.5, 144, -0.5, 143.5));
    hLayers[iLayer]->GetYaxis()->SetTitle("stack");
    hLayers[iLayer]->GetXaxis()->SetTitle("sector");

    auto xax = hLayers[iLayer].get()->GetXaxis();
    xax->SetBinLabel(8, "0");
    xax->SetBinLabel(24, "1");
    xax->SetBinLabel(38, "2");
    xax->SetBinLabel(52, "3");
    xax->SetBinLabel(68, "4");
    xax->SetTicks("-");
    xax->SetTickSize(0.01);
    xax->SetLabelSize(0.045);
    // xax->SetLabelOffset(0.01);
    //  xax->SetTitleOffset(0.8);
    auto yax = hLayers[iLayer].get()->GetYaxis();
    for (int iSec = 0; iSec < 18; ++iSec) {
      auto lbl = std::to_string(iSec);
      yax->SetBinLabel(iSec * 8 + 4, lbl.c_str());
    }
    yax->SetTicks("-");
    yax->SetTickSize(0.01);
    yax->SetLabelSize(0.045);
    yax->SetLabelOffset(0.01);
    yax->SetTitleOffset(1.4);
    hLayers[iLayer].get()->SetStats(0);
    hLayers[iLayer].get()->Draw("COLZ ");

    hMask[iLayer].reset(new TH2F(Form("layer%i_mask", iLayer), "", 76, -0.5, 75.5, 144, -0.5, 143.5));
    hMask[iLayer].get()->SetMarkerColor(kRed);
    hMask[iLayer].get()->SetMarkerSize(0.9);
    hMask[iLayer].get()->Draw("text same");

    getObjectsManager()->startPublishing(c[iLayer].get());
    drawLinesMCM(hLayers[iLayer].get(), 144, 3);
  }
}

void TrackletsTask::initialize(o2::framework::InitContext& /*ctx*/)
{
  ILOG(Info, Support) << "initialize TrackletsTask" << ENDM;

  buildHistograms();
  retrieveCCDBSettings();
}

void TrackletsTask::startOfActivity(Activity& activity)
{
  ILOG(Info, Support) << "startOfActivity " << activity.mId << ENDM;
  for (Int_t sm = 0; sm < o2::trd::constants::NSECTOR; ++sm) {
    moHCMCM[sm]->Reset();
  }
}

void TrackletsTask::startOfCycle()
{
  ILOG(Info, Support) << "startOfCycle" << ENDM;
}

void TrackletsTask::monitorData(o2::framework::ProcessingContext& ctx)
{
  for (auto&& input : ctx.inputs()) {
    if (input.header != nullptr && input.payload != nullptr) {

      auto digits = ctx.inputs().get<gsl::span<o2::trd::Digit>>("digits");
      auto tracklets = ctx.inputs().get<gsl::span<o2::trd::Tracklet64>>("tracklets");
      auto triggerrecords = ctx.inputs().get<gsl::span<o2::trd::TriggerRecord>>("triggers");

      for (auto& trigger : triggerrecords) {
        if (trigger.getNumberOfTracklets() == 0)
          continue; // bail if we have no digits in this trigger
        // now sort digits to det,row,pad
        mTrackletsPerEvent->Fill(trigger.getNumberOfTracklets());
        for (int currenttracklet = trigger.getFirstTracklet(); currenttracklet < trigger.getFirstTracklet() + trigger.getNumberOfTracklets() - 1; ++currenttracklet) {
          int detector = tracklets[currenttracklet].getDetector();
          int sm = detector / 30;
          int detLoc = detector % 30;
          int layer = detector % 6;
          int istack = detLoc / 6;
          int iChamber = sm * 30 + istack * o2::trd::constants::NLAYER + layer;
          int stackoffset = istack * o2::trd::constants::NSTACK * o2::trd::constants::NROBC1;
          if (istack >= 2) {
            stackoffset -= 2; // only 12in stack 2
          }

          // 8 rob x 16 mcm each per chamber
          //  5 stack(y), 6 layers(x)
          //  y=stack_rob, x=layer_mcm
          int x = o2::trd::constants::NMCMROB * layer + tracklets[currenttracklet].getMCM();
          int y = o2::trd::constants::NROBC1 * istack + tracklets[currenttracklet].getROB();
          if (mNoiseMap != nullptr && mNoiseMap->isTrackletFromNoisyMCM(tracklets[currenttracklet])) {
            moHCMCMn[sm]->Fill(x, y);
            mTrackletSlopen->Fill(tracklets[currenttracklet].getUncalibratedDy());
            mTrackletSlopeRawn->Fill(tracklets[currenttracklet].getSlope());
            mTrackletPositionn->Fill(tracklets[currenttracklet].getUncalibratedY());
            mTrackletPositionRawn->Fill(tracklets[currenttracklet].getPosition());
            mTrackletHCIDn->Fill(tracklets[currenttracklet].getHCID());
          } else {
            moHCMCM[sm]->Fill(x, y);
            mTrackletSlope->Fill(tracklets[currenttracklet].getUncalibratedDy());
            mTrackletSlopeRaw->Fill(tracklets[currenttracklet].getSlope());
            mTrackletPosition->Fill(tracklets[currenttracklet].getUncalibratedY());
            mTrackletPositionRaw->Fill(tracklets[currenttracklet].getPosition());
            mTrackletHCID->Fill(tracklets[currenttracklet].getHCID());
          }

          int rowGlb = istack < 3 ? tracklets[currenttracklet].getPadRow() + istack * 16 : tracklets[currenttracklet].getPadRow() + 44 + (istack - 3) * 16; // pad row within whole sector
          int colGlb = tracklets[currenttracklet].getColumn() + sm * 8 + sm * 4;
          hLayers[layer]->Fill(rowGlb, colGlb);
        }
      }
    }
  }

  fillMaskHisto();
}
void TrackletsTask::fillMaskHisto()
{
  auto chamberStatus = halfChamberMask();
  for (int iSec = 0; iSec < 18; ++iSec) {
    for (int iStack = 0; iStack < 5; ++iStack) {
      int rowMax = (iStack == 2) ? 12 : 16;
      for (int iLayer = 0; iLayer < 6; ++iLayer) {
        for (int iCol = 0; iCol < 8; ++iCol) {
          int side = (iCol < 4) ? 0 : 1;
          int det = iSec * 30 + iStack * 6 + iLayer;
          int hcid = (side == 0) ? det * 2 : det * 2 + 1;
          for (int iRow = 0; iRow < rowMax; ++iRow) {
            int rowGlb = iStack < 3 ? iRow + iStack * 16 : iRow + 44 + (iStack - 3) * 16; // pad row within whole sector
            int colGlb = iCol + iSec * 8;
            // bin number 0 is underflow
            rowGlb += 1;
            colGlb += 1;
            if (chamberStatus.isMasked(hcid)) {
              hMask[iLayer]->SetBinContent(rowGlb, colGlb, 1);
            }
          }
        }
      }
    }
  }
}

o2::trd::HalfChamberStatusQC TrackletsTask::halfChamberMask()
{

  o2::trd::HalfChamberStatusQC status;

  std::vector<std::tuple<int, int, int>> mC{ { 0, 0, 1 }, { 0, 0, 3 }, { 0, 1, 0 }, { 0, 4, 2 }, { 1, 2, 1 }, { 1, 3, 1 }, { 1, 4, 1 }, { 3, 0, 2 }, { 3, 0, 3 }, { 3, 1, 2 }, { 3, 2, 2 }, { 3, 2, 4 }, { 3, 2, 5 }, { 3, 3, 2 }, { 3, 3, 3 }, { 3, 3, 5 }, { 3, 4, 2 }, { 3, 4, 5 }, { 4, 1, 2 }, { 4, 1, 5 }, { 4, 2, 0 }, { 5, 1, 0 }, { 5, 1, 5 }, { 5, 2, 5 }, { 5, 2, 3 }, { 5, 3, 1 }, { 5, 4, 2 }, { 7, 0, 5 }, { 7, 4, 2 }, { 7, 4, 4 }, { 9, 2, 5 }, { 9, 3, 2 }, { 10, 4, 3 }, { 11, 0, 5 }, { 11, 2, 0 }, { 11, 2, 3 }, { 11, 3, 0 }, { 11, 4, 0 }, { 11, 4, 5 }, { 12, 4, 5 }, { 13, 0, 0 }, { 13, 1, 2 }, { 13, 2, 0 }, { 13, 2, 1 }, { 13, 2, 2 }, { 13, 2, 3 }, { 13, 2, 4 }, { 13, 2, 5 }, { 13, 4, 5 }, { 14, 2, 0 }, { 14, 2, 0 }, { 14, 2, 1 }, { 14, 2, 2 }, { 14, 2, 3 }, { 14, 2, 4 }, { 14, 2, 5 }, { 14, 3, 1 }, { 15, 0, 2 }, { 15, 0, 5 }, { 15, 1, 0 }, { 15, 1, 1 }, { 15, 2, 0 }, { 15, 2, 1 }, { 15, 2, 2 }, { 15, 2, 3 }, { 15, 2, 4 }, { 15, 3, 2 }, { 15, 3, 5 }, { 15, 4, 0 }, { 15, 4, 2 }, { 16, 3, 2 }, { 16, 3, 4 }, { 17, 0, 0 }, { 17, 0, 5 }, { 17, 1, 0 }, { 17, 2, 1 }, { 17, 2, 4 }, { 17, 3, 0 }, { 17, 3, 1 }, { 17, 4, 2 }, { 17, 4, 4 } };

  std::vector<std::tuple<int, int, int>> mHCA{ { 0, 2, 4 }, { 1, 0, 4 }, { 4, 0, 0 }, { 4, 0, 5 }, { 5, 3, 4 }, { 6, 0, 2 }, { 6, 0, 3 }, { 8, 0, 1 }, { 9, 2, 2 }, { 10, 1, 1 }, { 10, 3, 1 }, { 15, 1, 3 }, { 16, 0, 2 } };
  std::vector<std::tuple<int, int, int>> mHCB{ { 1, 1, 3 }, { 1, 2, 5 }, { 7, 4, 3 }, { 7, 4, 5 }, { 8, 0, 3 }, { 8, 2, 5 }, { 10, 4, 0 }, { 11, 1, 2 }, { 12, 2, 0 }, { 13, 1, 0 }, { 13, 1, 5 }, { 15, 0, 3 }, { 15, 3, 0 }, { 16, 0, 0 }, { 16, 1, 2 }, { 16, 4, 4 }, { 17, 0, 4 } };

  for (int mc = 0; mc < mC.size(); ++mc) {
    status.maskChamber(std::get<0>(mC[mc]), std::get<1>(mC[mc]), std::get<2>(mC[mc]));
  }
  for (int mc = 0; mc < mHCA.size(); ++mc) {
    status.maskHalfChamberA(std::get<0>(mHCA[mc]), std::get<1>(mHCA[mc]), std::get<2>(mHCA[mc]));
  }
  for (int mc = 0; mc < mHCB.size(); ++mc) {
    status.maskHalfChamberB(std::get<0>(mHCB[mc]), std::get<1>(mHCB[mc]), std::get<2>(mHCB[mc]));
  }
  return status;
}
void TrackletsTask::endOfCycle()
{
  ILOG(Info, Support) << "endOfCycle" << ENDM;
  // scale 2d mHCMCM plots so they all have the same max height.
  int max = 0;
  for (auto& hist : moHCMCM) {
    if (hist->GetMaximum() > max) {
      max = hist->GetMaximum();
    }
  }
  for (auto& hist : moHCMCM) {
    hist->SetMaximum(max);
  }
}

void TrackletsTask::endOfActivity(Activity& /*activity*/)
{
  ILOG(Info, Support) << "endOfActivity" << ENDM;
}

void TrackletsTask::reset()
{
  // clean all the monitor objects here

  ILOG(Info, Support) << "Resetting the histogram" << ENDM;
  mTrackletPosition.get()->Reset();
  mTrackletPositionRaw.get()->Reset();
  mTrackletSlope.get()->Reset();
  mTrackletSlopeRaw.get()->Reset();
  mTrackletHCID.get()->Reset();
  mTrackletsPerEvent.get()->Reset();
}

} // namespace o2::quality_control_modules::trd
