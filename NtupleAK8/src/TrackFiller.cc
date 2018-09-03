/*
 * TrackFiller.cc
 *
 *  Created on: May 24, 2017
 *      Author: hqu
 */

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DeepNTuples/NtupleAK8/interface/SVFiller.h"

#include <unordered_map>
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DeepNTuples/NtupleAK8/interface/TrackFiller.h"

namespace deepntuples {

void TrackFiller::readConfig(const edm::ParameterSet& iConfig, edm::ConsumesCollector&& cc) {
  vtxToken_ = cc.consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  svToken_ = cc.consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("SVs"));
}

void TrackFiller::readEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  iEvent.getByToken(vtxToken_, vertices);
  iEvent.getByToken(svToken_, SVs);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder_);
}

void TrackFiller::book() {

  data.add<int>("n_tracks", 0);
  data.add<float>("ntracks", 0);

  // adding displaced jet info
  data.add<float>("alpha_loose", -1);
  data.add<float>("alphaStar_loose", -1);
  data.add<float>("alpha_ratio_loose", -1);

  data.add<float>("alpha_tight", -1);
  data.add<float>("alphaStar_tight", -1);
  data.add<float>("alpha_ratio_tight", -1);

  data.add<float>("alpha_used", -1);
  data.add<float>("alphaStar_used", -1);
  data.add<float>("alpha_ratio_used", -1);

  data.add<float>("alpha_NoVert", -1);
  data.add<float>("alphaStar_NoVert", -1);
  data.add<float>("alpha_ratio_NoVert", -1);

  data.add<float>("JetDeltaR", -1);
  data.add<float>("Jet2DIPSigMedian", -1);
  data.add<float>("Jet3DIPSigMedian", -1);
  data.add<float>("Jet2DIPSig", -1);
  data.add<float>("Jet3DIPSig", -1);

  data.add<float>("SV_alpha_Max", -1);
  data.add<float>("SV_ThetaMedian", -1);
  //data.add<float>("alpha", -1);
  //data.add<float>("alphaMax", -1);
  //data.add<float>("alphaSV", -1);
  //data.add<float>("alphaSVMax", -1);

  // basic kinematics
  data.addMulti<float>("track_ptrel");
  data.addMulti<float>("track_erel");
  data.addMulti<float>("track_phirel");
  data.addMulti<float>("track_etarel");
  data.addMulti<float>("track_deltaR");
  data.addMulti<float>("track_puppiw");
  data.addMulti<float>("track_pt");
  data.addMulti<float>("track_mass");

  data.addMulti<float>("track_drminsv");
  data.addMulti<float>("track_drsubjet1");
  data.addMulti<float>("track_drsubjet2");

  data.addMulti<float>("track_charge");
  data.addMulti<float>("track_isMu");
  data.addMulti<float>("track_isEl");
  data.addMulti<float>("track_isChargedHad");

  // for charged
  data.addMulti<float>("track_VTX_ass");
  data.addMulti<float>("track_fromPV");
  data.addMulti<float>("track_lostInnerHits");

  // impact parameters
  data.addMulti<float>("track_dz");
  data.addMulti<float>("track_dzsig");
  data.addMulti<float>("track_dxy");
  data.addMulti<float>("track_dxysig");

  // track quality
  data.addMulti<float>("track_normchi2");
  data.addMulti<float>("track_quality");

  // track covariance
  data.addMulti<float>("track_dptdpt");
  data.addMulti<float>("track_detadeta");
  data.addMulti<float>("track_dphidphi");
  data.addMulti<float>("track_dxydxy");
  data.addMulti<float>("track_dzdz");
  data.addMulti<float>("track_dxydz");
  data.addMulti<float>("track_dphidxy");
  data.addMulti<float>("track_dlambdadz");

  // track btag info
  data.addMulti<float>("trackBTag_Momentum");
  data.addMulti<float>("trackBTag_Eta");
  data.addMulti<float>("trackBTag_EtaRel");
  data.addMulti<float>("trackBTag_PtRel");
  data.addMulti<float>("trackBTag_PPar");
  data.addMulti<float>("trackBTag_DeltaR");
  data.addMulti<float>("trackBTag_PtRatio");
  data.addMulti<float>("trackBTag_PParRatio");
  data.addMulti<float>("trackBTag_Sip2dVal");
  data.addMulti<float>("trackBTag_Sip2dSig");
  data.addMulti<float>("trackBTag_Sip3dVal");
  data.addMulti<float>("trackBTag_Sip3dSig");
  data.addMulti<float>("trackBTag_JetDistVal");
//  data.addMulti<float>("trackBTag_JetDistSig"); // always gives 0


}

float CalcMedian(std::vector<float> vec)
{
  size_t size = vec.size();
  if (size == 0){
    return 0;  // Undefined, really.
  }else{
    //std::sort(vec.begin(), vec.end());
    if (size % 2 == 0){
      return (vec[size / 2 - 1] + vec[size / 2]) / 2;
    }else {
      return vec[size / 2];
    }
  }
}

Measurement1D TrackFiller::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

bool TrackFiller::fill(const pat::Jet& jet, size_t jetidx, const JetHelper& jet_helper) {

  std::vector<const pat::PackedCandidate*> chargedPFCands;
  std::unordered_map<const pat::PackedCandidate*, TrackInfoBuilder> trackInfoMap;
  for (const auto * pfcand : jet_helper.getJetConstituents()){
    if (pfcand->charge() != 0) {
      chargedPFCands.push_back(pfcand);
      trackInfoMap[pfcand];
      trackInfoMap[pfcand].buildTrackInfo(builder_, *pfcand, jet, vertices->at(0));
    }
  }
  
  // sort by Sip2d significance
  std::sort(chargedPFCands.begin(), chargedPFCands.end(), [&](const pat::PackedCandidate *p1, const pat::PackedCandidate *p2){
    return trackInfoMap.at(p1).getTrackSip2dSig() > trackInfoMap.at(p2).getTrackSip2dSig();
  });
  
  // calculating the alpha for pv and sv
  //float sumJetPt  = 0;
  //float alpha     = 0;
  
  //const auto &lead_pv = vertices->at(0);
  ////bool in_lead_pv = false; 
  //// const auto &lead_sv = SVs->at(0);
  //// std::cout << " Jet Radius : " << jetR_ << std::endl;
  //// these are the sums for this specific vertex
  //for (const auto *tkr : chargedPFCands){
  //  if((tkr->pt() < 1.0)) continue; //apply a cut on the track pt 
  //  //float lead_pv_dz = std::abs(tkr->trackRef()->dz(lead_pv->position()));
  //  int vtx_i = 0 ;
  //  for (const auto &iv : *vertices){
  //    if( iv.isFake() || iv.ndof() < 4 ) continue;
  //    bool is_lead_pv  = (iv.position() - lead_pv.position()).r() < 0.02;
  //    std::cout << " -- " << vtx_i << " : "<< is_lead_pv << " -- trak dz(pv)" << tkr->dz(iv.position()) << std::endl;
  //    vtx_i ++ ;
  //  }
 //}

   
  //***********************************START alpha variables for vertices***********************************

  float alpha_den = 0;
  float alpha_num_loose = 0;
  float alphaStar_num_loose = 0;
  float alpha_num_tight = 0;
  float alphaStar_num_tight = 0;
  float alpha_num_used = 0;
  float alphaStar_num_used = 0;
  float alpha_num_NoVert = 0;
  float alphaStar_num_NoVert = 0;
  float alphaStar_num_Max_loose = 0;
  float alphaStar_num_Max_tight = 0;
  float alphaStar_num_Max_used = 0;
  float alphaStar_num_Max_NoVert = 0;
  float alpha_loose = 0;
  float alpha_tight = 0;
  float alpha_used = 0;
  float alpha_NoVert = 0;
  int sizeV = vertices->size();

  for (int vID =0 ; vID < sizeV ; vID++){
      alpha_num_loose = 0;
      alphaStar_num_loose = 0;
      alpha_num_tight = 0;
      alphaStar_num_tight = 0;
      alpha_num_used = 0;
      alphaStar_num_used = 0;
      alpha_num_NoVert = 0;
      alphaStar_num_NoVert = 0;
      alpha_den = 0;
      for (const auto *cpf : chargedPFCands){
          if (cpf->pt() < 1) continue ;
          alpha_den += cpf->pt();
          if (cpf->fromPV(vID) >= 1){
             if (vID == 0) alpha_num_loose += cpf->pt();
             if (vID != 0) alphaStar_num_loose += cpf->pt();
          }if (cpf->fromPV(vID) >= 2){
             if (vID == 0) alpha_num_tight += cpf->pt();
             if (vID != 0) alphaStar_num_tight += cpf->pt();
          }if (cpf->fromPV(vID) >= 3){
             if (vID == 0) alpha_num_used += cpf->pt();
             if (vID != 0) alphaStar_num_used += cpf->pt();
          }if (cpf->fromPV(vID) == 0){
             if (vID == 0) alpha_num_NoVert += cpf->pt();
             if (vID != 0) alphaStar_num_NoVert += cpf->pt();
          }
     }
     if (alphaStar_num_loose >= alphaStar_num_Max_loose) alphaStar_num_Max_loose = alphaStar_num_loose;
     if (alphaStar_num_tight >= alphaStar_num_Max_tight) alphaStar_num_Max_tight = alphaStar_num_tight;
     if (alphaStar_num_used >= alphaStar_num_Max_used) alphaStar_num_Max_used = alphaStar_num_used;
     if (alphaStar_num_NoVert >= alphaStar_num_Max_NoVert) alphaStar_num_Max_NoVert = alphaStar_num_NoVert;  
     if (alpha_den == 0) alpha_den = -1;
     if (vID == 0) {
        alpha_loose = alpha_num_loose / alpha_den;
        alpha_tight = alpha_num_tight / alpha_den;
        alpha_used = alpha_num_used / alpha_den;
        alpha_NoVert = alpha_num_NoVert / alpha_den;
     }
  }
  float alphaStar_loose = alphaStar_num_Max_loose / alpha_den;
  float alphaStar_tight = alphaStar_num_Max_tight / alpha_den;
  float alphaStar_used = alphaStar_num_Max_used / alpha_den;
  float alphaStar_NoVert = alphaStar_num_Max_NoVert / alpha_den;

  float alpha_ratio_loose = alphaStar_loose / alpha_loose;
  if (alpha_loose == 0 ) alpha_ratio_loose = 0;
  float alpha_ratio_tight = alphaStar_tight / alpha_tight;
  if (alpha_tight == 0 ) alpha_ratio_tight = 0;
  float alpha_ratio_used = alphaStar_used / alpha_used;
  if (alpha_used == 0 ) alpha_ratio_used = 0;
  float alpha_ratio_NoVert = alphaStar_NoVert / alpha_NoVert;
  if (alpha_NoVert == 0 ) alpha_ratio_NoVert = 0;

  //***********************************END alpha variables for vertices***********************************
  float den = 0;
  float JetDeltaR_num = 0;
  float IPSig2D_num = 0;
  float IPSig3D_num = 0;
  std::vector<float> IPSig2D_List ;
  std::vector<float> IPSig3D_List ; 
  for (const auto *cpf : chargedPFCands){
      if (cpf->pt() < 1) continue ;
      den += cpf->pt();
      const auto &trkinfo = trackInfoMap.at(cpf);
      JetDeltaR_num += trkinfo.getTrackDeltaR()  * cpf->pt();
      //2dIPSig
      //if ( trkinfo.getTrackSip2dVal() / trkinfo.getTrackSip2dSig() != trkinfo.getTrackSip2dVal() / trkinfo.getTrackSip2dSig()) || trkinfo.getTrackSip2dVal() / trkinfo.getTrackSip2dSig() == abs(inf) ) continue;
      if ( isnan(trkinfo.getTrackSip2dVal() / trkinfo.getTrackSip2dSig()) == 1 || isinf(trkinfo.getTrackSip2dVal() / trkinfo.getTrackSip2dSig()) == 1) continue; // {
      IPSig2D_num += trkinfo.getTrackSip2dVal() / trkinfo.getTrackSip2dSig();
      IPSig2D_List.push_back(IPSig2D_num);
      //}
      //3DIPSig
      if ( isnan(trkinfo.getTrackSip3dVal() / trkinfo.getTrackSip3dSig()) == 0 || isinf(trkinfo.getTrackSip3dVal() / trkinfo.getTrackSip3dSig()) == 0) {
      IPSig3D_num += trkinfo.getTrackSip3dVal() / trkinfo.getTrackSip3dSig();
      IPSig3D_List.push_back(IPSig3D_num);
      }
  }
  //Median Variables
  std::sort ( IPSig2D_List.begin(), IPSig2D_List.end() );
  std::sort ( IPSig3D_List.begin(), IPSig3D_List.end() );
  float Jet2DIPSigMedian = CalcMedian(IPSig2D_List);
  float Jet3DIPSigMedian = CalcMedian(IPSig3D_List);
  //Mean Variables
  float JetDeltaR = JetDeltaR_num / den;
  int denom = chargedPFCands.size();
  float Jet2DIPSig = IPSig2D_num / denom;
  float Jet3DIPSig = IPSig3D_num / denom;


  //Now Secondary Vertex stuff
  float SV_alpha_num = 0; 
  float SV_alpha_Max = 0;
  int counting = 0;

  const auto &lead_pv = vertices->at(0); 
  float dot = 0;
  float lenSq1 = 0;
  float lenSq2 = 0;
  std::vector<float> angle;

  for (const auto &sv : *SVs){
      SV_alpha_num = 0;
      counting += 1;   
      size_t SV_Size = sv.numberOfDaughters();
      for (size_t Did = 0 ; Did < SV_Size ; Did++){
          float PT = sv.daughter(Did)->pt();
          //std::cout << "Index:  " << Did << "   PT:   "  << PT <<  "X:   " << sv.position().x() << "   track X:    " << sv.daughter(Did)->vx() <<  "   Y:   " << sv.position().y() <<  "   track Y:    " << sv.daughter(Did)->vy() <<  "   Z:   "  << sv.position().z() << "   track Z:    " << sv.daughter(Did)->vz() <<  std::endl;
          double delr = reco::deltaR(sv.daughter(Did)->p4(), sv);
          if (PT < 1) continue;
          if (delr > 0.8) continue;
          SV_alpha_num += PT;

         //Theta Values now
         dot = ( ( sv.position().x() - lead_pv.x() ) * sv.daughter(Did)->p4().px() ) + ( ( sv.position().y() - lead_pv.y() ) * sv.daughter(Did)->p4().py() ) + ( ( sv.position().z() - lead_pv.z() ) * sv.daughter(Did)->p4().pz() ) ;
         lenSq1 = ( sv.position().x() - lead_pv.x() ) * ( sv.position().x() - lead_pv.x() ) + ( sv.position().y() - lead_pv.y() ) * ( sv.position().y() - lead_pv.y() ) + ( sv.position().z() - lead_pv.z() ) * ( sv.position().z() - lead_pv.z() ) ;  
         lenSq2 = sv.daughter(Did)->p4().px() * sv.daughter(Did)->p4().px() + sv.daughter(Did)->p4().py() * sv.daughter(Did)->p4().py() + sv.daughter(Did)->p4().pz() * sv.daughter(Did)->p4().pz() ;
         float Cosangle = acos(dot/sqrt(lenSq1*lenSq2)) ; 
         angle.push_back(Cosangle);
      }
      if (SV_alpha_num > SV_alpha_Max) SV_alpha_Max = SV_alpha_num;
  }

  for  (unsigned int j =0 ; j < angle.size() ; j++ ) {
       std::cout << angle[j] << "   " ;
  }
  std::cout << std::endl; 
  float MedianTheta = 0;
  std::sort( angle.begin(), angle.end() );
  MedianTheta = CalcMedian(angle);
  std::cout << "the Median Now:  " << MedianTheta << std::endl;  

  SV_alpha_Max = SV_alpha_Max / den;
 
  data.fill<float>("alpha_loose", alpha_loose);
  data.fill<float>("alphaStar_loose", alphaStar_loose);
  data.fill<float>("alpha_ratio_loose", alpha_ratio_loose);

  data.fill<float>("alpha_tight", alpha_tight);
  data.fill<float>("alphaStar_tight", alphaStar_tight);
  data.fill<float>("alpha_ratio_tight", alpha_ratio_tight);

  data.fill<float>("alpha_used", alpha_used);
  data.fill<float>("alphaStar_used", alphaStar_used);
  data.fill<float>("alpha_ratio_used", alpha_ratio_used);

  data.fill<float>("alpha_NoVert", alpha_NoVert);
  data.fill<float>("alphaStar_NoVert", alphaStar_NoVert);
  data.fill<float>("alpha_ratio_NoVert", alpha_ratio_NoVert); 
  
  data.fill<float>("JetDeltaR", JetDeltaR);
  data.fill<float>("Jet2DIPSigMedian", Jet2DIPSigMedian);
  data.fill<float>("Jet3DIPSigMedian", Jet3DIPSigMedian);
  data.fill<float>("Jet2DIPSig", Jet2DIPSig);
  data.fill<float>("Jet3DIPSig", Jet3DIPSig);

  data.fill<float>("SV_alpha_Max", SV_alpha_Max);
  data.fill<float>("SV_ThetaMedian", MedianTheta);

  data.fill<int>("n_tracks", chargedPFCands.size());
  data.fill<float>("ntracks", chargedPFCands.size());

  float etasign = jet.eta()>0 ? 1 : -1;

  for (const auto *cpf : chargedPFCands){

    // basic kinematics, valid for both charged and neutral
    data.fillMulti<float>("track_ptrel", cpf->pt()/jet.pt());
    data.fillMulti<float>("track_erel", cpf->energy()/jet.energy());
    data.fillMulti<float>("track_phirel", reco::deltaPhi(*cpf, jet));
    data.fillMulti<float>("track_etarel", etasign * (cpf->eta() - jet.eta()));
    data.fillMulti<float>("track_deltaR", reco::deltaR(*cpf, jet));
    data.fillMulti<float>("track_puppiw", cpf->puppiWeight());
    data.fillMulti<float>("track_pt", cpf->pt());
    data.fillMulti<float>("track_mass", cpf->mass());

    double minDR = 999;
    for (const auto &sv : *SVs){
      double dr = reco::deltaR(*cpf, sv);
      if (dr < minDR) minDR = dr;
    }
    data.fillMulti<float>("track_drminsv", minDR==999 ? -1 : minDR);

    const auto& subjets = jet_helper.getSubJets();
    data.fillMulti<float>("track_drsubjet1", subjets.size()>0 ? reco::deltaR(*cpf, *subjets.at(0)) : -1);
    data.fillMulti<float>("track_drsubjet2", subjets.size()>1 ? reco::deltaR(*cpf, *subjets.at(1)) : -1);

    data.fillMulti<float>("track_charge", cpf->charge());
    data.fillMulti<float>("track_isEl", std::abs(cpf->pdgId())==11);
    data.fillMulti<float>("track_isMu", std::abs(cpf->pdgId())==13);
    data.fillMulti<float>("track_isChargedHad", std::abs(cpf->pdgId())==211);

    // for charged
    data.fillMulti<float>("track_VTX_ass", cpf->pvAssociationQuality());
    data.fillMulti<float>("track_fromPV", cpf->fromPV());
    data.fillMulti<float>("track_lostInnerHits", cpf->lostInnerHits());

    // impact parameters
    data.fillMulti<float>("track_dz", catchInfs(cpf->dz()));
    data.fillMulti<float>("track_dzsig", catchInfs(cpf->dz()/cpf->dzError()));
    data.fillMulti<float>("track_dxy", catchInfs(cpf->dxy()));
    data.fillMulti<float>("track_dxysig", catchInfs(cpf->dxy()/cpf->dxyError()));

    const auto &trk = cpf->pseudoTrack();
    data.fillMulti<float>("track_normchi2", catchInfs(trk.normalizedChi2()));
    data.fillMulti<float>("track_quality", trk.qualityMask());

    // track covariance
    auto cov = [&](unsigned i, unsigned j) {
      return catchInfs(trk.covariance(i, j));
    };
    data.fillMulti<float>("track_dptdpt", cov(0,0));
    data.fillMulti<float>("track_detadeta", cov(1,1));
    data.fillMulti<float>("track_dphidphi", cov(2,2));
    data.fillMulti<float>("track_dxydxy", cov(3,3));
    data.fillMulti<float>("track_dzdz", cov(4,4));
    data.fillMulti<float>("track_dxydz", cov(3,4));
    data.fillMulti<float>("track_dphidxy", cov(2,3));
    data.fillMulti<float>("track_dlambdadz", cov(1,4));

    const auto &trkinfo = trackInfoMap.at(cpf);
    data.fillMulti<float>("trackBTag_Momentum", trkinfo.getTrackMomentum());
    data.fillMulti<float>("trackBTag_Eta", trkinfo.getTrackEta());
    data.fillMulti<float>("trackBTag_EtaRel", trkinfo.getTrackEtaRel());
    data.fillMulti<float>("trackBTag_PtRel", trkinfo.getTrackPtRel());
    data.fillMulti<float>("trackBTag_PPar", trkinfo.getTrackPPar());
    data.fillMulti<float>("trackBTag_DeltaR", trkinfo.getTrackDeltaR());
    data.fillMulti<float>("trackBTag_PtRatio", trkinfo.getTrackPtRatio());
    data.fillMulti<float>("trackBTag_PParRatio", trkinfo.getTrackPParRatio());
    data.fillMulti<float>("trackBTag_Sip2dVal", trkinfo.getTrackSip2dVal());
    data.fillMulti<float>("trackBTag_Sip2dSig", trkinfo.getTrackSip2dSig());
    data.fillMulti<float>("trackBTag_Sip3dVal", trkinfo.getTrackSip3dVal());
    data.fillMulti<float>("trackBTag_Sip3dSig", trkinfo.getTrackSip3dSig());
    data.fillMulti<float>("trackBTag_JetDistVal", trkinfo.getTrackJetDistVal());
    // data.fillMulti<float>("trackBTag_JetDistSig", trkinfo.getTrackJetDistSig());    
    
  }
  //const auto &pv = vertices->at(0);
  //std::cout << "primary vertex : " << pv.x() << std::endl;
  
  return true;
}

} /* namespace deepntuples */
