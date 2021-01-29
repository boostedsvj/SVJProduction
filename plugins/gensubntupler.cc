//framework headers
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

//analysis headers
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/getRef.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "fastjet/config.h"

//ROOT headers
#include <TTree.h>
#include <TLorentzVector.h>
 
//STL headers 
#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <utility>
#include <any>
using std::vector;
using std::string;
using std::pair;
using std::map;


class gensubntupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    public:
        explicit gensubntupler(const edm::ParameterSet&);
        ~gensubntupler() {}

    private:
        void beginJob() override;
        void doBeginRun_(const edm::Run&, const edm::EventSetup&) override {}
        void analyze(const edm::Event&, const edm::EventSetup&) override;
        void doEndRun_(const edm::Run&, const edm::EventSetup&) override {}
        void endJob() override {}

        map< string, vector<float> > floatVectorBranches_;
        void makeFloatVectorBranch(string name){
            floatVectorBranches_.emplace(name, vector<float>());
            tree_->Branch(
                (name).c_str(),
                "vector<float>",
                &floatVectorBranches_[name],
                32000, 0
                );            
            }

        map< string, vector<int> > intVectorBranches_;
        void makeIntVectorBranch(string name){
            intVectorBranches_.emplace(name, vector<int>());
            tree_->Branch(
                (name).c_str(),
                "vector<int>",
                &intVectorBranches_[name],
                32000, 0
                );            
            }

        void makeFourVectorBranch(string name){
            makeFloatVectorBranch(name + "_pt");
            makeFloatVectorBranch(name + "_eta");
            makeFloatVectorBranch(name + "_phi");
            makeFloatVectorBranch(name + "_energy");
            makeFloatVectorBranch(name + "_mass");
            }

        void fillFourVectorBranch(string name, const reco::Candidate & candidate){
            floatVectorBranches_[name + "_pt"].push_back(candidate.pt());
            floatVectorBranches_[name + "_eta"].push_back(candidate.eta());
            floatVectorBranches_[name + "_phi"].push_back(candidate.phi());
            floatVectorBranches_[name + "_energy"].push_back(candidate.energy());
            floatVectorBranches_[name + "_mass"].push_back(candidate.mass());
            }

        void clearVectorBranches(){
            for( auto& element : floatVectorBranches_ ) {element.second.clear();}
            for( auto& element : intVectorBranches_ ) {element.second.clear();}
            }

        edm::Service<TFileService> fs;
        TTree* tree_;

        vector<string> jetTags_;
        vector<string> jetTagsWithECFs_;
        map< string, edm::EDGetTokenT<vector<reco::GenJet>> > jetTokens_;

        vector<string> basicJetTags_;
        map< string, edm::EDGetTokenT<vector<reco::BasicJet>> > basicJetTokens_;

        vector<string> valueMapTags_;
        map< string, edm::EDGetTokenT<edm::ValueMap<float>> > valueMapTokens_;
    };


gensubntupler::gensubntupler(const edm::ParameterSet& iConfig) : 
    tree_(NULL),
    jetTags_({
        "ak15GenJetsPackedDark", "ak15GenJetsPackedNoNu",
        "ak15GenJetsAreaDark", "ak15GenJetsAreaNoNu",
        }),
    jetTagsWithECFs_({
        "ak15GenJetsPackedDark", "ak15GenJetsPackedNoNu",
        }),
    basicJetTags_({
        "ak15GenJetsSoftDropDark", "ak15GenJetsSoftDropNoNu",
        }),
    valueMapTags_({
        "ecfN1b1",
        "ecfN1b2",
        "ecfN2b1",
        "ecfN2b2",
        "ecfN3b1",
        "ecfN3b2",
        })
    {
    for (auto & tag : basicJetTags_) {
        basicJetTokens_.emplace(tag, consumes<vector<reco::BasicJet>>(edm::InputTag(tag)));
        jetTokens_.emplace(tag + "_Subjets", consumes<vector<reco::GenJet>>(edm::InputTag(tag, "Subjets")));
        }
    for (auto & tag : jetTags_) {
        jetTokens_.emplace(tag, consumes<vector<reco::GenJet>>(edm::InputTag(tag)));
        }
    for (auto & tag : jetTagsWithECFs_) {
        for (auto & valueMapTag : valueMapTags_) {
            valueMapTokens_.emplace(
                tag + "_" + valueMapTag,
                consumes<edm::ValueMap<float>>(edm::InputTag(tag, valueMapTag))
                );
            }
        }
    usesResource("TFileService");
    }

void gensubntupler::beginJob() {
    tree_ = fs->make<TTree>("tree","tree");
    for (auto & tag : jetTags_) {
        makeFourVectorBranch(tag);
        }
    for (auto & tag : jetTagsWithECFs_) {
        for (auto & valueMapTag : valueMapTags_) {
            makeFloatVectorBranch(tag + "_" + valueMapTag);
            }
        }
    for (auto & tag : basicJetTags_) {
        makeFourVectorBranch(tag);
        makeFourVectorBranch(tag + "_Subjets");
        }
    }

void gensubntupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    clearVectorBranches();

    // Make map tag --> handle
    map< string, edm::Handle<vector<reco::GenJet>> > jetHandles;
    for (auto & tag : jetTags_) {
        edm::EDGetTokenT<vector<reco::GenJet>> token = jetTokens_[tag];
        edm::Handle<vector<reco::GenJet>> handle;
        iEvent.getByToken(token, handle);
        jetHandles.emplace(tag, handle);
        }

    map< string, edm::Handle<vector<reco::BasicJet>> > basicJetHandles;
    for (auto & tag : basicJetTags_) {
        edm::Handle<vector<reco::BasicJet>> handle;
        iEvent.getByToken(basicJetTokens_[tag], handle);
        basicJetHandles.emplace(tag, handle);

        edm::Handle<vector<reco::GenJet>> subjetHandle;
        iEvent.getByToken(jetTokens_[tag + "_Subjets"], subjetHandle);
        jetHandles.emplace(tag + "_Subjets", subjetHandle);
        }

    // Make map (tag, valueMapTag) --> valueMapHandle
    map< string, edm::Handle<edm::ValueMap<float>> > valueMapHandles;
    for (auto & tag : jetTagsWithECFs_) {
        for (auto & valueMapTag : valueMapTags_) {
            iEvent.getByToken(
                valueMapTokens_[tag + "_" + valueMapTag],
                valueMapHandles[tag + "_" + valueMapTag]
                );
            }
        }

    // Fill the basic kinematic variables
    for (auto & tag : jetTags_) {
        for(const auto& jet : *(jetHandles[tag].product())){
            fillFourVectorBranch(tag, jet);
            }
        }

    // Fill the basic kinematic variables for candidates
    for (auto & tag : basicJetTags_) {
        for(const auto& basicJet : *(basicJetHandles[tag].product())){
            fillFourVectorBranch(tag, basicJet);
            }
        }

    // Getting values from the valueMap per jet (little painful)
    for (auto & tag : jetTagsWithECFs_) {
        for (unsigned int i_jet = 0; i_jet < jetHandles[tag].product()->size(); ++i_jet){
            auto jetRef = edm::getRef(jetHandles[tag], i_jet);
            for (auto & valueMapTag : valueMapTags_) {
                floatVectorBranches_[tag + "_" + valueMapTag].push_back(
                    (*valueMapHandles[tag + "_" + valueMapTag])[jetRef]
                    );
                }            
            }
        }

    tree_->Fill();
    }

DEFINE_FWK_MODULE(gensubntupler);