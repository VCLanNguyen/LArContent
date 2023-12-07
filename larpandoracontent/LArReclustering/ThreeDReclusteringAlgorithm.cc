/**
 *  @file   larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.cc
 *
 *  @brief  Implementation file for the reclustering algorithm that runs other algs.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Managers/ClusterManager.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h"

const float convertADCToMeV = 0.0075;
const int atomicNumberArgon = 18;
const float atomicMassArgon = 39.948;
const float criticalEnergyArgon = 32.84;
const float moliereRadiusCmArgon = 9.043; //cm

bool sortByCaloHits (pandora::CaloHitList a, pandora::CaloHitList b) { return (a.size()>b.size()); }

using namespace pandora;

namespace lar_content
{

ThreeDReclusteringAlgorithm::ThreeDReclusteringAlgorithm():
    m_visualDisplaysOn(false),
    m_hitThresholdForNewPfo(0),
    m_peakFinderThreshold(0.1),
    m_writeToTree(false),
    m_treeName("tree"),
    m_fileName("./FoMoutput.root"),
    m_eventId(-1)
{
}

ThreeDReclusteringAlgorithm::~ThreeDReclusteringAlgorithm()
{
       if(m_writeToTree){
	  std::cout <<">>>>>>>>>>>>>>>.Writing to tree"<<std::endl;
          PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
 	}
}

StatusCode ThreeDReclusteringAlgorithm::Run()
{
    ++m_eventId;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventId", m_eventId));

    std::cout << "--------------------------------------------- NEW EVENT ------------------------------------------------" << std::endl;
    if(m_visualDisplaysOn)PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    // Get shower pfos and then find the 3D cluster in the shower pfo.
    const PfoList *pShowerPfoList(nullptr);
    std::string initialPfosListName;
    m_newPfosListNameAllAfterReclustering = "newShowerParticles3D";
    std::string newPfosListNameUnchanged = "unchangedShowerParticles3D";

    if ((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, "ShowerParticles3D", pShowerPfoList)) && pShowerPfoList)
    {
        std::cout << "In this event there are " << pShowerPfoList->size() << " shower pfos." << std::endl;
        if (pShowerPfoList->size()==0) return STATUS_CODE_NOT_FOUND;
        PfoList unchangedPfoList;   

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Pfo>(*this, initialPfosListName));

	int nPfo = pShowerPfoList->size();

    	PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPfo", nPfo));

        int iPfo(0); //for debug
        for (const Pfo *const pShowerPfo : *pShowerPfoList)
        {
            std::cout << std::endl << "-----iPfo = " << iPfo << std::endl;
    	
	    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PfoId", iPfo));

            ClusterList clusterList3D;
            LArPfoHelper::GetThreeDClusterList(pShowerPfo, clusterList3D);

            if(!this->PassesCutsForReclustering(pShowerPfo)) continue; // this function just checks it's a shower at the moment

            // Get the longitudinal and transverse shower profiles
            if (clusterList3D.empty())
                continue;

            CaloHitList caloHitList3D;
            clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

            if(m_visualDisplaysOn)
            {
                PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(),&caloHitList3D, "ShowerPfoCaloHitsBeforeReclustering", RED);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            }

	    int NClustersBefore = clusterList3D.size();
	    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NClustersBefore", NClustersBefore));

            //Quality cuts
            if (caloHitList3D.size() < 2)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            //Some pfos are shower-like and yet include track-like 3D clusters. For the moment I don't want to deal with these.
            const ClusterList *pShowerClusters(nullptr);
            PandoraContentApi::GetList(*this, "ShowerClusters3D", pShowerClusters);
            if(!pShowerClusters) continue;
            if(pShowerClusters->end() == std::find(pShowerClusters->begin(), pShowerClusters->end(), clusterList3D.front())) continue;


            //Free the hits in this cluster, so that they are available for reclustering!
            //Ask to remove the 3D cluster from the parent pfo, so that it's not owned any more
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pShowerPfo, clusterList3D.front()));
            //Pop this cluster in a local clusterlist
            const ClusterList reclusterClusterList(1, clusterList3D.front());
            const TrackList reclusterTrackList; //dummy track list

            // Initialize reclustering with these local lists
            std::string currentClustersListName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentClustersListName));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, "ShowerClusters3D"));

            // Specify clusters and tracks to be used in reclustering
            std::string originalClustersListName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeReclustering(*this, reclusterTrackList, reclusterClusterList, originalClustersListName));

            //Call clustering loop
           std::string reclusterListName;
           ClusterList minimumFigureOfMeritClusterList;
           PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FindBestThreeDClusters(clusterList3D, reclusterListName, minimumFigureOfMeritClusterList));
          //minimumFigureOfMeritClusterList=clusterList3D;

	    int NClustersAfter = minimumFigureOfMeritClusterList.size();
	    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NClustersAfter", NClustersAfter));

           if(minimumFigureOfMeritClusterList==clusterList3D)
           {
               std::cout << "NO CHANGE! Keep the original pfo" << std::endl;
               PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pShowerPfo, clusterList3D.front()));
               PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, originalClustersListName));
               unchangedPfoList.push_back(pShowerPfo);
               
	       int IsRecluster = 0;
	       PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsRecluster", IsRecluster));
           }
           else
           {
               std::cout << "THE CLUSTER LIST CHANGED! I need to make " << minimumFigureOfMeritClusterList.size() << " new pfos." << std::endl;

               PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RebuildPfo(pShowerPfo, minimumFigureOfMeritClusterList));

               const PfoList *pDebugNewPfos(nullptr);
               PandoraContentApi::GetList(*this, m_newPfosListNameAllAfterReclustering, pDebugNewPfos);
		
	       int IsRecluster = 1;
	       PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsRecluster", IsRecluster));
               
	       //if(m_visualDisplaysOn)
               //{
               //   PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
               //   PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(),pDebugNewPfos, "newPfosListNameAllAfterReclustering", RED);
               //   PandoraMonitoringApi::ViewEvent(this->GetPandora());
               //}
           }
            iPfo++;
    	    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
        }

        if(unchangedPfoList.size()>0) PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<PfoList>(*this, "ShowerParticles3D", m_newPfosListNameAllAfterReclustering,  unchangedPfoList));

        const PfoList *pNewPfosListAllAfterReclustering;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList<PfoList>(*this, m_newPfosListNameAllAfterReclustering, pNewPfosListAllAfterReclustering));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_newPfosListNameAllAfterReclustering, "ShowerParticles3D"));
        
        /* 
        const PfoList *pShowerParticles3DDebug;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList<PfoList>(*this, "ShowerParticles3D", pShowerParticles3DDebug));
       
        if(m_visualDisplaysOn)
        {
           PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
           PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(),pShowerParticles3DDebug, "ShowerParticles3DDebug", RED);
           PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }*/
          
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, initialPfosListName));

    }

    //PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

    return STATUS_CODE_SUCCESS;
}

//This can then go in the SDK with the name/format PandoraContentApi::RebuildLArTPCPfo(Algorithm &, Pfo  *pPfoToRebuild, Cluster *pTemplate3DCluster, MapOfListNameInfo &)
StatusCode ThreeDReclusteringAlgorithm::RebuildPfo(const Pfo *pPfoToRebuild, ClusterList &newThreeDClusters)
{
   ClusterList clusterList2D;
   LArPfoHelper::GetTwoDClusterList(pPfoToRebuild, clusterList2D);
   
   std::map<int,const Cluster*> newClustersUMap, newClustersVMap,newClustersWMap;

   for(const Cluster *const pTwoDCluster : clusterList2D)
   {
       /*
       //Visualise 2D cluster before splitting into new 2D clusters
       if(m_visualDisplaysOn)
       {
           ClusterList twoDClusterToRecluster;
           twoDClusterToRecluster.push_back(pTwoDCluster);
           PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
           PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&twoDClusterToRecluster, "twoDClusterToRecluster", RED);
           PandoraMonitoringApi::ViewEvent(this->GetPandora());
       }*/
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfoToRebuild, pTwoDCluster));

       HitType hitType = LArClusterHelper::GetClusterHitType(pTwoDCluster);
       //std::cout << "Currently looking at a 2D cluster with type = " << hitType << std::endl;
       std::string clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

       std::string initialListName="", debugListName="";
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, initialListName));
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, debugListName));

       // Fragmentation initialisation
       std::string originalListName, fragmentListName;
       ClusterList originalClusterList;
       originalClusterList.push_back(pTwoDCluster);


       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

       const OrderedCaloHitList &twoDClusterOrderedCaloHitList(pTwoDCluster->GetOrderedCaloHitList());
       OrderedCaloHitList leftoverCaloHitList = twoDClusterOrderedCaloHitList;
       //std::cout << "Debug size of all hits for this 2D cluster = " << leftoverCaloHitList.size() << std::endl;

       int iCluster(0);
       //std::cout << "debug newThreeDClusters size = " << newThreeDClusters.size() << std::endl;
       for(const Cluster *const pNewCluster : newThreeDClusters)
       {
           PandoraContentApi::Cluster::Parameters parameters;
           CaloHitList newClusterCaloHitList3D;
           pNewCluster->GetOrderedCaloHitList().FillCaloHitList(newClusterCaloHitList3D);                          
           for(const CaloHit *const p3DCaloHit : newClusterCaloHitList3D)
           {
               for (const OrderedCaloHitList::value_type &mapEntry : twoDClusterOrderedCaloHitList)
               {
                   for (const CaloHit *const pCaloHit : *mapEntry.second)
                   {
                       if(pCaloHit==static_cast<const CaloHit *>(p3DCaloHit->GetParentAddress())) 
                       {
                         parameters.m_caloHitList.push_back(static_cast<const CaloHit *>(pCaloHit));
                         leftoverCaloHitList.Remove(pCaloHit);
                       }
                   }
               }
           }
           const Cluster *pNewTwoDCluster(nullptr);
           if (!parameters.m_caloHitList.empty())
           {
               PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewTwoDCluster));
           }
           if(pNewTwoDCluster!=nullptr && !parameters.m_caloHitList.empty() && hitType==TPC_VIEW_U) {newClustersUMap.insert(std::make_pair(iCluster,pNewTwoDCluster)); /*std::cout << "Inserted U cluster for 3D cluster id = " << iCluster << std::endl;*/}
           else if(pNewTwoDCluster!=nullptr && !parameters.m_caloHitList.empty() && hitType==TPC_VIEW_V) {newClustersVMap.insert(std::make_pair(iCluster,pNewTwoDCluster));/*std::cout << "Inserted V cluster for 3D cluster id = " << iCluster << std::endl;*/}

           else if(pNewTwoDCluster!=nullptr && !parameters.m_caloHitList.empty() && hitType==TPC_VIEW_W) {newClustersWMap.insert(std::make_pair(iCluster,pNewTwoDCluster));/*std::cout << "Inserted W cluster for 3D cluster id = " << iCluster << std::endl;*/}

           iCluster++;
       }

       //DEBUG 2D CLUSTERS
       //for (const auto &s : newClustersUMap) std::cout << "U cluster " << s.second << " IsAvailable? " << s.second->IsAvailable() <<std::endl;
       //for (const auto &s : newClustersVMap) std::cout << "V cluster " << s.second << " IsAvailable? " << s.second->IsAvailable() <<std::endl;
       //for (const auto &s : newClustersWMap) std::cout << "W cluster " << s.second << " IsAvailable? " << s.second->IsAvailable() <<std::endl;

       //Visualise 2D cluster after splitting into new 2D clusters
/*       ClusterList reclusteredTwoDClusters;
       if(hitType==TPC_VIEW_U) for (const auto &s : newClustersUMap)   reclusteredTwoDClusters.push_back(s.second);
       if(hitType==TPC_VIEW_V) for (const auto &s : newClustersVMap)   reclusteredTwoDClusters.push_back(s.second);
       if(hitType==TPC_VIEW_W) for (const auto &s : newClustersWMap)   reclusteredTwoDClusters.push_back(s.second);
       
       if(m_visualDisplaysOn)
       {
           PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
           PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&reclusteredTwoDClusters, "reclusteredTwoDClusters", AUTOITER);
           PandoraMonitoringApi::ViewEvent(this->GetPandora());
       }*/
       //Check the leftover caloHits. Attach to the nearest cluster in the new cluster list (newClustersUVect, newClustersVVect or newClustersWVect?
       //std::cout << "Number of leftover hits for this 2D cluster = " << leftoverCaloHitList.size() << std::endl;
       std::map<int,const Cluster*> clustersForLeftoverHitsMap;
       if(hitType==TPC_VIEW_U) clustersForLeftoverHitsMap = newClustersUMap;
       else if(hitType==TPC_VIEW_V) clustersForLeftoverHitsMap = newClustersVMap;
       else if(hitType==TPC_VIEW_W) clustersForLeftoverHitsMap = newClustersWMap;
       //std::cout << "hitType = " << hitType << " newClustersUMap size = " << newClustersUMap.size() << " newClustersVMap = " << newClustersVMap.size() << " newClustersWMap size = " << newClustersWMap.size() << std::endl;
       //Visualise leftover hits
/*       if(m_visualDisplaysOn)
       {
           CaloHitList leftoverCaloHitsToVisualise;
           leftoverCaloHitList.FillCaloHitList(leftoverCaloHitsToVisualise);                          
           PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
           PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(),&leftoverCaloHitsToVisualise, "leftoverCaloHits", AUTOITER);
           PandoraMonitoringApi::ViewEvent(this->GetPandora());
       }*/
       //std::cout << "debug leftover hits 0" << std::endl;
       if(clustersForLeftoverHitsMap.size())  //THE QUESTION REMAINS OF WHAT TO DO WITH LEFTOVER HITS IF TEHRE IS NO CLUSTER TO ATTACH THEM TO (THIS CONDITION FAILS)!!!
       {
       for(const OrderedCaloHitList::value_type &mapEntry : leftoverCaloHitList)
           {
               for (const CaloHit *const pCaloHit : *mapEntry.second)
               {
                   //std::cout << "New leftover 2D hit; now finding nearest cluster..." << std::endl;
                   const Cluster* pNearestCluster = nullptr;
                   double minimumDistance(std::numeric_limits<float>::max());
                   //std::cout << "debug" << std::endl;
                   //std::cout << "pCaloHitPosition = " << pCaloHit->GetPositionVector().GetX() << " " << pCaloHit->GetPositionVector().GetY() << " " << pCaloHit->GetPositionVector().GetZ() << std::endl;
                   //std::cout << "clustersForLeftoverHitsMap size = " << clustersForLeftoverHitsMap.size() << std::endl;
                   for(const auto & [clusterIndex, pNewTwoDCluster] : clustersForLeftoverHitsMap)
                   {
                       //std::cout << "clusterIndex = " << clusterIndex << std::endl;
                       //std::cout << "pNewTwoDCluster = " << pNewTwoDCluster << std::endl;
                       double dist = LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), pNewTwoDCluster);  
                       //std::cout << "distance = " << dist << std::endl;
                       if (dist<minimumDistance)
                       {
                           minimumDistance=dist;
                           pNearestCluster=pNewTwoDCluster;
                       }
                   }
                   //std::cout << "DEBUG minimumDistance = " << minimumDistance << " pNearestCluster = " << pNearestCluster << " is available? " << pNearestCluster->IsAvailable() << " pCaloHit = " << pCaloHit <<  std::endl;
                   PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this,pNearestCluster,pCaloHit));
                   //std::cout << "The nearest cluster has distance = " << minimumDistance << std::endl; 
               }
           }
       }
       //std::cout << "debug leftover hits 1" << std::endl;

       //Visualise 2D clusters after having added in the leftover hits
/*       if(m_visualDisplaysOn)
       {
           PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
           PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&reclusteredTwoDClusters, "reclusteredTwoDClustersWithLeftoverHits", AUTOITER);
           PandoraMonitoringApi::ViewEvent(this->GetPandora());
       }
*/
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, initialListName));

   }

   const PfoList *pNewPfoList(nullptr);
   std::string newPfoListName = "changedShowerParticles3D";
   PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNewPfoList, newPfoListName));

   std::string originalClusterListName="InitialCluster";
   int iCluster(0);

           //std::cout << "newClustersUMap.size = " << newClustersUMap.size() << std::endl;
           //std::cout << "newClustersVMap.size = " << newClustersVMap.size() << std::endl;
           //std::cout << "newClustersWMap.size = " << newClustersWMap.size() << std::endl;


   for(const Cluster *const pNewThreeDCluster : newThreeDClusters)
   {
           //std::cout << "trying to make pfo n. = " << iCluster << " for pNewThreeDCluster = " << pNewThreeDCluster << std::endl;

           PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
           //std::cout << "newClustersUMap.at(iCluster) = " << newClustersUMap.at(iCluster) << std::endl;
           //std::cout << "newClustersVMap.at(iCluster) = " << newClustersVMap.at(iCluster) << std::endl;
           //std::cout << "newClustersWMap.at(iCluster) = " << newClustersWMap.at(iCluster) << std::endl;
           //std::cout << "newClustersWMap.count(iCluster) = " << newClustersWMap.count(iCluster) << std::endl;
           //if(newClustersUMap.at(iCluster)->second) std::cout << "newClustersUMap.at(iCluster)->second->IsAvailable() = " << newClustersUMap.at(iCluster)->second->IsAvailable() << std::endl;
           const bool isAvailableU((newClustersUMap.count(iCluster)) && newClustersUMap.at(iCluster)->IsAvailable());
           const bool isAvailableV((newClustersVMap.count(iCluster)) && newClustersVMap.at(iCluster)->IsAvailable());
           const bool isAvailableW((newClustersWMap.count(iCluster)) && newClustersWMap.at(iCluster)->IsAvailable());
           //std::cout << "isAvailableU = " << isAvailableU << " isAvailableV = " << isAvailableV << " isAvailableW = " << isAvailableW << std::endl;
           CaloHitList clusterUHits, clusterVHits, clusterWHits;
           if(isAvailableU)newClustersUMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterUHits);
           if(isAvailableV)newClustersVMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterVHits);
           if(isAvailableW)newClustersWMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterWHits);
           //std::cout << "debug0" << std::endl;
           if(isAvailableU) pfoParameters.m_clusterList.push_back(newClustersUMap.at(iCluster));
           if(isAvailableV) pfoParameters.m_clusterList.push_back(newClustersVMap.at(iCluster));
           if(isAvailableW) pfoParameters.m_clusterList.push_back(newClustersWMap.at(iCluster));
           //std::cout << "debug1" << std::endl;
           pfoParameters.m_clusterList.push_back(pNewThreeDCluster);

           pfoParameters.m_particleId = pPfoToRebuild->GetParticleId(); // SHOWER, placeholder for now... Are the new clusters all showers???
           pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
           pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
           pfoParameters.m_energy = 0.f;
           pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);

           //std::cout << "debug2" << std::endl;
           const ParticleFlowObject *pNewPfo(nullptr);
           PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pNewPfo));

           ClusterList newClusterList2D;
           //std::cout << "debug3" << std::endl;
           LArPfoHelper::GetTwoDClusterList(pNewPfo, newClusterList2D);

           /*
           if(m_visualDisplaysOn)
           {
              PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
              PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&newClusterList2D, "Pfo2DCluster", RED);
              PandoraMonitoringApi::ViewEvent(this->GetPandora());
           }
	   */

           iCluster++;
           
   }
    std::cout << "Finished making the new pfos" << std::endl;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, newPfoListName, m_newPfosListNameAllAfterReclustering));
    std::cout << "Finished saving the new pfo list" << std::endl;

    return STATUS_CODE_SUCCESS;
}

StatusCode ThreeDReclusteringAlgorithm::FindBestThreeDClusters(ClusterList clusterList3D, std::string &reclusterListName, ClusterList &minimumFigureOfMeritClusterList)
{
    CaloHitList caloHitList3D;
    clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);
    float initialFigureOfMerit=this->GetFigureOfMerit(caloHitList3D);
    if(initialFigureOfMerit<0) {std::cout << "Could not calculate initial FOM!" << std::endl; return STATUS_CODE_FAILURE;} 

    //Call the reclustering algos that produce new cluster candidates
    float minimumFigureOfMerit(initialFigureOfMerit);
    std::cout << "Before Transverse FoM = " << initialFigureOfMerit << std::endl;

    //float RatioFoM = GetTransverseDepthRatioFigureOfMerit(caloHitList3D);
    //std::cout << "Before Depth Ratio FoM = " << RatioFoM << std::endl;

    float longFoM = GetLongitudinalProfileFigureOfMerit(caloHitList3D);
    std::cout << "Before Longitudinal FoM = " << longFoM << std::endl;

    float purity = GetCheatedFigureOfMerit(caloHitList3D);
    std::cout << "Before Purity = " << purity << std::endl;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "beforeTransverseFoM", minimumFigureOfMerit));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "beforeLongFoM", longFoM));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "beforePurity", purity));

    minimumFigureOfMeritClusterList = clusterList3D;
    std::vector<float> mainClusterFractionVector,  initialFigureOfMeritVector, newFigureOfMeritVector, nHitsInitialFomVector, nFinalClustersVector;
    const ClusterList *pReclusterList = NULL;

    //Begin looping through clustering algorithm   
    for (StringVector::const_iterator clusteringIter = m_clusteringAlgorithms.begin(), clusteringIterEnd = m_clusteringAlgorithms.end();
        clusteringIter != clusteringIterEnd; ++clusteringIter)
    {
        // Produce new cluster candidates
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, *clusteringIter, pReclusterList, reclusterListName));

        
        if (pReclusterList->empty())
            continue;
         
        //if(m_visualDisplaysOn)
        //{
        //    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
        //    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),pReclusterList, "ReclusteredClusters", BLUE);
        //    PandoraMonitoringApi::ViewEvent(this->GetPandora());
        //}

        
        //Loop over clusters and calculate new FOM including them if they have more than 10 hits
        std::vector<CaloHitList> newClustersCaloHitLists3D;
        ClusterList clustersFromReclustering;
        for(const Cluster *const pNewCluster : *pReclusterList)
        {
          CaloHitList newClusterCaloHitList3D;
          pNewCluster->GetOrderedCaloHitList().FillCaloHitList(newClusterCaloHitList3D);

          //std::cout << "newClusterCaloHitList3D.size() = " << newClusterCaloHitList3D.size() << std::endl;
          if((int)newClusterCaloHitList3D.size()<m_hitThresholdForNewPfo) continue; //This should remove empty clusters as well (isolated hits problem)

          newClustersCaloHitLists3D.push_back(newClusterCaloHitList3D);
          clustersFromReclustering.push_back(pNewCluster);
         
          //if(m_visualDisplaysOn)
          //{
          //    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &newClusterCaloHitList3D, "newClusterCaloHitList3D", AUTOITER));                 
          //    PandoraMonitoringApi::ViewEvent(this->GetPandora());
          //}

        }


        std::sort (newClustersCaloHitLists3D.begin(), newClustersCaloHitLists3D.end(), sortByCaloHits);
        if(!newClustersCaloHitLists3D.size()) continue; 

        float mainClusterFraction = (float)newClustersCaloHitLists3D.front().size()/caloHitList3D.size();
	
	int nNewClusters = newClustersCaloHitLists3D.size();
	std::cout << "N After clusters = " << nNewClusters << std::endl;

        float newFigureOfMerit = this->GetFigureOfMerit(caloHitList3D, newClustersCaloHitLists3D); //Cheated FOM for group of new clusters is the minimum cheated FOM among those clusters
        std::cout << "After Transverse FoM = " << newFigureOfMerit << std::endl;

        float newlongFoM = GetLongitudinalProfileFigureOfMerit(newClustersCaloHitLists3D);
        std::cout << "After Longitudinal FoM = " << newlongFoM << std::endl;
    
	float newpurity = GetCheatedFigureOfMerit(newClustersCaloHitLists3D);
        std::cout << "After Purity = " << newpurity << std::endl;

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nClustersAfterCheated", nNewClusters));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "afterTransverseFoM", newFigureOfMerit));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "afterLongFoM", newlongFoM));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "afterPurity", newpurity));

        mainClusterFractionVector.push_back(mainClusterFraction);  ///watch out, these are in the loop over many algorithms! if I add more clustering algos I will need to differentiate the entries in these

        //Will print these for study/debug purposes
        initialFigureOfMeritVector.push_back(initialFigureOfMerit);
        newFigureOfMeritVector.push_back(newFigureOfMerit);
        nHitsInitialFomVector.push_back(caloHitList3D.size());
        nFinalClustersVector.push_back(newClustersCaloHitLists3D.size());

        if(newFigureOfMerit<minimumFigureOfMerit)
        {
            minimumFigureOfMerit=newFigureOfMerit;
            minimumFigureOfMeritClusterList=clustersFromReclustering; 
        }
    }
    return STATUS_CODE_SUCCESS;
    //return minimumFigureOfMeritClusterList;
}


//Lateral profile at the shower maximum
float ThreeDReclusteringAlgorithm::GetLateralProfileAtShowerMaximum(float clusterEnergyInMeV, float radiusInCm){
    const float tau=1; //shower maximum
    float energy=clusterEnergyInMeV/criticalEnergyArgon;
    float radius=radiusInCm/moliereRadiusCmArgon;
    float z1 = 0.0251 + 0.00319*std::log(energy);
    float z2 = 0.1162 - 0.000381*atomicNumberArgon;
    float k1 = 0.659 - 0.00309*atomicNumberArgon;
    float k2 = 0.645;
    float k3 = -2.59;
    float k4 = 0.3585 + 0.0421*std::log(energy);
    float p1 = 2.632 - 0.00094*atomicNumberArgon;
    float p2 = 0.401 + 0.00187*atomicNumberArgon;
    float p3 = 1.313 - 0.0686*std::log(energy);

    float Rc = z1+ z2*tau;
    float Rt = k1*(std::exp(k3*(tau-k2))+std::exp(k4*(tau-k2)));
    float prob = p1*std::exp((p2-tau)/p3 - std::exp((p2-tau)/p3));
    float profile = 2*radius*(prob*Rc*Rc/std::pow((radius*radius+Rc*Rc),2)+(1-prob)*Rt*Rt/std::pow((radius*radius+Rt*Rt),2));
    return profile;

}

//Lateral profile at the shower maximum
float ThreeDReclusteringAlgorithm::GetLateralProfileAtShowerDepth(float clusterEnergyInMeV, float radiusInCm, float ratioMaxDepth){
    const float tau=ratioMaxDepth; //ratio to maximum shower depth
    float energy=clusterEnergyInMeV/criticalEnergyArgon;
    float radius=radiusInCm/moliereRadiusCmArgon;
    float z1 = 0.0251 + 0.00319*std::log(energy);
    float z2 = 0.1162 - 0.000381*atomicNumberArgon;
    float k1 = 0.659 - 0.00309*atomicNumberArgon;
    float k2 = 0.645;
    float k3 = -2.59;
    float k4 = 0.3585 + 0.0421*std::log(energy);
    float p1 = 2.632 - 0.00094*atomicNumberArgon;
    float p2 = 0.401 + 0.00187*atomicNumberArgon;
    float p3 = 1.313 - 0.0686*std::log(energy);

    float Rc = z1+ z2*tau;
    float Rt = k1*(std::exp(k3*(tau-k2))+std::exp(k4*(tau-k2)));
    float prob = p1*std::exp((p2-tau)/p3 - std::exp((p2-tau)/p3));
    float profile = 2*radius*(prob*Rc*Rc/std::pow((radius*radius+Rc*Rc),2)+(1-prob)*Rt*Rt/std::pow((radius*radius+Rt*Rt),2));
    return profile;

}


float ThreeDReclusteringAlgorithm::GetFigureOfMerit(std::string figureOfMeritName, CaloHitList mergedClusterCaloHitList3D)
{
    float figureOfMerit(-999);
    if(figureOfMeritName=="cheated") figureOfMerit=this->GetCheatedFigureOfMerit(mergedClusterCaloHitList3D);
    else if(figureOfMeritName=="transversecalo") figureOfMerit=this->GetTransverseProfileFigureOfMerit(mergedClusterCaloHitList3D);
    return figureOfMerit;
}

float ThreeDReclusteringAlgorithm::GetFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{
    std::vector<float> figureOfMeritVector;
    for (StringVector::const_iterator iter = m_figureOfMeritNames.begin(), iterEnd = m_figureOfMeritNames.end(); iter != iterEnd; ++iter)
    {
        figureOfMeritVector.push_back(this->GetFigureOfMerit(*iter,mergedClusterCaloHitList3D)); 
    }
    
    float figureOfMerit=*(std::min_element(figureOfMeritVector.begin(), figureOfMeritVector.end()));
    return figureOfMerit;
}


float ThreeDReclusteringAlgorithm::GetFigureOfMerit(std::string figureOfMeritName, CaloHitList mergedClusterCaloHitList3D, std::vector<CaloHitList> newClustersCaloHitLists3D)
{
    float figureOfMerit(-999);
    if(figureOfMeritName=="cheated")figureOfMerit=this->GetCheatedFigureOfMerit(newClustersCaloHitLists3D);
    else if(figureOfMeritName=="transversecalo") figureOfMerit=this->GetTransverseProfileFigureOfMerit(mergedClusterCaloHitList3D, newClustersCaloHitLists3D);
    return figureOfMerit;
}

float ThreeDReclusteringAlgorithm::GetFigureOfMerit(CaloHitList mergedClusterCaloHitList3D, std::vector<CaloHitList> newClustersCaloHitLists3D)
{
    std::vector<float> figureOfMeritVector;
    for (StringVector::const_iterator iter = m_figureOfMeritNames.begin(), iterEnd = m_figureOfMeritNames.end(); iter != iterEnd; ++iter)
    {
        figureOfMeritVector.push_back(this->GetFigureOfMerit(*iter,mergedClusterCaloHitList3D,newClustersCaloHitLists3D)); 
    }
    
    float figureOfMerit=*(std::min_element(figureOfMeritVector.begin(), figureOfMeritVector.end()));
    return figureOfMerit;
}

float ThreeDReclusteringAlgorithm::GetCheatedFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{
    std::map<int,int> mainMcParticleMap;

    for (const CaloHit *const pCaloHit : mergedClusterCaloHitList3D)
    {
        const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());

        //Find main contributing MC particle Id 
        int mainMcParticleIndex = this->GetMainMcParticleIndex(pParentCaloHit);

        std::map<int, int>::iterator it = mainMcParticleMap.find(mainMcParticleIndex);
 
        if (it != mainMcParticleMap.end()) {
            it->second++;
        }
        else {
            mainMcParticleMap.insert(std::make_pair(mainMcParticleIndex, 1));
        }
    }
    const MCParticleList *pMCParticleList(nullptr);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    MCParticleVector mcParticleVector(pMCParticleList->begin(),pMCParticleList->end());
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

    //std::cout << "Printing map for debug:" << std::endl;
    //for (auto i = mainMcParticleMap.begin(); i!= mainMcParticleMap.end(); i++) {
    //  std::cout << i->first << " : nHits = "<< i->second << " pdg code: " << mcParticleVector.at(i->first)->GetParticleId() << " n. parents = " << mcParticleVector.at(i->first)->GetParentList().size() << "  n. daughters = "  << mcParticleVector.at(i->first)->GetDaughterList().size() << " hierarchy tier = " << this->GetMCParticleHierarchyTier(mcParticleVector.at(i->first)) << " momentum = " << mcParticleVector.at(i->first)->GetMomentum().GetMagnitude() << std::endl;
    //}
    
    auto maxSharedHits = std::max_element(mainMcParticleMap.begin(), mainMcParticleMap.end(), [](const auto &x, const auto &y) {
                    return x.second < y.second;
                });
  
    std::cout << "Cluster size = " << mergedClusterCaloHitList3D.size() << std::endl;

    float mainMcParticleFraction = (float)maxSharedHits->second/mergedClusterCaloHitList3D.size();
    
    std::cout << "Main McParticle Fraction = " << mainMcParticleFraction << std::endl;

    return (mainMcParticleFraction);
    //return (1-mainMcParticleFraction);
}

int ThreeDReclusteringAlgorithm::GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit)
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
    MCParticleVector mcParticleVector(pMCParticleList->begin(),pMCParticleList->end());
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);
    int iMcPart(0); 
    for (const auto &weightMapEntry : pCaloHit->GetMCParticleWeightMap()) 
    { 
      if(weightMapEntry.second>0.5) 
      { 
        iMcPart=0;  
        for(const MCParticle *const pMCParticle: mcParticleVector) 
        { 
            if(pMCParticle==weightMapEntry.first) { break;} 
            iMcPart++; 
        }
      } 
    }
    return iMcPart;
}

int ThreeDReclusteringAlgorithm::GetMCParticleHierarchyTier(const pandora::MCParticle *const pMCParticle)
{
    MCParticleVector mcParticleHierarchyVector;
    int mcParticleHierarchyTier(0);
    if(pMCParticle->GetParentList().size())
    {
        mcParticleHierarchyTier++;
        const pandora::MCParticle *pParentMCParticle = nullptr;
        pParentMCParticle=LArMCParticleHelper::GetParentMCParticle(pMCParticle);
        while(pParentMCParticle->GetParentList().size()){
            pParentMCParticle=LArMCParticleHelper::GetParentMCParticle(pParentMCParticle);
            mcParticleHierarchyTier++;
        }
    }
    return mcParticleHierarchyTier;
}


float ThreeDReclusteringAlgorithm::GetLongitudinalProfileFigureOfMerit(std::vector<CaloHitList> newClustersCaloHitLists3D)
{
    float minimumFigureOfMerit(999999);

    for(CaloHitList newCaloHitList3D: newClustersCaloHitLists3D)
    {
        float currentFigureOfMerit = this->GetLongitudinalProfileFigureOfMerit(newCaloHitList3D);
        if(currentFigureOfMerit<minimumFigureOfMerit)minimumFigureOfMerit=currentFigureOfMerit;
    }

    return minimumFigureOfMerit;
}


float ThreeDReclusteringAlgorithm::GetCheatedFigureOfMerit(std::vector<CaloHitList> newClustersCaloHitLists3D)
{
    //float minimumFigureOfMerit(999999);
    float minimumFigureOfMerit(-999999);

    for(CaloHitList newCaloHitList3D: newClustersCaloHitLists3D)
    {
        float currentFigureOfMerit = this->GetCheatedFigureOfMerit(newCaloHitList3D);
        if(currentFigureOfMerit>minimumFigureOfMerit)minimumFigureOfMerit=currentFigureOfMerit;
        //if(currentFigureOfMerit<minimumFigureOfMerit)minimumFigureOfMerit=currentFigureOfMerit;
    }

    return minimumFigureOfMerit;
}

float ThreeDReclusteringAlgorithm::GetTransverseDepthRatioFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{
    std::cout << ">>>> I'm getting Transverse Depth Ratio Profile FOM BEFORE RECLUSTERING" << std::endl;

    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(mergedClusterCaloHitList3D, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);
 
    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));
    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

    // Now define ortho directions
    const CartesianVector seedDirection((axisDirection.GetX() < std::min(axisDirection.GetY(), axisDirection.GetZ())) ? CartesianVector(1.f, 0.f, 0.f) :
        (axisDirection.GetY() < std::min(axisDirection.GetX(), axisDirection.GetZ())) ? CartesianVector(0.f, 1.f, 0.f) : CartesianVector(0.f, 0.f, 1.f));

    const CartesianVector orthoDirection1(seedDirection.GetCrossProduct(axisDirection).GetUnitVector());

    const CartesianVector orthoDirection2(axisDirection.GetCrossProduct(orthoDirection1).GetUnitVector());
   
    //Observed transverse profile
    int transverseProfileNBins = 50;
    float transverseProfileLow = moliereRadiusCmArgon*-1;
    float transverseProfileHigh = moliereRadiusCmArgon;
    float transverseProfileBinSize = (transverseProfileHigh-transverseProfileLow)/transverseProfileNBins;
    TwoDHistogram observedTransverseProfileMin(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
   
    float minPos = 999999;
    float maxPos = -99999;
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
    {
        const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
        const float position3(hitCoordinate.GetDotProduct(axisDirection));
	if(position3 < minPos) minPos = position3;
	if(position3 > maxPos) maxPos = position3;
     }

    std::cout << "minPos = " << minPos << std::endl;
    std::cout << "maxPos = " << maxPos << std::endl;
    std::cout << "lenth = " << maxPos - minPos << std::endl;

    float showerDepth = (maxPos - minPos); 

    float clusterEnergyInMeV(0);

    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D) clusterEnergyInMeV += pCaloHit3D->GetInputEnergy()*convertADCToMeV;
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
    {
        const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
        const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
        const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
        const float position3(hitCoordinate.GetDotProduct(axisDirection));
        if (position1 < transverseProfileHigh && position1 > transverseProfileLow && position2 <transverseProfileHigh && position2 > transverseProfileLow && position3 < minPos + 30){
          observedTransverseProfileMin.Fill(position1, position2, pCaloHit3D->GetInputEnergy()*convertADCToMeV); // Units: MeV; note counting U, V and W parent hits here!
	  clusterEnergyInMeV += pCaloHit3D->GetInputEnergy()*convertADCToMeV;													 //
	}
     }
     observedTransverseProfileMin.Scale(clusterEnergyInMeV/observedTransverseProfileMin.GetCumulativeSum());

     //A 1 cluster expected profile
     TwoDHistogram expectedTransverseProfile_oneShower(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
     //Expected tranvserse profile (Grindhammer parametrisation)
     for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
     {
         float profileX=transverseProfileLow+iBinX*transverseProfileBinSize;
         for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
         {
            float profileY=transverseProfileLow+iBinY*transverseProfileBinSize;
            float profileRadius=std::sqrt(profileX*profileX+profileY*profileY);
            float profileValue=GetLateralProfileAtShowerDepth(clusterEnergyInMeV,profileRadius, 30/showerDepth);
            expectedTransverseProfile_oneShower.SetBinContent(iBinX, iBinY, profileValue);
         }
     }
     expectedTransverseProfile_oneShower.Scale(clusterEnergyInMeV/expectedTransverseProfile_oneShower.GetCumulativeSum());

     if(m_visualDisplaysOn)
     { 
         std::cout << "Before Clustering: Transverse Profile At Shower Depth Ratio " << 30/showerDepth << std::endl;
         PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedTransverseProfileMin, "COLZ");

         std::cout << "Before Clustering: Expected One Shower Transverse Profile " << std::endl;
         PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedTransverseProfile_oneShower, "COLZ");
     }

     //Calculate figure of merit for this cluster
     float squaredDiffSum(0);
     for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
     {
         for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
         {
           float diff = expectedTransverseProfile_oneShower.GetBinContent(iBinX, iBinY)-observedTransverseProfileMin.GetBinContent(iBinX,iBinY);
           float squaredDiff = diff*diff;     
           squaredDiffSum+=squaredDiff;
         }
     }
     float figureOfMerit = squaredDiffSum/clusterEnergyInMeV; 

     return figureOfMerit;
}

float ThreeDReclusteringAlgorithm::GetTransverseProfileFigureOfMerit(CaloHitList mergedClusterCaloHitList3D, std::vector<CaloHitList> newClustersCaloHitLists3D)
{
    std::cout << ">>>> I'm getting Transverse Profile FOM AFTER RECLUSTERING" << std::endl;

    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(mergedClusterCaloHitList3D, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);
 
    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));
 
    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

    // Now define ortho directions
    const CartesianVector seedDirection((axisDirection.GetX() < std::min(axisDirection.GetY(), axisDirection.GetZ())) ? CartesianVector(1.f, 0.f, 0.f) :
        (axisDirection.GetY() < std::min(axisDirection.GetX(), axisDirection.GetZ())) ? CartesianVector(0.f, 1.f, 0.f) : CartesianVector(0.f, 0.f, 1.f));

    const CartesianVector orthoDirection1(seedDirection.GetCrossProduct(axisDirection).GetUnitVector());

    const CartesianVector orthoDirection2(axisDirection.GetCrossProduct(orthoDirection1).GetUnitVector());
   
    //Estimate total cluster energy
    float clusterEnergyInMeV(0);
    //for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D) clusterEnergyInMeV += pCaloHit3D->GetInputEnergy()*convertADCToMeV;

    //Observed transverse profile
    int transverseProfileNBins = 50;
    float transverseProfileLow = moliereRadiusCmArgon*-1;
    float transverseProfileHigh = moliereRadiusCmArgon;
    float transverseProfileBinSize = (transverseProfileHigh-transverseProfileLow)/transverseProfileNBins;
    TwoDHistogram observedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);

    std::vector<float> caloHitProjectedPosition1Vect, caloHitProjectedPosition2Vect;

    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
    {
        const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
        const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
        const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
        observedTransverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
        //Also save hit projections into vectors
        caloHitProjectedPosition1Vect.push_back(position1);
        caloHitProjectedPosition2Vect.push_back(position2);
    }

    //Scale observed profile to total cluster energy
    observedTransverseProfile.Scale(clusterEnergyInMeV/observedTransverseProfile.GetCumulativeSum());

    //Expected transverse profile
    TwoDHistogram expectedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
    if(newClustersCaloHitLists3D.size()==1 || newClustersCaloHitLists3D.at(0)==mergedClusterCaloHitList3D)
    {
        //Expected tranvserse profile (Grindhammer parametrisation)
        for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
        {
            float profileX=transverseProfileLow+iBinX*transverseProfileBinSize;
            for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
            {
               float profileY=transverseProfileLow+iBinY*transverseProfileBinSize;
               float profileRadius=std::sqrt(profileX*profileX+profileY*profileY);
               float profileValue=GetLateralProfileAtShowerMaximum(clusterEnergyInMeV,profileRadius);
               expectedTransverseProfile.SetBinContent(iBinX, iBinY, profileValue);
            }
        }
        expectedTransverseProfile.Scale(clusterEnergyInMeV/expectedTransverseProfile.GetCumulativeSum());

        if(m_visualDisplaysOn)
        { 
            std::cout << "After Clustering: Expected Transverse Profile " << std::endl;
            PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedTransverseProfile, "COLZ");
        }
    }

       //Now, if I passed a new set of cluster hit lists in the second argument, I should separately calculate all of their transverse profiles projecting them onto the same plane as merged shower.
       else if(newClustersCaloHitLists3D.size()!=1 and newClustersCaloHitLists3D.at(0)!=mergedClusterCaloHitList3D)
     {
        std::cout << "After clustering, # of new clusters = " << newClustersCaloHitLists3D.size() << std::endl;
       
        std::vector<TwoDHistogram> newObservedTransverseProfiles;
        std::vector<double> newClusterEnergies;
        std::vector<float> newClustersCenterPositionsX, newClustersCenterPositionsY;
        for(CaloHitList newCaloHitList3D: newClustersCaloHitLists3D)
        {
            TwoDHistogram newObservedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);;
             
            double clusterEnergy(0);
            for (const CaloHit *const pCaloHit3D : newCaloHitList3D)
            {
                const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
                const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
                const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
                newObservedTransverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
                clusterEnergy+=pCaloHit3D->GetInputEnergy();
            }
            
            newClusterEnergies.push_back(clusterEnergy);
            newClustersCenterPositionsX.push_back(newObservedTransverseProfile.GetMeanX());
            newClustersCenterPositionsY.push_back(newObservedTransverseProfile.GetMeanY());
        }

        //Expected tranvserse profile (Grindhammer parametrisation as a combination of N shower profiles)
        for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
        {
            float profileX=transverseProfileLow+iBinX*transverseProfileBinSize;
            for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
            {
               float profileY=transverseProfileLow+iBinY*transverseProfileBinSize;
               float profileValue(0), shiftedRadius(0), shiftedX(0), shiftedY(0);
               for(std::vector<double>::size_type iCluster=0; iCluster<newClusterEnergies.size(); iCluster++){
                   shiftedX=profileX-newClustersCenterPositionsX.at(iCluster);
                   shiftedY=profileY-newClustersCenterPositionsY.at(iCluster);
                   shiftedRadius=std::sqrt(shiftedX*shiftedX+shiftedY*shiftedY);
                   profileValue+=GetLateralProfileAtShowerMaximum(newClusterEnergies.at(iCluster)*convertADCToMeV,shiftedRadius);
               }
               expectedTransverseProfile.SetBinContent(iBinX, iBinY, profileValue);
            }
        }

        expectedTransverseProfile.Scale(clusterEnergyInMeV/expectedTransverseProfile.GetCumulativeSum());

        if(m_visualDisplaysOn)
        { 
            std::cout << "After reclustering: Expected transverse profile of individual cluster " <<std::endl;
            PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedTransverseProfile, "COLZ");
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }
    }

    //Calculate figure of merit for this cluster
    float squaredDiffSum(0);
    for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
    {
        for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
        {
          float diff = expectedTransverseProfile.GetBinContent(iBinX, iBinY)-observedTransverseProfile.GetBinContent(iBinX,iBinY);
          float squaredDiff = diff*diff;     
          squaredDiffSum+=squaredDiff;
        }
    }
    float figureOfMerit = squaredDiffSum/clusterEnergyInMeV; 

    return figureOfMerit;
}

float ThreeDReclusteringAlgorithm::GetTransverseProfileFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{
    std::cout << ">>>> I'm getting Transverse Profile FOM BEFORE RECLUSTERING" << std::endl;

    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(mergedClusterCaloHitList3D, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);
 
    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));
    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

    // Now define ortho directions
    const CartesianVector seedDirection((axisDirection.GetX() < std::min(axisDirection.GetY(), axisDirection.GetZ())) ? CartesianVector(1.f, 0.f, 0.f) :
        (axisDirection.GetY() < std::min(axisDirection.GetX(), axisDirection.GetZ())) ? CartesianVector(0.f, 1.f, 0.f) : CartesianVector(0.f, 0.f, 1.f));

    const CartesianVector orthoDirection1(seedDirection.GetCrossProduct(axisDirection).GetUnitVector());

    const CartesianVector orthoDirection2(axisDirection.GetCrossProduct(orthoDirection1).GetUnitVector());
   
    //Estimate total cluster energy
    float clusterEnergyInMeV(0);
    //for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D) clusterEnergyInMeV += pCaloHit3D->GetInputEnergy()*convertADCToMeV;

    //Observed transverse profile
    int transverseProfileNBins = 50;
    float transverseProfileLow = moliereRadiusCmArgon*-1;
    float transverseProfileHigh = moliereRadiusCmArgon;
    float transverseProfileBinSize = (transverseProfileHigh-transverseProfileLow)/transverseProfileNBins;
    TwoDHistogram observedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);

    std::vector<float> caloHitProjectedPosition1Vect, caloHitProjectedPosition2Vect, caloHitEnergyVect;

    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
    {
        const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
        const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
        const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
        if (position1 < transverseProfileHigh && position1 > transverseProfileLow && position2 <transverseProfileHigh && position2 > transverseProfileLow){
          observedTransverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()*convertADCToMeV); // Units: MeV; note counting U, V and W parent hits here!
          
          //Also save hit projections into vectors
          caloHitProjectedPosition1Vect.push_back(position1);
          caloHitProjectedPosition2Vect.push_back(position2);
         
          // get energy vector for 2shower prediction
	  caloHitEnergyVect.push_back(pCaloHit3D->GetInputEnergy()*convertADCToMeV);
	  clusterEnergyInMeV += pCaloHit3D->GetInputEnergy()*convertADCToMeV;
        }
     }

     //Scale observed profile to total cluster energy
     observedTransverseProfile.Scale(clusterEnergyInMeV/observedTransverseProfile.GetCumulativeSum());

     //A 1 cluster expected profile
     TwoDHistogram expectedTransverseProfile_oneShower(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
     //Expected tranvserse profile (Grindhammer parametrisation)
     for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
     {
         float profileX=transverseProfileLow+iBinX*transverseProfileBinSize;
         for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
         {
            float profileY=transverseProfileLow+iBinY*transverseProfileBinSize;
            float profileRadius=std::sqrt(profileX*profileX+profileY*profileY);
            float profileValue=GetLateralProfileAtShowerMaximum(clusterEnergyInMeV,profileRadius);
            expectedTransverseProfile_oneShower.SetBinContent(iBinX, iBinY, profileValue);
         }
     }
     expectedTransverseProfile_oneShower.Scale(clusterEnergyInMeV/expectedTransverseProfile_oneShower.GetCumulativeSum());
    
     //------------------- Andy Chappel Implementation -------------//
     //Canvas canvas( transverseProfileNBins, transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileLow, transverseProfileHigh); 
    
     ////Fill canvas from 2D observed histogram   
     //for (int xp = 0; xp < observedTransverseProfile.GetNBinsX(); ++xp)
     //{
     //  for (int yp = 0; yp < observedTransverseProfile.GetNBinsY(); ++yp)
     //  {
     //    canvas.m_canvas[yp][xp] = observedTransverseProfile.GetBinContent(xp, yp);      
     //  }
     //}
 
     //CartesianPointVector vertices;
     //std::vector<float> energies;
     //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetVerticesFromCanvas(canvas, vertices, energies));
     //
     //TwoDHistogram peakFinderTransverseProfileCenters(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
     //
     //std::vector<float> caloHitProjectedPosition1VectPeak, caloHitProjectedPosition2VectPeak, caloHitEnergyVectPeak;

     //for (unsigned int i = 0; i < vertices.size(); i++)
     //{
     //  std::cout << "x = " <<vertices[i].GetX() << " y = " << vertices[i].GetY() << " E = " << energies[i] << std::endl;
     //  peakFinderTransverseProfileCenters.Fill(vertices[i].GetX(), vertices[i].GetY(), energies[i]);
     //  caloHitProjectedPosition1VectPeak.push_back(vertices[i].GetX());
     //  caloHitProjectedPosition2VectPeak.push_back(vertices[i].GetY());
     //  caloHitEnergyVectPeak.push_back(energies[i]);
     //}
     //
     //CartesianVector center1Peak(0,0,0), center2Peak(0,0,0);
     //float energy1Peak(-999999), energy2Peak(-999999);

     //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetKmeansPredictions(caloHitProjectedPosition1VectPeak, caloHitProjectedPosition2VectPeak, caloHitEnergyVectPeak, center1Peak, center2Peak, energy1Peak, energy2Peak)); 
     //std::cout << "center1 = " << center1Peak << " center2 = " << center2Peak << " energy1 = " << energy1Peak << " energy2 = " << energy2Peak << std::endl;

     ////Just visualise transverse profile for debugging kMeans
     //TwoDHistogram observedTransverseProfileCentersPeak(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
     //observedTransverseProfileCentersPeak.Fill(center1Peak.GetX(),center1Peak.GetY(),energy1Peak);
     //observedTransverseProfileCentersPeak.Fill(center2Peak.GetX(),center2Peak.GetY(),energy2Peak);
     //
     ////A 2 cluster expected profile
     //TwoDHistogram peakFinderExpectedTransverseProfile_twoShowers(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
     //
     ////Expected tranvserse profile (Grindhammer parametrisation)
     //for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
     //{
     //    float profileX=transverseProfileLow+iBinX*transverseProfileBinSize;
     //    for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
     //    {
     //       float profileY=transverseProfileLow+iBinY*transverseProfileBinSize;
     //       float shiftedX1=profileX-center1Peak.GetX();
     //       float shiftedY1=profileY-center1Peak.GetY();
     //       float shiftedRadius1=std::sqrt(shiftedX1*shiftedX1+shiftedY1*shiftedY1);

     //       float shiftedX2=profileX-center2Peak.GetX();
     //       float shiftedY2=profileY-center2Peak.GetY();
     //       float shiftedRadius2=std::sqrt(shiftedX2*shiftedX2+shiftedY2*shiftedY2);

     //       float profileValue=GetLateralProfileAtShowerMaximum(energy1Peak,shiftedRadius1)+GetLateralProfileAtShowerMaximum(energy2Peak,shiftedRadius2);
     //       //float profileValue=GetLateralProfileAtShowerMaximum(clusterEnergyInMeV,shiftedRadius1)+GetLateralProfileAtShowerMaximum(clusterEnergyInMeV,shiftedRadius2);
     //       peakFinderExpectedTransverseProfile_twoShowers.SetBinContent(iBinX, iBinY, profileValue);
     //    }
     //}

     //peakFinderExpectedTransverseProfile_twoShowers.Scale(clusterEnergyInMeV/peakFinderExpectedTransverseProfile_twoShowers.GetCumulativeSum());
     //------------------- End of Andy Chappel Implementation -------------//

     //-------------------kMeans Implementation for 2 showers predictions -------------//

     //Now I want to calculate a 2 cluster expected profile, calculate FOM in both cases, and take the ratio, and this will be the final FOM
     //I need to use something like kMeans to predict shower centers and energies based on observed profile!
     CartesianVector center1(0,0,0), center2(0,0,0);
     float energy1(-999999), energy2(-999999);

     PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetKmeansPredictions(caloHitProjectedPosition1Vect, caloHitProjectedPosition2Vect, caloHitEnergyVect, center1, center2, energy1, energy2)); 
     //std::cout << "center1 = " << center1 << " center2 = " << center2 << " energy1 = " << energy1 << " energy2 = " << energy2 << std::endl;

     //Just visualise transverse profile for debugging kMeans
     TwoDHistogram observedTransverseProfileCenters(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
     observedTransverseProfileCenters.Fill(center1.GetX(),center1.GetY(),energy1);
     observedTransverseProfileCenters.Fill(center2.GetX(),center2.GetY(),energy2);

     //A 2 cluster expected profile
     TwoDHistogram expectedTransverseProfile_twoShowers(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
     
     //Expected tranvserse profile (Grindhammer parametrisation)
     for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
     {
         float profileX=transverseProfileLow+iBinX*transverseProfileBinSize;
         for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
         {
            float profileY=transverseProfileLow+iBinY*transverseProfileBinSize;
	    float shiftedX1=profileX-center1.GetX();
            float shiftedY1=profileY-center1.GetY();
            float shiftedRadius1=std::sqrt(shiftedX1*shiftedX1+shiftedY1*shiftedY1);

	    float shiftedX2=profileX-center2.GetX();
            float shiftedY2=profileY-center2.GetY();
            float shiftedRadius2=std::sqrt(shiftedX2*shiftedX2+shiftedY2*shiftedY2);

            float profileValue=GetLateralProfileAtShowerMaximum(energy1,shiftedRadius1)+GetLateralProfileAtShowerMaximum(energy2,shiftedRadius2);
            //float profileValue=GetLateralProfileAtShowerMaximum(clusterEnergyInMeV,shiftedRadius1)+GetLateralProfileAtShowerMaximum(clusterEnergyInMeV,shiftedRadius2);
            expectedTransverseProfile_twoShowers.SetBinContent(iBinX, iBinY, profileValue);
         }
     }

     expectedTransverseProfile_twoShowers.Scale(clusterEnergyInMeV/expectedTransverseProfile_twoShowers.GetCumulativeSum());
     
     //------------------kMeans Implementation ends-------------------// 
     if(m_visualDisplaysOn)
     { 
        std::cout << std::endl;
        std::cout << "Observed Transverse Profile" << std::endl;
        PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedTransverseProfile, "COLZ");

        std::cout << "Expected 1 shower Transverse Profile" << std::endl;
        PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedTransverseProfile_oneShower, "COLZ");

        std::cout << "kMeans 2 centers" << std::endl;
        PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedTransverseProfileCenters, "COLZ");

        std::cout << "Expected 2 showers Transverse Profile" << std::endl;
        PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedTransverseProfile_twoShowers, "COLZ");

    	//std::cout << "Peak Finder centers" << std::endl;
        //PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), peakFinderTransverseProfileCenters, "COLZ");

        //std::cout << "Peak Finder + kMeans 2 centers" << std::endl;
        //PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedTransverseProfileCentersPeak, "COLZ");
        //
        //std::cout << "Peak Finder Expected 2 showers Transverse Profile" << std::endl;
        //PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), peakFinderExpectedTransverseProfile_twoShowers, "COLZ");
        std::cout << std::endl;
     }
 
     //Calculate figure of merit for this cluster
     float squaredDiffSumOne(0);
     float squaredDiffSumTwo(0);
     //float squaredDiffSumTwoPeak(0);
     for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
     {
         for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
         {
           float diffOne = expectedTransverseProfile_oneShower.GetBinContent(iBinX, iBinY)-observedTransverseProfile.GetBinContent(iBinX,iBinY);
           float diffTwo = expectedTransverseProfile_twoShowers.GetBinContent(iBinX, iBinY)-observedTransverseProfile.GetBinContent(iBinX,iBinY);
           //float diffTwoPeak = peakFinderExpectedTransverseProfile_twoShowers.GetBinContent(iBinX, iBinY)-observedTransverseProfile.GetBinContent(iBinX,iBinY);

           float squaredDiffOne = diffOne*diffOne;     
           float squaredDiffTwo = diffTwo*diffTwo;     
           //float squaredDiffTwoPeak = diffTwoPeak*diffTwoPeak;     

           squaredDiffSumOne+=squaredDiffOne;
           squaredDiffSumTwo+=squaredDiffTwo;
           //squaredDiffSumTwoPeak+=squaredDiffTwoPeak;
         }
     }

     float figureOfMeritOne = squaredDiffSumOne/clusterEnergyInMeV;   
     float figureOfMeritTwo = squaredDiffSumTwo/clusterEnergyInMeV;   
     //float figureOfMeritTwoPeak = squaredDiffSumTwoPeak/clusterEnergyInMeV;   
   
     std::cout << "Before Clustering: FoM One Shower = " << figureOfMeritOne << std::endl;
     std::cout << "Before Clustering: FoM Two Shower = " << figureOfMeritTwo << std::endl;
     //std::cout << "Before Clustering: FoM Two Shower Peak Finder = " << figureOfMeritTwoPeak << std::endl;
 
     return figureOfMeritTwo;
}

float ThreeDReclusteringAlgorithm::GetLongitudinalProfileFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{
    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(mergedClusterCaloHitList3D, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));

    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

   // Observed longitudinal profile
    Histogram observedLongitudinalProfile(140, 0., 140.);

    const float convertGeVToMeV(1000.f);
    const float convertCmToX0(1.f / 14.f);
    float clusterEnergyInMeV(0.f);

    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
    {
        const CaloHit *const pParentCaloHit(static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress()));

        if (TPC_VIEW_W != pParentCaloHit->GetHitType())
            continue;

        clusterEnergyInMeV += convertADCToMeV * pParentCaloHit->GetInputEnergy(); // Used later on

        const float longitudinalCoordInCm((pCaloHit3D->GetPositionVector() - axisIntercept).GetDotProduct(axisDirection));
        observedLongitudinalProfile.Fill(longitudinalCoordInCm * convertCmToX0, convertADCToMeV * pParentCaloHit->GetInputEnergy());
    }

    // Expected longitudinal profile
    Histogram expectedLongitudinalProfile(140, 0., 140.);

    const float clusterEnergyInGeV(clusterEnergyInMeV / convertGeVToMeV);
    const float longProfileCriticalEnergy(0.08f);
    const float longProfileParameter0(1.25f);
    const float longProfileParameter1(0.5f);

    const double a(longProfileParameter0 + longProfileParameter1 * std::log(clusterEnergyInGeV / longProfileCriticalEnergy));
    const double gammaA(std::exp(lgamma(a)));

    float t(0.f);
    for (int iBin = 0; iBin < expectedLongitudinalProfile.GetNBinsX(); ++iBin)
    {
        t += expectedLongitudinalProfile.GetXBinWidth();
        expectedLongitudinalProfile.Fill(t, convertGeVToMeV * clusterEnergyInGeV / 2. * std::pow(t / 2.f, static_cast<float>(a - 1.)) *
            std::exp(-t / 2.) * expectedLongitudinalProfile.GetXBinWidth() / gammaA);
    }

    //if(m_visualDisplaysOn) {
    //    std::cout << "Observed longitudinal energy profile " << std::endl;
    //    PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedLongitudinalProfile, "");
    //    std::cout << "Expected longitudinal energy profile " << std::endl;
    //    PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedLongitudinalProfile, "");
    //}

    //Calculate figure of merit for this cluster
    float squaredDiffSum(0);
    for (int iBinX=0; iBinX < observedLongitudinalProfile.GetNBinsX(); iBinX++) 
    {
          float diff = expectedLongitudinalProfile.GetBinContent(iBinX)-observedLongitudinalProfile.GetBinContent(iBinX);
          float squaredDiff = diff*diff;     
          squaredDiffSum+=squaredDiff;
    }

    float figureOfMerit = squaredDiffSum/clusterEnergyInMeV;   

    return figureOfMerit; //placeholder for fom

}

bool ThreeDReclusteringAlgorithm::PassesCutsForReclustering(const pandora::ParticleFlowObject *const pShowerPfo)
{
    if (LArPfoHelper::IsShower(pShowerPfo)) return true;
    return false;
}

////Andy Chappel Find Peak code
//StatusCode ThreeDReclusteringAlgorithm::GetVerticesFromCanvas(const Canvas &canvas, CartesianPointVector &vertices, std::vector<float> &energies) const
//{
//    //dx, dy should be the bin width of 2d hist
//    const double dx{2}, dy{2};
//    //need to find 2 peaks, with some seperation
//
//    float maxIntensity{0.f};
//    for (int xp = 0; xp < canvas.m_width; ++xp)
//    {
//        for (int yp = 0; yp < canvas.m_height; ++yp)
//        {
//            const float localIntensity{canvas.m_canvas[yp][xp]};
//            if (localIntensity > maxIntensity)
//                maxIntensity = localIntensity;
//        }
//    }
//    const float threshold{maxIntensity * m_peakFinderThreshold};
//
//    std::vector<std::vector<std::pair<int, int>>> peaks;
//    for (int xp = 0; xp < canvas.m_width; ++xp)
//    {
//        for (int yp = 0; yp < canvas.m_height; ++yp)
//        {
//            const float localIntensity{canvas.m_canvas[yp][xp]};
//
//            std::vector<std::pair<int, int>> peak;
//            bool hasLowNeighbour{false};
//            for (int dr = -1; dr <= 1; ++dr)
//            {
//                for (int dc = -1; dc <=1; ++dc)
//                {
//                    if (dr == 0 && dc == 0)
//                        continue;
//                    const int r{yp + dr}, c{xp + dc};
//                    if (r < 0 || r >= canvas.m_height || c < 0 || c >= canvas.m_width)
//                        continue;
//
//                    const float neighborIntensity{canvas.m_canvas[r][c]};
//                    if (localIntensity > neighborIntensity)
//                    {
//                        hasLowNeighbour = true;
//                    }
//                    else if (localIntensity < neighborIntensity)
//                    {
//                        hasLowNeighbour = false;
//                        break;
//                    }
//                }
//            }
//            if (hasLowNeighbour && localIntensity > threshold)
//                this->GrowPeak(canvas, xp, yp, localIntensity, peak);
//            if (!peak.empty())
//	    {
//                peaks.emplace_back(peak);
//	    }
//        }
//    }
//
//    std::cout << "Found #peaks = " << peaks.size() << std::endl;
//    
//    for (const auto peak : peaks)
//    {
//        int row{0}, col{0};
//        for (const auto pixel : peak)
//        {
//            row += pixel.second;
//            col += pixel.first;
//        }
//        row /= peak.size();
//        col /= peak.size();
//
//        const float x{static_cast<float>(col * dx + canvas.m_xMin)};
//        const float y{static_cast<float>(row * dy + canvas.m_yMin)};
//        const float e{canvas.m_canvas[row][col]};
//
//        //std::cout << "x = " << x << " y = " << y << std::endl;
//        //std::cout << "energy = " << e << std::endl;
//
//        CartesianVector pt(x, y, 0);
//        vertices.emplace_back(pt);
//        energies.emplace_back(e);
//    }
//
//    return STATUS_CODE_SUCCESS;
//}
//
//bool ThreeDReclusteringAlgorithm::GrowPeak(const Canvas &canvas, int col, int row, float intensity, std::vector<std::pair<int, int>> &peak) const
//{
//    if (col < 0 || col >= canvas.m_width || row < 0 || row >= canvas.m_height || canvas.m_visited[row][col] || canvas.m_canvas[row][col] < intensity)
//        return false;
//
//    // Check that no adjacent pixel is larger than this one
//    for (int dc = -1; dc <= 1; ++dc)
//    {
//        const int c{col + dc};
//        if (c < 0 || c >= canvas.m_width)
//            continue;
//
//        for (int dr = -1; dr <= 1; ++dr)
//        {
//            const int r{row + dr};
//            if (r < 0 || r >= canvas.m_height)
//                continue;
//
//            if (dr == 0 && dc == 0)
//                continue;
//
//            const float neighborIntensity{canvas.m_canvas[r][c]};
//            if (neighborIntensity > intensity)
//                return false;
//        }
//    }
//
//    // Need to check we aren't growing into a higher peak, if we are restart from the current pixel
//    float localIntensity{canvas.m_canvas[row][col]};
//    if (localIntensity > intensity)
//    {
//        intensity = localIntensity;
//        for (const auto pixel : peak)
//            canvas.m_visited[pixel.second][pixel.first] = false;
//        peak.clear();
//        this->GrowPeak(canvas, col, row, intensity, peak);
//        return true;
//    }
//
//    // Add pixel to the peak
//    canvas.m_visited[row][col] = true;
//    peak.emplace_back(std::make_pair(col, row));
//
//    for (int dc = -1; dc <= 1; ++dc)
//    {
//        for (int dr = -1; dr <= 1; ++dr)
//        {
//            if (dr == 0 && dc == 0)
//                continue;
//            bool reset{this->GrowPeak(canvas, col + dc, row + dr, intensity, peak)};
//            // If we started growing a non-peak region, stop looking relative to the previous peak
//            if (reset)
//                return reset;
//        }
//    }
//
//    return false;
//}

//A kmeans algorithm implementation from the internet. For now K=2!!
StatusCode ThreeDReclusteringAlgorithm::GetKmeansPredictions(std::vector<float> caloHitProjectedPosition1Vect, std::vector<float> caloHitProjectedPosition2Vect, std::vector<float> caloHitEnergyVect, CartesianVector &center1, CartesianVector &center2, float &energy1, float &energy2)
{

    //int max_iterations(20);
    //Use the value of K random data points as the centers
    std::vector<float> InitialCenterPositions1, InitialCenterPositions2;
    int rand1 = rand()%caloHitProjectedPosition1Vect.size();
    InitialCenterPositions1.push_back(caloHitProjectedPosition1Vect.at(rand1));
    InitialCenterPositions2.push_back(caloHitProjectedPosition2Vect.at(rand1));
    int rand2 = rand1;
    while(rand2 == rand1) rand2 = rand()%caloHitProjectedPosition1Vect.size();
    InitialCenterPositions1.push_back(caloHitProjectedPosition1Vect.at(rand2));
    InitialCenterPositions2.push_back(caloHitProjectedPosition2Vect.at(rand2));

    //float InitialError(std::numeric_limits<float>::max());
    float error(999999), centerDistance(999999)/*, cluster1Position1Sum(0), cluster1Position2Sum(0),cluster2Position1Sum(0),cluster2Position2Sum(0)*/;
    float cluster1Position1(0), cluster1Position2(0), cluster2Position1(0), cluster2Position2(0);
    int nIterations(0);

    while(error > 2 && nIterations < 20 && centerDistance > 0.001)
    {
        //std::cout << "error at start of while loop = " << error << std::endl;
        error=0; centerDistance=0; energy1=0; energy2=0;
        cluster1Position1=0;
        cluster1Position2=0;
        cluster2Position1=0;
        cluster2Position2=0;
        //Add points to nearest center
        //std::vector<float> cluster1Position1, cluster1Position2, cluster2Position1, cluster2Position2;
        int nHitsInCluster1(0), nHitsInCluster2(0);
        for(long unsigned int iPoint=0; iPoint<caloHitProjectedPosition1Vect.size(); iPoint++)
        {
          float distFromCenter1 = std::sqrt((caloHitProjectedPosition1Vect.at(iPoint)-InitialCenterPositions1.at(0))*(caloHitProjectedPosition1Vect.at(iPoint)-InitialCenterPositions1.at(0))+(caloHitProjectedPosition2Vect.at(iPoint)-InitialCenterPositions2.at(0))*(caloHitProjectedPosition2Vect.at(iPoint)-InitialCenterPositions2.at(0)));
          
          float distFromCenter2 = std::sqrt((caloHitProjectedPosition1Vect.at(iPoint)-InitialCenterPositions1.at(1))*(caloHitProjectedPosition1Vect.at(iPoint)-InitialCenterPositions1.at(1))+(caloHitProjectedPosition2Vect.at(iPoint)-InitialCenterPositions2.at(1))*(caloHitProjectedPosition2Vect.at(iPoint)-InitialCenterPositions2.at(1)));

          //std::cout << "dist1 = " << distFromCenter1 << " dist2 = " << distFromCenter2 << std::endl;
          if(distFromCenter1<distFromCenter2)
          {
              error += distFromCenter1*distFromCenter1;
              //I sum positions here so that I can calculate new centers later
              //cluster1Position1+=caloHitProjectedPosition1Vect.at(iPoint);
              //cluster1Position2+=caloHitProjectedPosition2Vect.at(iPoint);
              //energy1+=caloHitEnergyVect.at(iPoint);
              //nHitsInCluster1++;
           
              //try weighted by energy of that position
              int weight = caloHitEnergyVect.at(iPoint);   
              cluster1Position1+=caloHitProjectedPosition1Vect.at(iPoint)*weight;
              cluster1Position2+=caloHitProjectedPosition2Vect.at(iPoint)*weight;
              energy1+=caloHitEnergyVect.at(iPoint);
              nHitsInCluster1 += weight;
          }
          else
          {
              error += distFromCenter2*distFromCenter2;
              //I sum positions here so that I can calculate new centers later
              //cluster2Position1+=caloHitProjectedPosition1Vect.at(iPoint);
              //cluster2Position2+=caloHitProjectedPosition2Vect.at(iPoint);
              //energy2+=caloHitEnergyVect.at(iPoint);
              //nHitsInCluster2++;
              
	      //try weighted by energy of that position
              int weight = caloHitEnergyVect.at(iPoint);   
              cluster2Position1+=caloHitProjectedPosition1Vect.at(iPoint)*weight;
              cluster2Position2+=caloHitProjectedPosition2Vect.at(iPoint)*weight;
              energy2+=caloHitEnergyVect.at(iPoint);
              nHitsInCluster2 += weight;
          }
        
        //std::cout << "iPoint = " << iPoint << std::endl;
        //std::cout << "og x = " << caloHitProjectedPosition1Vect.at(iPoint) << " og y = " << caloHitProjectedPosition1Vect.at(iPoint) << " og e = " << caloHitEnergyVect.at(iPoint) << std::endl;
        //std::cout << "cluster1 x = " << cluster1Position1 << " cluster1 y = " << cluster1Position2 << " nHits cluster1 = " << nHitsInCluster1 << std::endl;
        //std::cout << "cluster2 x = " << cluster2Position1 << " cluster2 y = " << cluster2Position2 << " nHits cluster2 = " << nHitsInCluster2 << std::endl;
        }
    
        //Update new clusters centers
        cluster1Position1=cluster1Position1/nHitsInCluster1;
        cluster1Position2=cluster1Position2/nHitsInCluster1;
        cluster2Position1=cluster2Position1/nHitsInCluster2;
        cluster2Position2=cluster2Position2/nHitsInCluster2;

        float centerDistance1=std::sqrt((InitialCenterPositions1.at(0)-cluster1Position1)*(InitialCenterPositions1.at(0)-cluster1Position1)+(InitialCenterPositions1.at(1)-cluster1Position2)*(InitialCenterPositions1.at(1)-cluster1Position2));
        float centerDistance2=std::sqrt((InitialCenterPositions2.at(0)-cluster2Position1)*(InitialCenterPositions2.at(0)-cluster2Position1)+(InitialCenterPositions2.at(1)-cluster2Position2)*(InitialCenterPositions2.at(1)-cluster2Position2));
        centerDistance=std::max(centerDistance1,centerDistance2);

        InitialCenterPositions1.at(0)=cluster1Position1;
        InitialCenterPositions1.at(1)=cluster1Position2;
        InitialCenterPositions2.at(0)=cluster2Position1;
        InitialCenterPositions2.at(1)=cluster2Position2;
        
        //std::cout <<"Iteration " << nIterations << ": Error = " << error << " new cluster 1 center = " << InitialCenterPositions1.at(0) << " " << InitialCenterPositions1.at(1) << " new cluster 2 center = " << InitialCenterPositions2.at(0) << " " << InitialCenterPositions2.at(1) << " maximum distance between old and new centers = " << centerDistance << std::endl;
        //std::cout << std::endl;

        nIterations++;   
    }
        //std::cout << "AT THE END Error = " << error << " new cluster 1 center = " << InitialCenterPositions1.at(0) << " " << InitialCenterPositions1.at(1) << " new cluster 2 center = " << InitialCenterPositions2.at(0) << " " << InitialCenterPositions2.at(1) << " maximum distance between old and new centers = " << centerDistance << std::endl;

        center1.SetValues(InitialCenterPositions1.at(0),InitialCenterPositions1.at(1),0);
        center2.SetValues(InitialCenterPositions2.at(0),InitialCenterPositions2.at(1),0);
    
    //std::cout << "------end kmean calculation-------"<<std::endl;
    //std::cout << std::endl;

   return STATUS_CODE_SUCCESS;
}

StatusCode ThreeDReclusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "FigureOfMeritNames", m_figureOfMeritNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "HitThresholdForNewPfo", m_hitThresholdForNewPfo));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualDisplaysOn", m_visualDisplaysOn));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PeakFinderThreshold", m_peakFinderThreshold));
    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "WriteToTree", m_writeToTree));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "ClusteringAlgorithms", m_clusteringAlgorithms));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

//ThreeDReclusteringAlgorithm::Canvas::Canvas(const int width, const int height, const float xMin, const float xMax, const float yMin, const float yMax) :
//    m_width{width},
//    m_height{height},
//    m_xMin{xMin},
//    m_xMax{xMax},
//    m_yMin{yMin},
//    m_yMax{yMax}
//{
//    m_canvas = new float*[m_height];
//    m_visited = new bool*[m_height];
//    for (int r = 0; r < m_height; ++r)
//    {
//        m_canvas[r] = new float[m_width]{};
//        m_visited[r] = new bool[m_width]{false};
//    }
//}
//
////-----------------------------------------------------------------------------------------------------------------------------------------
//
//ThreeDReclusteringAlgorithm::Canvas::~Canvas()
//{
//    for (int r = 0; r < m_height; ++r)
//    {
//        delete[] m_canvas[r];
//        delete[] m_visited[r];
//    }
//    delete[] m_canvas;
//    delete[] m_visited;
//}

//-----------------------------------------------------------------------------------------------------------------------------------------

}// namespace lar_content
