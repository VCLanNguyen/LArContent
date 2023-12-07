/**
 *  @file   larpandoracontent/LArReclustering/CheatedThreeDClusteringAlgorithm.cc
 *
 *  @brief  Implementation file for the reclustering algorithm that uses transverse calorimetric profiles.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArReclustering/CheatedThreeDClusteringAlgorithm.h"


using namespace pandora;

namespace lar_content
{

CheatedThreeDClusteringAlgorithm::CheatedThreeDClusteringAlgorithm()
{
}

StatusCode CheatedThreeDClusteringAlgorithm::Run()
{

    //Access the clusters
    //Split the clusters using truth info
    //Return the new cluster list

    //InitializeReclustering in has set the hits to be reclustered as current
    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    if (pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;


    ClusterVector clusterVector;
/*    if(m_drawProfiles)
    {
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));   
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), pCaloHitList, "CurrentClusterHits", BLUE));
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }*/

    
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {

        const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
        //if (TPC_VIEW_W != pParentCaloHit->GetHitType()) 
        //    continue; 
        //Find main contributing MC particle Id 
        int mainMcParticleIndex = this->GetMainMcParticleIndex(pParentCaloHit);


        //Loop over clusterVector entries
        const Cluster *pMainMCParticleCluster(nullptr);
        for (ClusterVector::const_iterator cluIter = clusterVector.begin(), cluIterEnd = clusterVector.end(); cluIter != cluIterEnd; ++cluIter)
        {
             const Cluster *const pCluster = *cluIter;
             //I only need to check the parent MC particle of one hit in the cluster. I'll take the front 
             CaloHitList clusterCaloHitList;;
             pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);
             const CaloHit *const pFrontHitParentCaloHit = static_cast<const CaloHit *>(clusterCaloHitList.front()->GetParentAddress());
             if(this->GetMainMcParticleIndex(pFrontHitParentCaloHit)==mainMcParticleIndex) pMainMCParticleCluster=pCluster;
             
        }

        //If pMainMCParticle is null (a cluster with matching MC particle hits is not found), create one
        if(!pMainMCParticleCluster)      
        {
            const Cluster *pCluster = nullptr;
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.push_back(pCaloHit);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
            clusterVector.push_back(pCluster);
        }
        else
        {
            //Attach calo hit to pMainMCParticleCluster
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pMainMCParticleCluster, pCaloHit));
        }

    }

    //Now I want to display all of these new clusters!
//    if(m_drawProfiles)PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));   
    for(const Cluster* pNewCluster: clusterVector)
    { 
      ClusterList newClusters;
      newClusters.push_back(pNewCluster);  
      //if(m_drawProfiles)PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &newClusters, "newClusters", AUTOITER));
    }
//    if(m_drawProfiles)PandoraMonitoringApi::ViewEvent(this->GetPandora());

    return STATUS_CODE_SUCCESS;
}

//Get vector of MC particles, sort them, loop into the vector and find in map (see "SortMCParticle")
int CheatedThreeDClusteringAlgorithm::GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit)
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


StatusCode CheatedThreeDClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DrawProfiles", m_drawProfiles));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
