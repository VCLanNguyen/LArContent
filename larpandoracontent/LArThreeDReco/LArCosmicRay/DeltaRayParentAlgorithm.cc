/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayParentAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray parent algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayParentAlgorithm.h"

using namespace pandora;

namespace lar_content
{

DeltaRayParentAlgorithm::DeltaRayParentAlgorithm() :
    m_distanceForMatching(5.f)    
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayParentAlgorithm::Run()
{
    const PfoList *pMuonPfoList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_muonPfoListName, pMuonPfoList));

    if (!pMuonPfoList || pMuonPfoList->empty())
        return STATUS_CODE_SUCCESS;
    
    const PfoList *pDeltaRayPfoList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_deltaRayPfoListName, pDeltaRayPfoList));

    if (!pDeltaRayPfoList || pDeltaRayPfoList->empty())
        return STATUS_CODE_SUCCESS;

    PfoLengthMap pfoLengthMap;
    this->InitialisePfoLengthMap(pMuonPfoList, pDeltaRayPfoList, pfoLengthMap);
    
    PfoVector deltaRayPfoVector(pDeltaRayPfoList->begin(), pDeltaRayPfoList->end());
    std::sort(deltaRayPfoVector.begin(), deltaRayPfoVector.end(), LArPfoHelper::SortByNHits);
    
    for (const ParticleFlowObject *const pPfo : deltaRayPfoVector)
    {
        if (!pPfo->GetParentPfoList().empty())
            continue;

        const ParticleFlowObject *pParentPfo(nullptr);
        this->FindParentPfo(pfoLengthMap, pPfo, pParentPfo);

        if (!pParentPfo)
            continue;
        
        this->AssignToParentPfo(pMuonPfoList, pDeltaRayPfoList, pPfo, pParentPfo, pfoLengthMap);
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayParentAlgorithm::InitialisePfoLengthMap(const PfoList *const muonPfoList, const PfoList *const deltaRayPfoList, PfoLengthMap &pfoLengthMap) const
{
    for (const ParticleFlowObject *const pPfo : *muonPfoList)
        pfoLengthMap[pPfo] = LArPfoHelper::GetTwoDLengthSquared(pPfo);

    for (const ParticleFlowObject *const pPfo : *deltaRayPfoList)
        pfoLengthMap[pPfo] = LArPfoHelper::GetTwoDLengthSquared(pPfo);    
}    
    
//------------------------------------------------------------------------------------------------------------------------------------------    

void DeltaRayParentAlgorithm::FindParentPfo(const PfoLengthMap &pfoLengthMap, const ParticleFlowObject *const pPfo, const ParticleFlowObject *&pParentPfo) const
{
    PfoVector allPfoVector;

    for (auto &entry : pfoLengthMap) 
      allPfoVector.push_back(entry.first);

    const PfoLengthMap::const_iterator currentIter(pfoLengthMap.find(pPfo));

    if (currentIter == pfoLengthMap.end())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float lengthSquared(currentIter->second);    
    float bestDistance(m_distanceForMatching);    

    for (const ParticleFlowObject *const pTestParent : allPfoVector)
    {
        if (pTestParent == pPfo)
            continue;
        
        const PfoLengthMap::const_iterator testIter(pfoLengthMap.find(pTestParent));

        if (testIter == pfoLengthMap.end())
            throw StatusCodeException(STATUS_CODE_FAILURE);
        
        if (testIter->second < lengthSquared)
            continue;

        float distance(std::numeric_limits<float>::max());

	    if (LArPfoHelper::GetTwoDSeparation(pPfo, pTestParent, distance) == STATUS_CODE_NOT_FOUND)
            continue;

        if (distance < bestDistance)
        {
            pParentPfo = pTestParent;
            bestDistance = distance;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayParentAlgorithm::AssignToParentPfo(const PfoList *const pMuonPfoList, const PfoList *const pDeltaRayPfoList, const ParticleFlowObject *const pPfo,
    const ParticleFlowObject *const pParentPfo, PfoLengthMap &pfoLengthMap) const
{
    if (std::find(pMuonPfoList->begin(), pMuonPfoList->end(), pParentPfo) != pMuonPfoList->end())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pParentPfo, pPfo));

    if (std::find(pDeltaRayPfoList->begin(), pDeltaRayPfoList->end(), pParentPfo) != pDeltaRayPfoList->end())
    {
        ClusterList pfoClusters;
        LArPfoHelper::GetTwoDClusterList(pPfo, pfoClusters);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo, m_deltaRayPfoListName));

        for (const Cluster *const pCluster : pfoClusters)
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pParentPfo, pCluster));

        this->UpdatePfoLengthMap({pPfo, pParentPfo}, pParentPfo, pfoLengthMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayParentAlgorithm::UpdatePfoLengthMap(const PfoList &pfosToRemove, const ParticleFlowObject *const pPfoToAdd, PfoLengthMap &pfoLengthMap) const
{
    for (const ParticleFlowObject *const pPfoToRemove : pfosToRemove)
    {
        const PfoLengthMap::const_iterator iter(pfoLengthMap.find(pPfoToRemove));

        if(iter == pfoLengthMap.end())
            throw StatusCodeException(STATUS_CODE_FAILURE);
    
        pfoLengthMap.erase(iter);
    }

    pfoLengthMap[pPfoToAdd] = LArPfoHelper::GetTwoDLengthSquared(pPfoToAdd);
}    
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode DeltaRayParentAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MuonPfoListName", m_muonPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "DeltaRayPfoListName", m_deltaRayPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceForMatching", m_distanceForMatching));    
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
