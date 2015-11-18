/**
 *  @file   LArContent/src/LArVertex/VertexSelectionAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex selection algorithm class.
 * 
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArUtility/KDTreeLinkerAlgoT.h"

#include "LArVertex/VertexSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

VertexSelectionAlgorithm::VertexSelectionAlgorithm() :
    m_replaceCurrentVertexList(true),
    m_beamMode(false),
    m_selectSingleVertex(true),
    m_kernelEstimateSigma(0.026f),
    m_maxOnHitDisplacement(1.f),
    m_maxHitVertexDisplacement1D(100.f),
    m_maxTopScoreCandidates(50),
    m_maxTopScoreSelections(3),
    m_maxBeamTopScoreCandidates(50),
    m_maxBeamTopScoreSelections(3),
    m_minCandidateDisplacement(2.f),
    m_minCandidateScoreFraction(0.5f),
    m_minBeamCandidateScoreFraction(0.5f),
    m_nDecayLengthsInZSpan(2.f),
    m_bestScoreMultiplier(0.75f),
    m_bestBeamScoreMultiplier(1.1f),
    m_mustUseBeamScoreMultiplier(1.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::Run()
{
    const VertexList *pInputVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pInputVertexList));

    HitKDTree2D kdTreeU, kdTreeV, kdTreeW;
    this->InitializeKDTrees(kdTreeU, kdTreeV, kdTreeW);

    VertexScoreList vertexScoreList;
    float minZCoordinate(std::numeric_limits<float>::max()), maxZCoordinate(-std::numeric_limits<float>::max());

    for (const Vertex *const pVertex : *pInputVertexList)
    {
        try
        {
            float figureOfMerit(this->GetFigureOfMerit(pVertex, kdTreeU, kdTreeV, kdTreeW));
            vertexScoreList.push_back(VertexScore(pVertex, figureOfMerit));

            if (pVertex->GetPosition().GetZ() < minZCoordinate)
                minZCoordinate = pVertex->GetPosition().GetZ();

            if (pVertex->GetPosition().GetZ() > maxZCoordinate)
                maxZCoordinate = pVertex->GetPosition().GetZ();
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }

    std::sort(vertexScoreList.begin(), vertexScoreList.end());
    const float zSpan(maxZCoordinate - minZCoordinate);
    const float decayConstant((zSpan < std::numeric_limits<float>::epsilon()) ? 0.f : (m_nDecayLengthsInZSpan / zSpan));

    VertexScoreList selectedVertexScoreList;
    this->SelectTopScoreVertices(vertexScoreList, selectedVertexScoreList);
    this->SelectTopScoreBeamVertices(vertexScoreList, minZCoordinate, decayConstant, selectedVertexScoreList);

    VertexList finalVertexList;
    this->SelectFinalVertices(selectedVertexScoreList, minZCoordinate, decayConstant, finalVertexList);

    if (!finalVertexList.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_outputVertexListName, finalVertexList));

        if (m_replaceCurrentVertexList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::InitializeKDTrees(HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const
{
    this->InitializeKDTree(m_inputCaloHitListNameU, kdTreeU);
    this->InitializeKDTree(m_inputCaloHitListNameV, kdTreeV);
    this->InitializeKDTree(m_inputCaloHitListNameW, kdTreeW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::InitializeKDTree(const std::string &caloHitListName, HitKDTree2D &kdTree) const
{
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, caloHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexSelectionAlgorithm: unable to find calo hit list " << caloHitListName << std::endl;

        return;
    }

    HitKDNode2DList hitKDNode2DList;
    KDTreeBox hitsBoundingRegion2D = fill_and_bound_2d_kd_tree(this, *pCaloHitList, hitKDNode2DList, true);
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::GetFigureOfMerit(const Vertex *const pVertex, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const
{
    KernelEstimate kernelEstimateU(m_kernelEstimateSigma);
    KernelEstimate kernelEstimateV(m_kernelEstimateSigma);
    KernelEstimate kernelEstimateW(m_kernelEstimateSigma);

    if ((!kdTreeU.empty() && !this->IsVertexOnHit(pVertex, TPC_VIEW_U, kdTreeU)) ||
        (!kdTreeV.empty() && !this->IsVertexOnHit(pVertex, TPC_VIEW_V, kdTreeV)) ||
        (!kdTreeW.empty() && !this->IsVertexOnHit(pVertex, TPC_VIEW_W, kdTreeW)))
    {
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
    }

    this->FillKernelEstimate(pVertex, TPC_VIEW_U, kdTreeU, kernelEstimateU);
    this->FillKernelEstimate(pVertex, TPC_VIEW_V, kdTreeV, kernelEstimateV);
    this->FillKernelEstimate(pVertex, TPC_VIEW_W, kdTreeW, kernelEstimateW);

    return this->GetFigureOfMerit(kernelEstimateU, kernelEstimateV, kernelEstimateW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::GetFigureOfMerit(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
    const KernelEstimate &kernelEstimateW) const
{
    float figureOfMerit(0.f);
    figureOfMerit += this->GetFigureOfMerit(kernelEstimateU);
    figureOfMerit += this->GetFigureOfMerit(kernelEstimateV);
    figureOfMerit += this->GetFigureOfMerit(kernelEstimateW);

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::GetFigureOfMerit(const KernelEstimate &kernelEstimate) const
{
    float figureOfMerit(0.f);
    const KernelEstimate::ContributionList &contributionList(kernelEstimate.GetContributionList());

    for (const KernelEstimate::ContributionList::value_type &contribution : contributionList)
    {
        figureOfMerit += contribution.second * kernelEstimate.Sample(contribution.first);
    }

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::IsVertexOnHit(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

    CaloHitList nearbyHits;
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxOnHitDisplacement, m_maxOnHitDisplacement);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    return (!found.empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::FillKernelEstimate(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree, KernelEstimate &kernelEstimate) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

    CaloHitList nearbyHits;
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxHitVertexDisplacement1D, m_maxHitVertexDisplacement1D);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    for (const auto &hit : found)
    {
        const CartesianVector displacement(hit.data->GetPositionVector() - vertexPosition2D);
        const float magnitude(displacement.GetMagnitude());

        if (magnitude < std::numeric_limits<float>::epsilon())
            continue;

        const float phi(this->atan2Fast(displacement.GetZ(), displacement.GetX()));
        const float weight(1.f / std::sqrt(magnitude));
        kernelEstimate.AddContribution(phi, weight);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::SelectTopScoreVertices(const VertexScoreList &vertexScoreList, VertexScoreList &selectedVertexScoreList) const
{
    // ATTN Assumes sorted vertex score list
    unsigned int nVerticesConsidered(0);

    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        if (++nVerticesConsidered > m_maxTopScoreCandidates)
            break;

        if (selectedVertexScoreList.size() >= m_maxTopScoreSelections)
            break;

        if (!selectedVertexScoreList.empty() && !this->AcceptVertexLocation(iter->GetVertex(), selectedVertexScoreList))
            continue;

        if (!selectedVertexScoreList.empty() && !this->AcceptVertexScore(iter->GetScore(), m_minCandidateScoreFraction, selectedVertexScoreList))
            continue;

        selectedVertexScoreList.push_back(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::SelectTopScoreBeamVertices(const VertexScoreList &vertexScoreList, const float minZCoordinate,
    const float decayConstant, VertexScoreList &selectedVertexScoreList) const
{
    if (!m_beamMode)
        return;

    VertexScoreList beamVertexScoreList;

    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        const float beamScore(iter->GetScore() * std::exp(-(iter->GetVertex()->GetPosition().GetZ() - minZCoordinate) * decayConstant));
        beamVertexScoreList.push_back(VertexScore(iter->GetVertex(), beamScore));
    }

    unsigned int nVerticesConsidered(0), nVerticesAdded(0);
    std::sort(beamVertexScoreList.begin(), beamVertexScoreList.end());

    for (VertexScoreList::const_iterator iter = beamVertexScoreList.begin(), iterEnd = beamVertexScoreList.end(); iter != iterEnd; ++iter)
    {
        if (++nVerticesConsidered > m_maxBeamTopScoreCandidates)
            break;

        if (nVerticesAdded >= m_maxBeamTopScoreSelections)
            break;

        if (!selectedVertexScoreList.empty() && !this->AcceptVertexLocation(iter->GetVertex(), selectedVertexScoreList))
            continue;

        const float nonBeamScore(iter->GetScore() / std::exp(-(iter->GetVertex()->GetPosition().GetZ() - minZCoordinate) * decayConstant));

        if (!selectedVertexScoreList.empty() && !this->AcceptVertexScore(nonBeamScore, m_minBeamCandidateScoreFraction, selectedVertexScoreList))
            continue;

        ++nVerticesAdded;
        selectedVertexScoreList.push_back(VertexScore(iter->GetVertex(), nonBeamScore));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::AcceptVertexLocation(const Vertex *const pVertex, const VertexScoreList &vertexScoreList) const
{
    const CartesianVector position(pVertex->GetPosition());

    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        if (iter->GetVertex() == pVertex)
            return false;

        const float displacement3D((position - iter->GetVertex()->GetPosition()).GetMagnitude());

        if (displacement3D < m_minCandidateDisplacement)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::AcceptVertexScore(const float score, const float minScoreFraction, const VertexScoreList &vertexScoreList) const
{
    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        if (score < minScoreFraction * iter->GetScore())
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::SelectFinalVertices(const VertexScoreList &vertexScoreList, const float minZCoordinate, const float decayConstant,
    VertexList &finalVertexList) const
{
    const Vertex *pFinalVertex(NULL);
    float bestScore(0.f), bestBeamScore(0.f);

    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        if (!m_selectSingleVertex)
        {
            finalVertexList.insert(iter->GetVertex());
            continue;
        }

        if (!m_beamMode && (iter->GetScore() < bestScore))
            continue;

        const float beamScore(iter->GetScore() * std::exp(-(iter->GetVertex()->GetPosition().GetZ() - minZCoordinate) * decayConstant));

        if (m_beamMode && (beamScore < m_bestBeamScoreMultiplier * bestBeamScore))
            continue;

        if (m_beamMode && (iter->GetScore() < m_bestScoreMultiplier * bestScore) && (beamScore < m_mustUseBeamScoreMultiplier * bestBeamScore))
            continue;

        pFinalVertex = iter->GetVertex();
        bestScore = iter->GetScore();
        bestBeamScore = beamScore;
    }

    if (m_selectSingleVertex && (NULL != pFinalVertex))
        finalVertexList.insert(pFinalVertex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::atan2Fast(const float y, const float x) const
{
    const float ONE_QTR_PI(0.25f * M_PI);
    const float THR_QTR_PI(0.75f * M_PI);

    const float abs_y(std::max(std::fabs(y), std::numeric_limits<float>::epsilon()));
    const float abs_x(std::fabs(x));

    const float r((x < 0.f) ? (x + abs_y) / (abs_y + abs_x) : (abs_x - abs_y) / (abs_x + abs_y));
    const float angle(((x < 0.f) ? THR_QTR_PI : ONE_QTR_PI) + (0.1963f * r * r - 0.9817f) * r);

    return ((y < 0.f) ? -angle : angle); // negate if in quad III or IV
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::KernelEstimate::Sample(const float x) const
{
    const float bandwidth(this->GetBandwidth());
    const float weightSum(this->GetWeightSum());

    if ((bandwidth < std::numeric_limits<float>::epsilon()) || (weightSum < std::numeric_limits<float>::epsilon()))
        return 0.f;

    const ContributionList &contributionList(this->GetContributionList());
    ContributionList::const_iterator lowerIter(contributionList.lower_bound(x - 3.f * bandwidth));
    ContributionList::const_iterator upperIter(contributionList.upper_bound(x + 3.f * bandwidth));

    float sample(0.f);
    const float gaussConstant(1.f / std::sqrt(2.f * M_PI * bandwidth * bandwidth));

    for (ContributionList::const_iterator iter = lowerIter; iter != upperIter; ++iter)
    {
        const float deltaSigma((x - iter->first) / bandwidth);
        const float gaussian(gaussConstant * std::exp(-0.5f * deltaSigma * deltaSigma));
        sample += iter->second * gaussian;
    }

    return (sample / weightSum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::KernelEstimate::GetBandwidth() const
{
    if (!m_bandWidth.IsInitialized())
        m_bandWidth = 1.06f * m_sigma * std::pow(m_weightSum, -0.2f);

    return m_bandWidth.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::KernelEstimate::AddContribution(const float x, const float weight)
{
    m_contributionList.insert(ContributionList::value_type(x, weight));
    m_weightSum += weight;
    m_bandWidth.Reset();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameU", m_inputCaloHitListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameV", m_inputCaloHitListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameW", m_inputCaloHitListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BeamMode", m_beamMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectSingleVertex", m_selectSingleVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KernelEstimateSigma", m_kernelEstimateSigma));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxOnHitDisplacement", m_maxOnHitDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitVertexDisplacement1D", m_maxHitVertexDisplacement1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTopScoreCandidates", m_maxTopScoreCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTopScoreSelections", m_maxTopScoreSelections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxBeamTopScoreCandidates", m_maxBeamTopScoreCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxBeamTopScoreSelections", m_maxBeamTopScoreSelections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateDisplacement", m_minCandidateDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateScoreFraction", m_minCandidateScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinBeamCandidateScoreFraction", m_minBeamCandidateScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NDecayLengthsInZSpan", m_nDecayLengthsInZSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BestScoreMultiplier", m_bestScoreMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BestBeamScoreMultiplier", m_bestBeamScoreMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MustUseBeamScoreMultiplier", m_mustUseBeamScoreMultiplier));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
