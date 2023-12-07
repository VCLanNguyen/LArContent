/**
 *  @file   larpandoracontent/LArReclustering/CheatedThreeDClusteringAlgorithm.h
 *
 *  @brief  Header file for the reclustering algorithm that uses transverse calorimetric profiles.
 *
 *  $Log: $
 */

#ifndef LAR_CHEATED_THREE_D_CLUSTERING_ALGORITHM_H
#define LAR_CHEATED_THREE_D_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{
/**
 *  @brief  RecursivePfoMopUpAlgorithm class
 */
class CheatedThreeDClusteringAlgorithm : public pandora::Algorithm
{
public:

    CheatedThreeDClusteringAlgorithm();

private:

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    int GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit);

    std::string m_mcParticleListName; ///< The mc particle list name 
    bool m_drawProfiles; //Boolean to enable and disable displaying transverse profiles
};

} // namespace lar_content

#endif // #endif LAR_TRANSVERSE_CALO_RECLUSTERING_ALGORITHM_H
