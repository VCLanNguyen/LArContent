/**
 *  @file   larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h
 *
 *  @brief  Header file for the reclustering algorithm that runs other algs.
 *
 *  $Log: $
 */

#ifndef LAR_THREE_D_RECLUSTERING_ALGORITHM_H
#define LAR_THREE_D_RECLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{
/**
 *  @brief  RecursivePfoMopUpAlgorithm class
 */
class ThreeDReclusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDReclusteringAlgorithm();

   /**
    *  @brief  Destructor
    */
    ~ThreeDReclusteringAlgorithm();

private:

    class Canvas
    {
    public:
        /**
         *  @brief Default constructor
         */
        Canvas(const int width, const int height,
		const float xMin, const float xMax, const float yMin, const float yMax);

        virtual ~Canvas();

        float **m_canvas;
        bool **m_visited;
        const int m_width;
        const int m_height;
        const float m_xMin;
        const float m_xMax;
        const float m_yMin;
        const float m_yMax;
    };

    pandora::StatusCode Run();

    /**
     *  @brief 
     *
     *  @param 
     *
     *  @return List of 
     */
    pandora::StatusCode FindBestThreeDClusters(pandora::ClusterList clusterList3D, std::string &reclusterListName, pandora::ClusterList &minimumFigureOfMeritClusterList);

    /**
     *  @brief Get the theoretical lateral profile value at the shower maximum for a given energy and radius (Grindhammer)
     *
     *  @param energy the cluster energy in MeV
     *  @param radiusInCm the radius at which to calculate the profile value
     *
     *  @return List of PfoMergeStats for each Pfo
     */
    float GetLateralProfileAtShowerMaximum(float clusterEnergyInMeV, float radiusInCm);
    float GetLateralProfileAtShowerDepth(float clusterEnergyInMeV, float radiusInCm, float ratioMaxDepth);


    float GetFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D);
    
    float GetFigureOfMerit(std::string figureOfMeritName, pandora::CaloHitList mergedClusterCaloHitList3D);

    float GetFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D, std::vector<pandora::CaloHitList> newClustersCaloHitList3D);

    float GetFigureOfMerit(std::string figureOfMeritName, pandora::CaloHitList mergedClusterCaloHitList3D, std::vector<pandora::CaloHitList> newClustersCaloHitLists3D);

    /** 
     *  @brief Use this to find cheated FOM by finding the shower purity as the fraction of hits that the main Mc particle is contributing to
     * @param Merged cluster hit list
     * @return The figure of merit
     */
    float GetCheatedFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D);
    
    float GetCheatedFigureOfMerit(std::vector<pandora::CaloHitList> newClustersCaloHitList3D);



    /** 
     *  @brief Use this to find FOM by comparing observed transverse profiles to prediction for one shower
     * @param Merged cluster hit list
     * @return The figure of merit
     */
    float GetTransverseProfileFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D);

    float GetTransverseDepthRatioFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D);

    float GetLongitudinalProfileFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D);

    float GetLongitudinalProfileFigureOfMerit(std::vector<pandora::CaloHitList> newClustersCaloHitList3D);

    /** 
     *  @brief Use this to find FOM by comparing observed transverse profile (merged) to prediction produced as combination of a list of underlying clusters
     * @param Merged cluster hit list
     * @param Reclustered clusters hit lists vector

     * @return The figure of merit
     */
    float GetTransverseProfileFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D, std::vector<pandora::CaloHitList> newClustersCaloHitList3D); //Eventually, perhaps move each of these in its own class, and let choose which ones to calculate via XML, like the clustering algos themselves


    //Decide whether to try reclustering for this pfo
    bool PassesCutsForReclustering(const pandora::ParticleFlowObject *const pPfo);

    int GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit);
    int GetMCParticleHierarchyTier(const pandora::MCParticle *const pMCParticle);

    //This can then go in the SDK with the name/format PandoraContentApi::RebuildLArTPCPfo(Algorithm &, Pfo  *pPfoToRebuild, Cluster *pTemplate3DCluster, MapOfListNameInfo &)
    pandora::StatusCode RebuildPfo(const pandora::Pfo *pPfoToRebuild, pandora::ClusterList &newThreeDClusters); 

    // Get kMeans prediction for k = 2   
    pandora::StatusCode GetKmeansPredictions(std::vector<float> caloHitProjectedPosition1Vect, std::vector<float> caloHitProjectedPosition2Vect, std::vector<float> caloHitEnergyVect, pandora::CartesianVector &center1, pandora::CartesianVector &center2, float &energy1, float &energy2);

    //Grow Peak function
    /*
     *  @brief  Create a list of vertices from a canvas
     *
     *  @param  canvases The input canvases
     *  @param  positionVector The output vector of wire plane positions
     *
     *  @return The StatusCode resulting from the function
     **/
    pandora::StatusCode GetVerticesFromCanvas(const Canvas &canvas, pandora::CartesianPointVector &vertices, std::vector<float> &energies) const;
    
    /**
     *  @brief  Determine if the pixel under consideration is part of a peak and grow that peak to include all connected pixels of equal value
     *
     *  @param  canvas The canvas within which peaks are sought
     *  @param  col The column of the pixel under consideration
     *  @param  row The row of the pixel under consideration
     *  @param  intensity The target intensity of the candidate peak
     *  @param  peak The output vector of pixels constituting the peak under consideration
     *
     *  @return true if we found a better peak while growing the current region, false otherwise
     */
    bool GrowPeak(const Canvas &canvas, int col, int row, float intensity, std::vector<std::pair<int, int>> &peak) const; 


    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_pfoListName; ///< The input pfo list name (e.g. list of neutrino or testbeam pfos)
    pandora::StringVector m_figureOfMeritNames; ///what figure(s) of merit to use
    std::string m_newPfosListNameAllAfterReclustering;
    bool m_visualDisplaysOn; //Boolean to enable and disable displaying transverse profiles
    int m_hitThresholdForNewPfo; //Minimum nr. of hits to form new 3Dcluster and pfo
    pandora::StringVector   m_clusteringAlgorithms; ///< The ordered list of clustering algorithms to be used
    std::string m_mcParticleListName; ///< The mc particle list name

    float m_peakFinderThreshold; //Set threshold for peak intensity

    bool m_writeToTree;
    std::string m_treeName; 
    std::string m_fileName;
    double m_eventId;
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_RECLUSTERING_ALGORITHM_H
