package com.hartwig.hmftools.svannotation.analysis;

import java.util.List;

import com.google.common.collect.Lists;

public class SvClusteringMethods {

    public SvClusteringMethods()
    {

    }

    // unused for now..

    /*

    public void runClusteringByLocals() {
        // first find outer-most clusters at CRMS level and identify others which will span CRMS

        // LOGGER.debug("sample({}) isolating inter-chromosomal variants", mSampleId);

        // first remove cross-chromosomal from initial consideration
        // List<SvClusterData> crossChromosomal = findCrossChromosomalSVs();

        LOGGER.debug("sample({}) clustering contained variants", mSampleId);

        createOutermostChromosomalClusters();

        LOGGER.debug("sample({}) assigning sub-clusters", mSampleId);

        // now need to assign and create sub-clusters within these outer-most clusters
        for(SvCluster cluster : mClusters) {
            assignSubClusters(cluster);
        }

        LOGGER.debug("sample({}) matching {} overlapping variants", mSampleId, mUnassignedVariants.size());

        // final create links for the cross-chromosomal SVs and other overlapping SVs
        assignCrossClusterSVs();
    }

    private List<SvClusterData> findCrossChromosomalSVs()
    {
        // first remove cross-chromosomal from initial consideration
        List<SvClusterData> crossChromosomal = Lists.newArrayList();

        int currentIndex = 0;
        while (currentIndex < mUnassignedVariants.size()) {
            SvClusterData currentVar = mUnassignedVariants.get(currentIndex);

            // first skip over cross-CRMS for now
            if (!currentVar.isLocal()) {
                crossChromosomal.add(currentVar);

                LOGGER.debug("cross-CRMS SV: {}", currentVar.posId());

                mUnassignedVariants.remove(currentIndex); // index will point at next
                continue;
            } else {
                ++currentIndex;
            }
        }

        return crossChromosomal;

    }

    private void createOutermostChromosomalClusters()
    {
        int currentIndex = 0;

        List<SvClusterData> subSVsToRemove = Lists.newArrayList();

        while(currentIndex < mUnassignedVariants.size()) {
            SvClusterData currentVar = mUnassignedVariants.get(currentIndex);

            if(!currentVar.isLocal())
            {
                LOGGER.debug("skipping inter-chromosomal SV: {}", currentVar.posId());
                ++currentIndex;
                continue;
            }

            boolean spansOtherSVs = false;
            List<SvClusterData> subSVs = Lists.newArrayList();

            for (SvClusterData otherVar : mUnassignedVariants)
            {
                if(otherVar.id().equals(currentVar.id()))
                {
                    continue;
                }

                if(mClusteringUtils.isLocalOverlap(currentVar, otherVar))
                {
                    LOGGER.debug("local overlap SVs: v1({}) and v2({})", currentVar.posId(), otherVar.posId());
                    spansOtherSVs = true;
                }

                if (mClusteringUtils.isWithin(currentVar, otherVar))
                {
                    if(currentVar.addSubSV(otherVar)) {
                        LOGGER.debug("wholy contained: outer({}) and inner({})", currentVar.posId(), otherVar.posId());
                    }

                }
                else if (mClusteringUtils.isWithin(otherVar, currentVar))
                {
                    if(otherVar.addSubSV(currentVar)) {
                        LOGGER.debug("wholy contained: outer({}) and inner({})", otherVar.posId(), currentVar.posId());
                    }
                }
            }

            if(!currentVar.isSubSV() && (currentVar.hasSubSVs() || !spansOtherSVs))
            {
                LOGGER.debug("adding outer-most CRMS SV: {}", currentVar.posId());

                // make a new cluster
                SvCluster newCluster = new SvCluster(getNextClusterId(), mClusteringUtils);

                newCluster.setSpanningSV(currentVar);
                mClusters.add(newCluster);

                mUnassignedVariants.remove(currentIndex); // index will point at next
            }
            else
            {
                if(currentVar.isSubSV())
                    subSVsToRemove.add(currentVar);

                ++currentIndex;
            }
        }

        for(SvClusterData variant : subSVsToRemove)
        {
            // remove from unassigned
            mUnassignedVariants.remove(variant);
        }
    }

    private void assignSubClusters(SvCluster cluster)
    {
        LOGGER.debug("cluster({}) sub-clustering {} variants", cluster.getClusterId(), cluster.getSVs().size());

        // create sub-clusters for any SVs which have sub SVs
        for(SvClusterData currentVar : cluster.getSVs())
        {
            if(currentVar == cluster.getSpanningSV())
                continue;

            if(!currentVar.hasSubSVs())
            {
                if(currentVar.areClustersSet())
                {
                    // keep with the smallest cluster
                    long clusterSpan = currentVar.getStartCluster().getSpan();
                    long newSpan = cluster.getSpan();

                    if(newSpan == -1 || newSpan > clusterSpan)
                    {
                        continue;
                    }

                }

                LOGGER.debug("variant({}) assigned to cluster({}) span({} -> {})",
                        currentVar.posId(), cluster.getClusterId(), cluster.getSpanningSV().position(true), cluster.getSpanningSV().position(false));

                currentVar.setStartCluster(cluster);
                currentVar.setEndCluster(cluster);
                continue;
            }
            else {
                // also check that this variant is a sub-variant of another variant at this level
                boolean skipClustering = false;
                for(SvClusterData other : cluster.getSVs())
                {
                    if(currentVar.equals(other) || other == cluster.getSpanningSV())
                        continue;

                    if(other.getSubSVs().contains(currentVar))
                    {
                        LOGGER.debug("skipping sub-cluster for variant({}) since part of other({})", currentVar.id(), other.id());
                        skipClustering = true;
                        break;
                    }
                }

                if(skipClustering)
                    continue;
            }

            // create a cluster, set its parent cluster and then call iteratively
            SvCluster newCluster = new SvCluster(getNextClusterId(), mClusteringUtils);

            LOGGER.debug("cluster({}) creating new subCluster({})", cluster.getClusterId(), newCluster.getClusterId());

            // add the spanning SV and its sub-SVs
            newCluster.setSpanningSV(currentVar);
            cluster.addSubCluster(newCluster);

            // and call recursively to create lower-level clusters
            assignSubClusters(newCluster);
        }
    }

    private void assignCrossClusterSVs()
    {
        // for all remaining SVs, including those which cross chromosomes,
        // assign them to the most precise cluster (ie that with the smallest positional range)
        for(SvClusterData currentVar : mUnassignedVariants) {

            for(int i = 0; i < 2; ++i) {

                boolean isStart = (i == 0);
                String isStartStr = isStart ? "start" : "end";

                // first start position, the end
                boolean found = false;
                for (SvCluster cluster : mClusters) {
                    SvCluster matchedCluster = cluster.findClosestCluster(currentVar.chromosome(isStart), currentVar.position(isStart));
                    if (matchedCluster != null) {
                        LOGGER.debug("variant: id({}) {}({}:{}) matched with cluster({}) on spanningVariant({})",
                                currentVar.id(), isStartStr, currentVar.chromosome(isStart), currentVar.position(isStart),
                                matchedCluster.getClusterId(), matchedCluster.getSpanningSV().posId());

                        if(isStart)
                            currentVar.setStartCluster(matchedCluster);
                        else
                            currentVar.setEndCluster(matchedCluster);

                        found = true;
                        break;
                    }
                }

                if (!found) {
                    LOGGER.debug("variant: id({}) {}({}:{}) unmatched",
                            currentVar.id(), isStartStr, currentVar.chromosome(isStart), currentVar.position(isStart));
                }
            }
        }
    }

    */
}
