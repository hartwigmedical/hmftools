package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.PERMITED_DUP_BE_DISTANCE;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.haveLinkedAssemblies;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_DIFF;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_INFER_ONLY;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_DB;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ClusterAnalyser {

    final SvClusteringConfig mConfig;
    final SvUtilities mUtils;
    SvClusteringMethods mClusteringMethods;

    private ChainFinder mChainFinder;
    private LinkFinder mLinkFinder;

    public static int MIN_TEMPLATED_INSERTION_LENGTH = 30;
    private static int MAX_TEMPLATED_INSERTION_LENGTH = 500;
    public static double MAX_COPY_NUMBER_DIFF = 0.5;
    public static double MAX_COPY_NUMBER_DIFF_PERC = 0.1;
    public static int CLUSTER_SIZE_ANALYSIS_LIMIT = 200;
    public static int SMALL_CLUSTER_SIZE = 3;

    public static String TRANS_TYPE_TRANS = "TRANS";
    public static String TRANS_TYPE_SPAN = "SPAN";

    private static final Logger LOGGER = LogManager.getLogger(ClusterAnalyser.class);

    public ClusterAnalyser(final SvClusteringConfig config, final SvUtilities utils, SvClusteringMethods clusteringMethods)
    {
        mConfig = config;
        mUtils = utils;
        mClusteringMethods = clusteringMethods;
        mLinkFinder = new LinkFinder(mConfig, mUtils, mClusteringMethods);
        mChainFinder = new ChainFinder(mUtils);
        mChainFinder.setLogVerbose(mConfig.LogVerbose);
    }

    public void findSimpleCompleteChains(final String sampleId, List<SvCluster> clusters)
    {
        // for small clusters, try to find a full chain through all SVs
        for(SvCluster cluster : clusters)
        {
            if(cluster.getCount() == 1 || cluster.getCount() > SMALL_CLUSTER_SIZE)
                continue;

            if(!cluster.isConsistent())
                continue;

            // first establish links between SVs (eg TIs and DBs)
            mLinkFinder.findLinkedPairs(sampleId, cluster);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(sampleId, cluster);

            if(cluster.isFullyChained())
            {
                LOGGER.debug("sample({}) cluster({}) fully chained with {} SVs", sampleId, cluster.getId(), cluster.getCount());
            }
        }
    }

    public void findLinksAndChains(final String sampleId, List<SvCluster> clusters)
    {
        for(SvCluster cluster : clusters)
        {
            applyCopyNumberReplication(sampleId, cluster);

            if(cluster.getUniqueSvCount() < 2)
                continue;

            // first establish links between SVs (eg TIs and DBs)
            mLinkFinder.findLinkedPairs(sampleId, cluster);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(sampleId, cluster);
        }

        // now look at merging unresolved & inconsistent clusters where they share the same chromosomal arms
        mergeInconsistentClusters(clusters);

        for(SvCluster cluster : clusters)
        {
            if(cluster.hasSubClusters())
            {
                cluster.setDesc(cluster.getClusterTypesAsString());

                // repeat the search for inferred links now that additional SVs have been merged in
                List<SvLinkedPair> newLinkedPairs = mLinkFinder.createInferredLinkedPairs(sampleId, cluster, true);
                cluster.getInferredLinkedPairs().addAll(newLinkedPairs);

                createCopyNumberSegments(sampleId, cluster);

                findChains(sampleId, cluster);
            }

            mLinkFinder.resolveTransitiveSVs(sampleId, cluster);
            cacheFinalLinkedPairs(sampleId, cluster);
        }

        // any clusters which were merged to resolve a collection of them, but
        // which did not lead to any longer chains, are now de-merged
        demergeClusters(clusters);
    }

    private void findChains(final String sampleId, SvCluster cluster)
    {
        mChainFinder.initialise(sampleId, cluster);
        // mChainFinder.setRequireFullChains(true);

        boolean hasFullChain = mChainFinder.formClusterChains();

        if(!hasFullChain)
            return;

        cluster.setIsFullyChained(true);

        // remove any inferred link which isn't in the full chain
        final SvChain fullChain = cluster.getChains().get(0);

        List<SvLinkedPair> inferredLinkedPairs = cluster.getInferredLinkedPairs();
        inferredLinkedPairs.clear();

        for(SvLinkedPair pair : fullChain.getLinkedPairs())
        {
            if(pair.isInferred())
                inferredLinkedPairs.add(pair);
        }
    }

    private void applyCopyNumberReplication(final String sampleId, SvCluster cluster)
    {
        // first remove any replication previously performed before clusters were merged
        cluster.removeReplicatedSvs();
        mClusteringMethods.applyCopyNumberReplication(sampleId, cluster);
    }

    private void applyCopyNumberReplicationOld(final String sampleId, SvCluster cluster)
    {
        List<SvClusterData> newSVs = Lists.newArrayList();

        for(final SvClusterData var : cluster.getSVs())
        {
            int calcCopyNumber = var.impliedCopyNumber(true);

            int chromosomeCopyNumber = mClusteringMethods.getChrCopyNumberMap().get(var.chromosome(true));

            chromosomeCopyNumber /= 2;

            if(chromosomeCopyNumber >= 2)
                calcCopyNumber /= chromosomeCopyNumber;

            if(calcCopyNumber <= 1)
                continue;

            LOGGER.debug("sample({}) replicating SV({}) {} times", sampleId, var.posId(), calcCopyNumber-1);

            for(int i = 1; i < calcCopyNumber; ++i)
            {
                SvClusterData newVar = new SvClusterData(var);
                newSVs.add(newVar);
            }
        }

        for(SvClusterData var : newSVs)
        {
            cluster.addVariant(var);
        }
    }

    private void mergeInconsistentClusters(List<SvCluster> clusters)
    {
        // it's possible that to resolve arms and more complex arrangements, that clusters not merged
        // by proximity of overlaps must be put together to solve inconsistencies (ie loose ends)

        // merge any cluster which is itself not consistent and has breakends on the same arm as another
        int index1 = 0;
        while(index1 < clusters.size())
        {
            SvCluster cluster1 = clusters.get(index1);

            if(cluster1.isConsistent())
            {
                ++index1;
                continue;
            }

            if(!cluster1.isFullyChained() && cluster1.getUniqueSvCount() > 1 && !cluster1.hasSubClusters())
            {
                ++index1;
                continue;
            }

            boolean cluster1Merged = false;
            SvCluster newCluster = null;

            int index2 = index1 + 1;
            while(index2 < clusters.size())
            {
                SvCluster cluster2 = clusters.get(index2);

                if(cluster2.isConsistent())
                {
                    ++index2;
                    continue;
                }

                if(!cluster2.isFullyChained() && cluster2.getUniqueSvCount() > 1 && !cluster1.hasSubClusters())
                {
                    ++index2;
                    continue;
                }

                boolean cluster2Merged = false;

                boolean foundConnection = false;

                for(int be1 = SVI_START; be1 <= SVI_END; ++be1)
                {
                    final String c1LinkingArm = cluster1.linkingArm(isStart(be1));

                    if(c1LinkingArm.isEmpty())
                        continue;

                    final String c1LinkingChr = cluster1.linkingChromosome(isStart(be1));

                    for(int be2 = SVI_START; be2 <= SVI_END; ++be2)
                    {
                        final String c2LinkingArm = cluster2.linkingArm(isStart(be2));

                        if(c2LinkingArm.isEmpty())
                            continue;

                        if(c1LinkingArm.equals(c2LinkingArm) && c1LinkingChr.equals(cluster2.linkingChromosome(isStart(be2))))
                        {
                            LOGGER.debug("inconsistent cluster({}) and cluster({}) linked on chrArm({}:{})",
                                    cluster1.getId(), cluster2.getId(), c1LinkingChr, c1LinkingArm);

                            foundConnection = true;
                            break;
                        }
                    }

                    if(foundConnection)
                        break;
                }

                if(!foundConnection)
                {
                    ++index2;
                    continue;
                }

                if(cluster1.hasSubClusters())
                {
                    LOGGER.debug("cluster({} svs={}) merges in cluster({} svs={})",
                            cluster1.getId(), cluster1.getCount(), cluster2.getId(), cluster2.getCount());

                    cluster2Merged = true;
                    cluster1.addSubCluster(cluster2);
                }
                else
                {
                    cluster1Merged = true;
                    cluster2Merged = true;

                    newCluster = new SvCluster(mClusteringMethods.getNextClusterId(), mUtils);

                    LOGGER.debug("new cluster({}) from merge of cluster({} svs={}) and cluster({} svs={})",
                            newCluster.getId(), cluster1.getId(), cluster1.getCount(), cluster2.getId(), cluster2.getCount());

                    newCluster.addSubCluster(cluster1);
                    newCluster.addSubCluster(cluster2);
                }

                if(cluster2Merged)
                    clusters.remove(index2);
                else
                    ++index2;

                if(cluster1Merged)
                    break;
            }

            if(cluster1Merged && newCluster != null)
            {
                // cluster has been replaced with a commbined one
                clusters.remove(index1);
                clusters.add(index1, newCluster);
            }
            else
            {
                ++index1;
            }
        }
    }

    private void demergeClusters(List<SvCluster> clusters)
    {
        // de-merge any clusters which didn't form longer chains
        int clusterCount = clusters.size();
        for(int i = 0; i < clusterCount;)
        {
            SvCluster cluster = clusters.get(i);

            if (!cluster.hasSubClusters())
            {
                ++i;
                continue;
            }

            int mainChainCount = cluster.getMaxChainCount();
            int maxSubClusterChainCount = 0;

            for(final SvCluster subCluster : cluster.getSubClusters())
            {
                maxSubClusterChainCount = max(maxSubClusterChainCount, subCluster.getMaxChainCount());
            }

            if(mainChainCount > maxSubClusterChainCount)
            {
                ++i;
                continue;
            }

            for(final SvCluster subCluster : cluster.getSubClusters())
            {
                clusters.add(subCluster);
            }

            if(mainChainCount > 0)
            {
                LOGGER.debug("removed cluster({}) since maxChainCount({}) less than subclusters({}) maxSubClusterChainCount({})",
                        cluster.getId(), mainChainCount, cluster.getSubClusters().size(), maxSubClusterChainCount);
            }

            clusters.remove(i);
        }
    }

    private void createCopyNumberSegments(final String sampleId, final SvCluster cluster)
    {
        int cnId = 0;
        final List<SvClusterData> clusterSVs = cluster.getSVs();

        Map<String, List<SvCNData>> chrCNDataMap = new HashMap();

        for(Map.Entry<String, List<SvBreakend>> entry : mClusteringMethods.getChrBreakendMap().entrySet())
        {
            final String chromosome = entry.getKey();

            List<SvBreakend> breakendList = entry.getValue();

            List<SvCNData> copyNumberData = Lists.newArrayList();

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvClusterData var = breakend.getSV();

                double copyNumber = var.copyNumber(breakend.usesStart());
                double copyNumberChange = var.copyNumberChange(breakend.usesStart());
                double prevCopyNumber = copyNumber - copyNumberChange;

                int adjCopyNumber = (int)round(copyNumber/2 - 1);

                //LOGGER.debug(String.format("sample(%s) chr(%s) seg %d: copyNumber(%d %.2f prev=%.2f chg=%.2f)",
                //        sampleId, chromosome, i, adjCopyNumber, copyNumber, prevCopyNumber, copyNumberChange));

                if(clusterSVs.contains(var))
                {
                    SvCNData cnData = new SvCNData(cnId++, chromosome, breakend.position(), 0, "", "", 0, 0, 0, adjCopyNumber, "");
                    copyNumberData.add(cnData);
                }
            }

            chrCNDataMap.put(chromosome, copyNumberData);
        }

        cluster.setChrCNData(chrCNDataMap);
    }

    private void cacheFinalLinkedPairs(final String sampleId, SvCluster cluster)
    {
        List<SvLinkedPair> linkedPairs;

        if(cluster.isFullyChained())
        {
            linkedPairs = cluster.getChains().get(0).getLinkedPairs();
        }
        else
        {
            linkedPairs = Lists.newArrayList();

            // add all chained links
            for(final SvChain chain : cluster.getChains())
            {
                linkedPairs.addAll(chain.getLinkedPairs());
            }

            // any any unchained assembly links
            for(final SvLinkedPair pair : cluster.getAssemblyLinkedPairs())
            {
                if(!linkedPairs.contains(pair))
                    linkedPairs.add(pair);
            }

            // finally add any other potential inferred links which don't clash with existing links
            for(final SvLinkedPair pair : cluster.getInferredLinkedPairs())
            {
                if(linkedPairs.contains(pair))
                    continue;

                boolean hasClash = false;
                for(final SvLinkedPair existingLink : linkedPairs)
                {
                    if (existingLink.hasLinkClash(pair))
                    {
                        hasClash = true;
                        break;
                    }
                }

                if(!hasClash)
                    linkedPairs.add(pair);
            }

            // mark the resultant set of inferred links
            for (SvLinkedPair pair : linkedPairs)
            {
                if(pair.isInferred())
                {
                    pair.first().setAssemblyMatchType(ASSEMBLY_MATCH_INFER_ONLY, pair.firstLinkOnStart());
                    pair.second().setAssemblyMatchType(ASSEMBLY_MATCH_INFER_ONLY, pair.secondLinkOnStart());
                }
            }
        }

        if(!linkedPairs.isEmpty())
        {
            LOGGER.info("sample({}) cluster({}: {} count={}) has {} linked pairs",
                    sampleId, cluster.getId(), cluster.getDesc(), cluster.getUniqueSvCount(), linkedPairs.size());

            if(LOGGER.isDebugEnabled())
            {
                for (final SvLinkedPair pair : linkedPairs)
                {
                    LOGGER.debug("linked {} {} pair length({}) variants({})",
                            pair.isInferred() ? "inferred" : "assembly", pair.linkType(), pair.length(), pair.toString());
                }
            }
        }

        cluster.setLinkedPairs(linkedPairs);
    }

    public static void reduceInferredToShortestLinks(List<SvLinkedPair> linkedPairs, List<SvChain> chains)
    {
        // any linked pair used in a chain must be kept
        List<SvLinkedPair> reqLinkedPairs = Lists.newArrayList();
        for(SvChain chain : chains)
        {
            reqLinkedPairs.addAll(chain.getLinkedPairs());
        }

        // now remove mutually exclusive linked sections by using the shortest first (they are already ordered)
        int i = 0;
        while(i < linkedPairs.size())
        {
            final SvLinkedPair pair = linkedPairs.get(i);

            if(reqLinkedPairs.contains(pair))
            {
                ++i;
                continue;
            }

            boolean removeFirst = false;

            // search for another pair with a matching breakend
            for(int j = i+1; j < linkedPairs.size();)
            {
                final SvLinkedPair pair2 = linkedPairs.get(j);

                if((pair.first().equals(pair2.first()) && pair.firstLinkOnStart() == pair2.firstLinkOnStart())
                || (pair.first().equals(pair2.second()) && pair.firstLinkOnStart() == pair2.secondLinkOnStart())
                || (pair.second().equals(pair2.first()) && pair.secondLinkOnStart() == pair2.firstLinkOnStart())
                || (pair.second().equals(pair2.second()) && pair.secondLinkOnStart() == pair2.secondLinkOnStart()))
                {
                    // to avoid logging unlikely long TIs
                    //LOGGER.debug("removing duplicate linked pair({} len={}) vs shorter({} len={})",
                    //        pair2.toString(), pair2.length(), pair.toString(), pair.length());

                    if(reqLinkedPairs.contains(pair2))
                    {
                        removeFirst = true;
                        break;
                    }

                    // remove the longer pair with a duplicate breakend
                    linkedPairs.remove(j);
                }
                else
                {
                    ++j;
                }
            }

            if(removeFirst)
            {
                // LOGGER.debug("duplicate linked pair({} len={})", pair.toString(), pair.length(), pair.first().getTransSvLinks());
                linkedPairs.remove(i);
            }
            else
            {
                ++i;
            }
        }
    }

    public void logStats()
    {
        mChainFinder.getContinuousFinderPc().logStats(false);
        mChainFinder.getRecursiveFinderPc().logStats(false);
    }

}
