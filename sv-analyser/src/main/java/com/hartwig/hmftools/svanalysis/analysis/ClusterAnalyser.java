package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.svanalysis.analysis.SvCluster.findCluster;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.findVariantById;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.haveLinkedAssemblies;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_INFER_ONLY;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ClusterAnalyser {

    final SvClusteringConfig mConfig;
    final SvUtilities mUtils;
    SvClusteringMethods mClusteringMethods;

    List<SvCluster> mClusters;
    private ChainFinder mChainFinder;
    private LinkFinder mLinkFinder;

    public static int SMALL_CLUSTER_SIZE = 3;

    private static final Logger LOGGER = LogManager.getLogger(ClusterAnalyser.class);

    public ClusterAnalyser(final SvClusteringConfig config, final SvUtilities utils, SvClusteringMethods clusteringMethods)
    {
        mConfig = config;
        mUtils = utils;
        mClusteringMethods = clusteringMethods;
        mClusters = null;
        mLinkFinder = new LinkFinder(mConfig, mUtils, mClusteringMethods);
        mChainFinder = new ChainFinder(mUtils);
        mChainFinder.setLogVerbose(mConfig.LogVerbose);
    }

    public void setClusters(List<SvCluster> clusters)
    {
        mClusters = clusters;
    }

    public void findSimpleCompleteChains(final String sampleId)
    {
        // for small clusters, try to find a full chain through all SVs
        for(SvCluster cluster : mClusters)
        {
            if(cluster.getCount() == 1 || cluster.getCount() > SMALL_CLUSTER_SIZE)
                continue;

            if(!cluster.isConsistent() || cluster.hasVariedCopyNumber())
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

    public void findLinksAndChains(final String sampleId)
    {
        for(SvCluster cluster : mClusters)
        {
            cluster.setUniqueBreakends();

            if(cluster.isSimpleSingleSV() || cluster.isFullyChained())
                continue;

            if(cluster.getUniqueSvCount() < 2)
            {
                setClusterArmBoundaries(cluster);
                continue;
            }

            // first establish links between SVs (eg TIs and DBs)
            mLinkFinder.findLinkedPairs(sampleId, cluster);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(sampleId, cluster);

            setClusterArmBoundaries(cluster);
            // setClusterResolvedState(cluster);

            cluster.logDetails();
        }

        // now look at merging unresolved & inconsistent clusters where they share the same chromosomal arms
        List<SvCluster> mergedClusters = mergeInconsistentClusters(mClusters);

        if(!mergedClusters.isEmpty())
        {
            for(SvCluster cluster : mergedClusters)
            {
                cluster.setDesc(cluster.getClusterTypesAsString());

                // need to be careful replicating already replicated SVs..
                // especially those already in linked chains
                // may be only replicate stand-alone SVs at this point which have now been clustered
                replicateMergedClusterSVs(sampleId, cluster);

                // repeat the search for inferred links now that additional SVs have been merged in but only on unlinked SVs
                List<SvLinkedPair> newLinkedPairs = mLinkFinder.createInferredLinkedPairs(cluster, cluster.getUnlinkedSVs(), true);

                cluster.getInferredLinkedPairs().addAll(newLinkedPairs);

                // createCopyNumberSegments(sampleId, cluster);

                findChains(sampleId, cluster);
            }

            // any clusters which were merged to resolve a collection of them, but
            // which did not lead to any longer chains, are now de-merged
            demergeClusters(mergedClusters);
        }

        for(SvCluster cluster : mClusters)
        {
            if(cluster.hasSubClusters()) // these haven't been logged
                cluster.logDetails();

            mLinkFinder.resolveTransitiveSVs(sampleId, cluster);
            cacheFinalLinkedPairs(sampleId, cluster);
        }
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

    private void replicateMergedClusterSVs(final String sampleId, SvCluster cluster)
    {
        // first remove any replication previously performed before clusters were merged
        if(!cluster.hasSubClusters())
            return;

        if(!cluster.hasVariedCopyNumber())
            return;

        int minCopyNumber = cluster.getMinCopyNumber();

        for(SvCluster subCluster : cluster.getSubClusters())
        {
            if(subCluster.hasReplicatedSVs())
                continue;

            // for now to difficult to consider the impact on these
            if(!subCluster.getChains().isEmpty() || !subCluster.getLinkedPairs().isEmpty())
                continue;

            int clusterCount = subCluster.getCount();

            for(int i = 0; i < clusterCount; ++i)
            {
                SvVarData var = subCluster.getSVs().get(i);
                int calcCopyNumber = var.impliedCopyNumber(true);

                if(calcCopyNumber <= minCopyNumber)
                    continue;

                int svMultiple = calcCopyNumber / minCopyNumber;

                LOGGER.debug("sample({}) replicating SV({}) {} times, copyNumChg({} vs min={})",
                        sampleId, var.posId(), svMultiple, calcCopyNumber, minCopyNumber);

                var.setReplicatedCount(svMultiple);

                // add to the parent cluster only for now
                for(int j = 1; j < svMultiple; ++j)
                {
                    SvVarData newVar = new SvVarData(var);
                    cluster.addVariant(newVar);
                }
            }

        }

        // cluster.removeReplicatedSvs();
    }

    private void setClusterResolvedState(SvCluster cluster)
    {
        if (cluster.isSimpleSVs())
        {
            cluster.setResolved(true, "Simple");
            return;
        }

        // next simple reciprocal inversions and translocations
        if (cluster.getCount() == 2 && cluster.isConsistent() && cluster.getChains().size() == 1)
        {
            cluster.setResolved(true, "Reciprocal");
            return;
        }

        // next clusters with which start and end on the same arm, have the same start and end orientation
        // and the same start and end copy number

    }

    private void setClusterArmBoundaries(SvCluster cluster)
    {
        // for each arm group within the cluster, find the bounding SV breakends
        // excluding any which are part of a continuous chain (ends excluded)
        final List<SvVarData> unlinkedSVs = cluster.getUnlinkedSVs();
        final List<SvChain> chains = cluster.getChains();

        for(final SvArmGroup armGroup : cluster.getArmGroups())
        {
            armGroup.setBreakend(null, true);
            armGroup.setBreakend(null, false);

            SvVarData startVar = null;
            SvVarData endVar = null;
            long startPosition = -1;
            long endPosition = 0;

            for(final SvVarData var : armGroup.getSVs())
            {
                if(var.isReplicatedSv())
                    continue;

                boolean checkStart = false;
                boolean checkEnd = false;

                if(unlinkedSVs.contains(var))
                {
                    if(var.type() == BND)
                    {
                        if(var.chromosome(true).equals(armGroup.chromosome()))
                            checkStart = true;
                        else if(var.chromosome(false).equals(armGroup.chromosome()))
                            checkEnd = true;
                    }
                    else
                    {
                        checkStart = true;
                        checkEnd = true;
                    }
                }
                else
                {
                    // check chain ends for a match with this SV
                    // translocations are skipped if their open end is on another arm
                    for(final SvChain chain : chains)
                    {
                        if(chain.getFirstSV().equals(var)
                        && (!var.isTranslocation() || var.chromosome(chain.firstLinkOpenOnStart()).equals(armGroup.chromosome())))
                        {
                            if(chain.firstLinkOpenOnStart())
                                checkStart = true;
                            else
                                checkEnd = true;
                        }

                        if(chain.getLastSV().equals(var)
                        && (!var.isTranslocation() || var.chromosome(chain.lastLinkOpenOnStart()).equals(armGroup.chromosome())))
                        {
                            if(chain.lastLinkOpenOnStart())
                                checkStart = true;
                            else
                                checkEnd = true;
                        }
                    }
                }

                if(checkStart && (startPosition < 0 || var.position(true) < startPosition))
                {
                    startVar = var;
                    startPosition = var.position(true);
                }

                if(checkEnd && !var.isNullBreakend() && var.position(false) > endPosition)
                {
                    endVar = var;
                    endPosition = var.position(false);
                }
            }

            if(startVar != null)
            {
                boolean useStart = startVar.type() == BND ? startVar.chromosome(true).equals(armGroup.chromosome()) : true;
                armGroup.setBreakend(new SvBreakend(startVar, useStart), true);
            }

            if(endVar != null)
            {
                boolean useStart = endVar.type() == BND ? endVar.chromosome(true).equals(armGroup.chromosome()) : false;
                armGroup.setBreakend(new SvBreakend(endVar, useStart), false);
            }

            if(cluster.getCount() > 1)
            {
                LOGGER.debug("cluster({}) arm({}) consistent({}) start({}) end({})",
                        cluster.getId(), armGroup.id(), armGroup.isConsistent(),
                        armGroup.getBreakend(true) != null ? armGroup.getBreakend(true).toString() : "null",
                        armGroup.getBreakend(false) != null ? armGroup.getBreakend(false).toString() : "null");
            }
        }
    }

    private boolean isConsistentCluster(final SvCluster cluster)
    {
        if(cluster.isSimpleSVs())
            return true;

        if(!cluster.isConsistent())
            return false;

        if(cluster.isFullyChained())
            return true;

        // other wise check whether all remaining SVs are either in consistent chains
        // or themselves consistent
        for(final SvChain chain : cluster.getChains())
        {
            if(!chain.isConsistent())
                return false;
        }

        List<SvVarData> varList = Lists.newArrayList();
        for(final SvVarData var : cluster.getUnlinkedSVs())
        {
            varList.clear();
            varList.add(var);
            if(calcConsistency(varList) != 0)
                return false;
        }

        return true;
    }

    private List<SvCluster> mergeInconsistentClusters(List<SvCluster> clusters)
    {
        // it's possible that to resolve arms and more complex arrangements, clusters not merged
        // by proximity of overlaps must be put together to solve inconsistencies (ie loose ends)

        List<SvCluster> mergedClusters = Lists.newArrayList();

        // merge any cluster which is itself not consistent and has breakends on the same arm as another
        int index1 = 0;
        while(index1 < clusters.size())
        {
            SvCluster cluster1 = clusters.get(index1);

            if(isConsistentCluster(cluster1))
            // if(cluster1.isSimpleSVs() || (cluster1.isConsistent() && cluster1.isFullyChained()))
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

                if(isConsistentCluster(cluster2))
                // if(cluster2.isSimpleSVs() || (cluster2.isConsistent() && cluster2.isFullyChained()))
                {
                    ++index2;
                    continue;
                }

                boolean foundConnection = canMergeClustersOnOverlaps(cluster1, cluster2);

                if(!foundConnection)
                {
                    ++index2;
                    continue;
                }

                boolean cluster2Merged = false;

                if(cluster1.hasSubClusters())
                {
                    LOGGER.debug("cluster({} svs={}) merges in cluster({} svs={})",
                            cluster1.getId(), cluster1.getCount(), cluster2.getId(), cluster2.getCount());

                    cluster2Merged = true;
                    cluster1.addSubCluster(cluster2);
                    setClusterArmBoundaries(cluster1);
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
                    mergedClusters.add(newCluster);
                    setClusterArmBoundaries(newCluster);
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

        return mergedClusters;
    }

    private boolean canMergeClustersOnOverlaps(SvCluster cluster1, SvCluster cluster2)
    {
        // checks for overlapping breakends in inconsistent matching arms
        final List<SvArmGroup> armGroups1 = cluster1.getArmGroups();
        final List<SvArmGroup> armGroups2 = cluster2.getArmGroups();

        for (SvArmGroup armGroup1 : armGroups1)
        {
            if(armGroup1.isConsistent())
                continue;

            for (SvArmGroup armGroup2 : armGroups2)
            {
                if(armGroup2.isConsistent())
                    continue;

                if(!armGroup1.matches(armGroup2))
                    continue;

                // for now merge any inconsistent arm
                LOGGER.debug("inconsistent cluster({}) and cluster({}) linked on chrArm({})",
                        cluster1.getId(), cluster2.getId(), armGroup1.id());

                return true;

                /*
                for(int i = 0; i < 2; ++i)
                {
                    final SvArmGroup group1 = (i==0) ? armGroup1 : armGroup2;
                    final SvArmGroup group2 = (i==0) ? armGroup2 : armGroup1;

                    if(group2.getBreakend(true) == null || group2.getBreakend(false) == null)
                        continue;

                    for(int be = SVI_START; be <= SVI_END; ++be)
                    {
                        boolean useStart = isStart(be);

                        if(group1.getBreakend(useStart) == null)
                            continue;

                        if(group1.getBreakend(useStart).position() >= group2.getBreakend(true).position()
                        && group1.getBreakend(useStart).position() <= group2.getBreakend(false).position())
                        {
                            LOGGER.debug("inconsistent cluster({}) and cluster({}) linked on chrArm({})",
                                    cluster1.getId(), cluster2.getId(), armGroup1.id());

                            return true;
                        }
                    }
                }
                    */
            }
        }

        return false;
    }

    private boolean canMergeClustersOnSameArm(SvCluster cluster1, SvCluster cluster2)
    {
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

                    return true;
                }
            }
        }

        return false;
    }

    private void demergeClusters(List<SvCluster> clusters)
    {
        // de-merge any clusters which didn't form longer chains
        int clusterCount = clusters.size();
        for(int i = 0; i < clusterCount;)
        {
            SvCluster cluster = clusters.get(i);

            if (!cluster.hasSubClusters() || cluster.isFullyChained())
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

    public void markFoldbacks()
    {
        for(final Map.Entry<String, List<SvBreakend>> entry : mClusteringMethods.getChrBreakendMap().entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            for(int i = 1; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvBreakend prevBreakend = breakendList.get(i - 1);

                checkFoldbackBreakends(breakend, prevBreakend);
            }
        }
    }

    private void checkFoldbackBreakends(SvBreakend be1, SvBreakend be2)
    {
        // consecutive breakends, same orientation, same var or part of a chain
        if(be1.orientation() != be2.orientation())
            return;

        final SvVarData var1 = be1.getSV();
        final SvVarData var2 = be2.getSV();

        if(var1.type() == INS || var2.type() == INS)
            return;

        // skip unclustered DELs & DUPs, reciprocal INV or reciprocal BNDs
        final SvCluster cluster1 = findCluster(var1, mClusters);

        if(cluster1.isSimpleSVs() || (cluster1.getCount() == 2 && cluster1.isConsistent()))
            return;

        final SvCluster cluster2 = findCluster(var2, mClusters);

        if(cluster2.isSimpleSVs() || (cluster2.getCount() == 2 && cluster2.isConsistent()))
            return;

        if(!var1.equals(var2))
        {
            // must be same cluster and part of the same chain
            if(cluster1 != cluster2)
                return;

            if(var1.getReplicatedCount() != var2.getReplicatedCount())
                return;

            final SvChain chain1 = cluster1.findChain(var1);
            final SvChain chain2 = cluster2.findChain(var1);

            if(chain1 == null || chain2 == null || chain1 != chain2)
                return;
        }

        boolean v1Start = be1.usesStart();
        boolean v2Start = be2.usesStart();

        // check copy numbers match
        double cn1 = 0;
        double cn2 = 0;
        if((be1.orientation() == 1 && be1.position() < be2.position()) || (be1.orientation() == -1 && be1.position() > be2.position()))
        {
            // be1 is facing away from be2, so need to take its copy number on the break side
            cn1 = var1.copyNumber(v1Start) - var1.copyNumberChange(v1Start);
            cn2 = var2.copyNumber(v2Start);
        }
        else
        {
            cn2 = var2.copyNumber(v2Start) - var2.copyNumberChange(v2Start);
            cn1 = var1.copyNumber(v1Start);
        }

        if(!copyNumbersEqual(cn1, cn2))
            return;

        int length = (int)abs(be1.position() - be2.position());

        // if either variant already has foldback info set, favour
        // a) simple inversions then
        // b) shortest length

        boolean skipFoldback = false;
        if(!var1.getFoldbackLink(v1Start).isEmpty() && !var1.equals(var2) && var1.getFoldbackLen(v1Start) < length)
        {
            skipFoldback = true;
        }
        else if(!var2.getFoldbackLink(v2Start).isEmpty() && !var1.equals(var2) && var2.getFoldbackLen(v2Start) < length)
        {
            skipFoldback = true;
        }

        if(skipFoldback)
            return;

        if(!var1.getFoldbackLink(v1Start).isEmpty())
            clearFoldbackInfo(var1.getFoldbackLink(v1Start), var1.id(), cluster1, v1Start);

        if(!var2.getFoldbackLink(v2Start).isEmpty())
            clearFoldbackInfo(var2.getFoldbackLink(v2Start), var2.id(), cluster2, v2Start);

        var1.setFoldbackLink(v1Start, var2.id(), length);
        var2.setFoldbackLink(v2Start, var1.id(), length);

        if(var1.equals(var2))
        {
            LOGGER.debug(String.format("cluster(%s) foldback inversion SV(%s) length(%d) copyNumber(%.3f)",
                    cluster1.getId(), var1.posId(), length, cn1));
        }
        else
        {
            LOGGER.debug(String.format("cluster(%s) foldback be1(%s) be2(%s) length(%d) copyNumber(%.3f)",
                    cluster1.getId(), be1.toString(), be2.toString(), length, cn1));
        }
    }

    private void clearFoldbackInfo(final String varId, final String matchVarId, SvCluster cluster, boolean useStart)
    {
        SvVarData var = findVariantById(varId, cluster.getSVs());

        if(var == null)
            return;

        if(var.getFoldbackLink(true).equals(matchVarId))
            var.setFoldbackLink(true, "", -1);
        else
            var.setFoldbackLink(false, "", -1);
    }

    private void createCopyNumberSegments(final String sampleId, final SvCluster cluster)
    {
        int cnId = 0;
        final List<SvVarData> clusterSVs = cluster.getSVs();

        Map<String, List<SvCNData>> chrCNDataMap = new HashMap();

        for(Map.Entry<String, List<SvBreakend>> entry : mClusteringMethods.getChrBreakendMap().entrySet())
        {
            final String chromosome = entry.getKey();

            List<SvBreakend> breakendList = entry.getValue();

            List<SvCNData> copyNumberData = Lists.newArrayList();

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvVarData var = breakend.getSV();

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
            LOGGER.debug("sample({}) cluster({}: {} count={}) has {} linked pairs",
                    sampleId, cluster.getId(), cluster.getDesc(), cluster.getUniqueSvCount(), linkedPairs.size());

            /*
            if(LOGGER.isDebugEnabled())
            {
                for (final SvLinkedPair pair : linkedPairs)
                {
                    LOGGER.debug("linked {} {} pair length({}) variants({})",
                            pair.isInferred() ? "inferred" : "assembly", pair.linkType(), pair.length(), pair.toString());
                }
            }
            */
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
