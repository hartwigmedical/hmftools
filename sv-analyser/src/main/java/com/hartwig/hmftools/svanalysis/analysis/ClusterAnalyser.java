package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_COMMON_ARMS;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_FOLDBACKS;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_SOLO_SINGLE;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.addClusterReason;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.checkClusterDuplicates;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.isConsistentCluster;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_COMPLEX_CHAIN;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SGL_PAIR_DEL;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SGL_PAIR_DUP;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SIMPLE_CHAIN;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SGL_PAIR_INS;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SIMPLE_SV;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_SGL;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.findVariantById;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_INFER_ONLY;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.haveSameChrArms;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ClusterAnalyser {

    final SvClusteringConfig mConfig;
    final SvUtilities mUtils;
    SvClusteringMethods mClusteringMethods;

    String mSampleId;
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
        mSampleId = "";
        mLinkFinder = new LinkFinder(mConfig, mUtils, mClusteringMethods);
        mChainFinder = new ChainFinder(mUtils);
        mChainFinder.setLogVerbose(mConfig.LogVerbose);
    }

    public void setClusterData(final String sampleId, List<SvCluster> clusters)
    {
        mSampleId = sampleId;
        mClusters = clusters;
    }

    public void findSimpleCompleteChains()
    {
        // for small clusters, try to find a full chain through all SVs
        for(SvCluster cluster : mClusters)
        {
            if(cluster.getCount() == 1 && cluster.isSimpleSVs())
            {
                setClusterResolvedState(cluster);
                continue;
            }

            // skip more complicated clusters for now
            if(cluster.getCount() > SMALL_CLUSTER_SIZE || !cluster.isConsistent() || cluster.hasVariedCopyNumber())
                continue;

            // first establish links between SVs (eg TIs and DBs)
            mLinkFinder.findLinkedPairs(mSampleId, cluster);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster);

            cacheFinalLinkedPairs(cluster);

            setClusterResolvedState(cluster);

            if(isConsistentCluster(cluster))
            {
                LOGGER.debug("sample({}) cluster({}) simple and consistent with {} SVs", mSampleId, cluster.getId(), cluster.getCount());
            }
        }
    }

    public void findLinksAndChains()
    {
        for(SvCluster cluster : mClusters)
        {
            if(cluster.isSimpleSingleSV())
                continue;

            if(cluster.isFullyChained())
                continue;

            if(cluster.getUniqueSvCount() < 2)
            {
                setClusterArmBoundaries(cluster);
                continue;
            }

            // first establish links between SVs (eg TIs and DBs)
            mLinkFinder.findLinkedPairs(mSampleId, cluster);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster);

            setClusterArmBoundaries(cluster);

            cacheFinalLinkedPairs(cluster);

            setClusterResolvedState(cluster);

            cluster.logDetails();
        }

        // INVs and other SV-pairs which make foldbacks are now used in the inconsistent clustering logic
        markFoldbacks();

        // now look at merging unresolved & inconsistent clusters where they share the same chromosomal arms
        List<SvCluster> mergedClusters = Lists.newArrayList();

        mergedClusters = mergeInconsistentClusters(mergedClusters);
        boolean foundMerges = !mergedClusters.isEmpty();

        while(foundMerges)
        {
            List<SvCluster> newMergedClusters = mergeInconsistentClusters(mergedClusters);
            foundMerges = !newMergedClusters.isEmpty();

            for(SvCluster cluster : newMergedClusters)
            {
                if(!mergedClusters.contains(cluster))
                    mergedClusters.add(cluster);
            }
        }

        if(!mergedClusters.isEmpty())
        {
            checkClusterDuplicates(mergedClusters);

            for(SvCluster cluster : mergedClusters)
            {
                cluster.setDesc(cluster.getClusterTypesAsString());

                // need to be careful replicating already replicated SVs..
                // especially those already in linked chains
                // may be only replicate stand-alone SVs at this point which have now been clustered
                replicateMergedClusterSVs(cluster);

                // repeat the search for inferred links now that additional SVs have been merged in but only on unlinked SVs
                List<SvLinkedPair> newLinkedPairs = mLinkFinder.createInferredLinkedPairs(cluster, cluster.getUnlinkedSVs(), true);

                cluster.getInferredLinkedPairs().addAll(newLinkedPairs);

                // createCopyNumberSegments(sampleId, cluster);

                findChains(cluster);
            }

            /*
            // any clusters which were merged to resolve a collection of them, but
            // which did not lead to any longer chains, are now de-merged
            demergeClusters(mergedClusters);
            */

            for(SvCluster cluster : mergedClusters)
            {
                cacheFinalLinkedPairs(cluster);

                setClusterResolvedState(cluster);

                cluster.logDetails();
            }
        }

        // checkClusterDuplicates(mClusters);

        for(SvCluster cluster : mClusters)
        {
            mLinkFinder.resolveTransitiveSVs(mSampleId, cluster);
            // cacheFinalLinkedPairs(cluster);

            // reportChainFeatures(cluster);
        }
    }

    private void findChains(SvCluster cluster)
    {
        mChainFinder.initialise(mSampleId, cluster);
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

    private void replicateMergedClusterSVs(SvCluster cluster)
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

                LOGGER.debug("cluster({}) replicating SV({}) {} times, copyNumChg({} vs min={})",
                        cluster.getId(), var.posId(), svMultiple, calcCopyNumber, minCopyNumber);

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
        if(cluster.isResolved())
            return;

        if (cluster.isSimpleSVs())
        {
            if(cluster.getCount() == 1)
            {
                cluster.setResolved(true, RESOLVED_TYPE_SIMPLE_SV);
                return;
            }

            if(cluster.getTypeCount(DEL) + cluster.getTypeCount(DUP) == 2)
            {
                if(mClusteringMethods.markDelDupPairTypes(cluster, mSampleId))
                    return;
            }

            cluster.setResolved(true, RESOLVED_TYPE_SIMPLE_SV);
            return;
        }

        // next simple reciprocal inversions and translocations
        if (cluster.getCount() == 2 && cluster.isConsistent())
        {
            if(cluster.getTypeCount(BND) == 2)
            {
                mClusteringMethods.markBndPairTypes(cluster, mSampleId);
            }
            else if(cluster.getTypeCount(INV) == 2)
            {
                mClusteringMethods.markInversionPairTypes(cluster, mSampleId);
            }
            else if(cluster.getTypeCount(SGL) == 2)
            {
                if(cluster.getLinkedPairs().size() == 1 && cluster.getLinkedPairs().get(0).linkType() == LINK_TYPE_SGL)
                {
                    cluster.setResolved(true, RESOLVED_TYPE_SGL_PAIR_INS);
                }
                else
                {
                    final SvVarData var1 = cluster.getSVs().get(0);
                    final SvVarData var2 = cluster.getSVs().get(1);

                    if(var1.orientation(true) != var2.orientation(true))
                    {
                        boolean ssFirst = var1.position(true) < var2.position(true);
                        boolean ssPosOrientation = var1.orientation(true) == 1;

                        if(ssFirst == ssPosOrientation)
                            cluster.setResolved(true, RESOLVED_TYPE_SGL_PAIR_DEL);
                        else
                            cluster.setResolved(true, RESOLVED_TYPE_SGL_PAIR_DUP);
                    }
                }
            }

            return;
        }

        // next clusters with which start and end on the same arm, have the same start and end orientation
        // and the same start and end copy number
        if(isConsistentCluster(cluster) || cluster.isFullyChained())
        {
            boolean isResolved = cluster.isConsistent();

            if(!cluster.hasReplicatedSVs())
                cluster.setResolved(isResolved, RESOLVED_TYPE_SIMPLE_CHAIN);
            else
                cluster.setResolved(isResolved, RESOLVED_TYPE_COMPLEX_CHAIN);

            return;
        }

        if(cluster.hasLinkingLineElements())
        {
            // skip further classification for now
            cluster.setResolved(false, "Line");
            return;
        }

        // cluster remains largely unresolved..
        if(!cluster.getChains().isEmpty())
        {
            if(!cluster.hasReplicatedSVs())
                cluster.setResolved(false, "SimplePartialChain");
            else
                cluster.setResolved(false, "ComplexPartialChain");
        }
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

                for(int be = SVI_START; be <= SVI_END; ++be)
                {
                    if((!checkStart && be == SVI_START) || (!checkEnd && be == SVI_END))
                        continue;

                    boolean useStart = isStart(be);

                    long position = var.position(useStart);

                    if(startPosition < 0 || position < startPosition)
                    {
                        startVar = var;
                        startPosition = position;
                    }

                    if(position > endPosition)
                    {
                        endVar = var;
                        endPosition = position;
                    }
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

            /*
            if(cluster.getCount() > 1)
            {
                LOGGER.debug("cluster({}) arm({}) consistent({}) start({}) end({}) posBoundaries({} -> {})",
                        cluster.getId(), armGroup.id(), armGroup.isConsistent(),
                        armGroup.getBreakend(true) != null ? armGroup.getBreakend(true).toString() : "null",
                        armGroup.getBreakend(false) != null ? armGroup.getBreakend(false).toString() : "null",
                        armGroup.posStart(), armGroup.posEnd());
            }
            */
        }
    }

    private List<SvCluster> mergeInconsistentClusters(List<SvCluster> existingMergedClusters)
    {
        // it's possible that to resolve arms and more complex arrangements, clusters not merged
        // by proximity of overlaps must be put together to solve inconsistencies (ie loose ends)
        List<SvCluster> mergedClusters = Lists.newArrayList();

        // merge any cluster which is itself not consistent and has breakends on the same arm as another
        int index1 = 0;
        while(index1 < mClusters.size())
        {
            SvCluster cluster1 = mClusters.get(index1);

            if(cluster1.isResolved())
            {
                ++index1;
                continue;
            }

            boolean isConsistent1 = isConsistentCluster(cluster1);
            // boolean hasLongDelDup1 = mClusteringMethods.clusterHasLongDelDup(cluster1);
            boolean hasOpenSingle1 = hasOneOpenSingleVariant(cluster1);
            boolean isSoloSingle1 = cluster1.getCount() == 1 && cluster1.getSVs().get(0).type() == SGL;

            if(isConsistent1  && !hasOpenSingle1 && !isSoloSingle1) // && !hasLongDelDup1
            {
                ++index1;
                continue;
            }

            boolean cluster1Merged = false;
            SvCluster newCluster = null;

            int index2 = index1 + 1;
            while(index2 < mClusters.size())
            {
                SvCluster cluster2 = mClusters.get(index2);

                if(cluster2.isResolved())
                {
                    ++index2;
                    continue;
                }

                boolean isConsistent2 = isConsistentCluster(cluster2);

                boolean foundConnection = false;

                if(!isConsistent1 && !isConsistent2)
                {
                    foundConnection = canMergeClustersOnFoldbacks(cluster1, cluster2);

                    if(!foundConnection)
                        foundConnection = canMergeClustersOnCommonArms(cluster1, cluster2);
                }

                /*
                if(!foundConnection && hasLongDelDup1)
                {
                    foundConnection = mClusteringMethods.canMergeClustersOnLongDelDups(cluster1, cluster2);
                }

                if(!foundConnection && mClusteringMethods.clusterHasLongDelDup(cluster2))
                {
                    foundConnection = mClusteringMethods.canMergeClustersOnLongDelDups(cluster2, cluster1);
                }
                */

                if(!foundConnection)
                {
                    boolean hasOpenSingle2 = hasOneOpenSingleVariant(cluster2);
                    boolean isSoloSingle2 = cluster2.getCount() == 1 && cluster2.getSVs().get(0).type() == SGL;

                    if(hasOpenSingle1 && isSoloSingle2 && canMergeOpenSingles(cluster1, cluster2))
                    {
                        foundConnection = true;
                        addClusterReason(cluster2, CLUSTER_REASON_SOLO_SINGLE, "");
                    }

                    if(!foundConnection && hasOpenSingle2 && isSoloSingle1 && canMergeOpenSingles(cluster2, cluster1))
                    {
                        foundConnection = true;
                        addClusterReason(cluster1, CLUSTER_REASON_SOLO_SINGLE, "");
                    }

                    if(!foundConnection && isSoloSingle1 && isSingleClosestTI(cluster2, cluster2.getSVs().get(0)))
                    {
                        foundConnection = true;
                        addClusterReason(cluster1, CLUSTER_REASON_SOLO_SINGLE, "");
                    }

                    if(!foundConnection && isSoloSingle2 && isSingleClosestTI(cluster1, cluster2.getSVs().get(0)))
                    {
                        foundConnection = true;
                        addClusterReason(cluster2, CLUSTER_REASON_SOLO_SINGLE, "");
                    }
                }

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

                    if(cluster2.hasSubClusters())
                    {
                        for(SvCluster subCluster : cluster2.getSubClusters())
                        {
                            cluster1.addSubCluster(subCluster);
                        }

                        existingMergedClusters.remove(cluster2);
                    }
                    else
                    {
                        cluster1.addSubCluster(cluster2);
                    }

                    setClusterArmBoundaries(cluster1);
                }
                else
                {
                    cluster1Merged = true;
                    cluster2Merged = true;

                    newCluster = new SvCluster(mClusteringMethods.getNextClusterId());

                    LOGGER.debug("new cluster({}) from merge of cluster({} svs={}) and cluster({} svs={})",
                            newCluster.getId(), cluster1.getId(), cluster1.getCount(), cluster2.getId(), cluster2.getCount());

                    newCluster.addSubCluster(cluster1);

                    if(cluster2.hasSubClusters())
                    {
                        for(SvCluster subCluster : cluster2.getSubClusters())
                        {
                            newCluster.addSubCluster(subCluster);
                        }

                        existingMergedClusters.remove(cluster2);
                    }
                    else
                    {
                        newCluster.addSubCluster(cluster2);
                    }

                    mergedClusters.add(newCluster);
                    setClusterArmBoundaries(newCluster);
                }

                if(cluster2Merged)
                    mClusters.remove(index2);
                else
                    ++index2;

                if(cluster1Merged)
                    break;
            }

            if(cluster1Merged && newCluster != null)
            {
                // cluster has been replaced with a combined one
                mClusters.remove(index1);
                mClusters.add(index1, newCluster);
            }
            else
            {
                ++index1;
            }
        }

        return mergedClusters;
    }

    private boolean hasOneOpenSingleVariant(final SvCluster cluster)
    {
        if(cluster.getTypeCount(SGL) != 1)
            return false;

        if(!cluster.getUnlinkedSVs().isEmpty())
            return false;

        for(final SvChain chain : cluster.getChains())
        {
            if(chain.getFirstSV().type() == SGL || chain.getLastSV().type() == SGL)
                return true;
        }

        return false;
    }

    private boolean canMergeOpenSingles(final SvCluster otherCluster, final SvCluster soloSingleCluster)
    {
        // first cluster has the open single in a chain, the second is a solo single
        // can be merged if the open single is the closest cluster to this cluster
        final SvArmGroup singleArmGroup = soloSingleCluster.getArmGroups().get(0);
        final SvVarData soloVar = soloSingleCluster.getSVs().get(0);

        final List<SvBreakend> breakendList = mClusteringMethods.getChrBreakendMap().get(singleArmGroup.chromosome());

        for (int i = 1; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);

            if(breakend.getSV().equals(soloVar))
            {
                // check if the preceding or next breakend is in the cluster with the open single chain
                if((i > 0 && otherCluster.getSVs().contains(breakendList.get(i - 1).getSV()))
                || (i < breakendList.size() - 1 && otherCluster.getSVs().contains(breakendList.get(i + 1).getSV())))
                {
                    LOGGER.debug("inconsistent cluster({}) and cluster({}) have adjacent SVs with solo-single",
                            otherCluster.getId(), soloSingleCluster.getId());

                    return true;
                }

                break;
            }
        }

        return false;
    }

    private boolean isSingleClosestTI(final SvCluster cluster, final SvVarData soloVar)
    {
        // first cluster has the open single in a chain, the second is a solo single
        // can be merged if the open single is the closest cluster to this cluster
        final SvArmGroup singleArmGroup = cluster.getArmGroups().get(0);

        final List<SvBreakend> breakendList = mClusteringMethods.getChrBreakendMap().get(singleArmGroup.chromosome());

        for (int i = 1; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);

            if(breakend.getSV().equals(soloVar))
            {
                // check if the preceding or next breakend is in the cluster with the open single chain
                final SvBreakend prevBreakend = (i > 0) ? breakendList.get(i - 1) : null;
                final SvBreakend nextBreakend = (i < breakendList.size() - 1) ? breakendList.get(i + 1) : null;
                boolean checkPrev = prevBreakend != null && cluster.getSVs().contains(prevBreakend.getSV());
                boolean checkNext = nextBreakend != null && cluster.getSVs().contains(nextBreakend.getSV());

                if((checkPrev && canFormTIWithChainEnd(cluster.getChains(), prevBreakend.getSV(), soloVar))
                || (checkNext && canFormTIWithChainEnd(cluster.getChains(), nextBreakend.getSV(), soloVar)))
                {
                    LOGGER.debug("inconsistent cluster({}) can form TI with soloSingle({})",
                            cluster.getId(), soloVar.id());
                    return true;
                }

                return false;
            }
        }

        return false;
    }

    private boolean canFormTIWithChainEnd(final List<SvChain> chains, final SvVarData chainVar, final SvVarData soloVar)
    {
        for(final SvChain chain : chains)
        {
            if (chain.getFirstSV().equals(chainVar) && areLinkedSection(chainVar, soloVar, chain.firstLinkOpenOnStart(), true))
            {
                return true;
            }
            else if (chain.getLastSV().equals(chainVar) && areLinkedSection(chain.getLastSV(), soloVar, chain.lastLinkOpenOnStart(), true))
            {
                return true;
            }
        }

        return false;
    }

    private boolean canMergeClustersOnFoldbacks(final SvCluster cluster1, final SvCluster cluster2)
    {
        final List<SvVarData> cluster1Foldbacks = cluster1.getFoldbacks();
        final List<SvVarData> cluster2Foldbacks = cluster2.getFoldbacks();

        for (final SvVarData var1 : cluster1Foldbacks)
        {
            if (var1.inLineElement())
                continue;

            for (int be1 = SVI_START; be1 <= SVI_END; ++be1)
            {
                boolean v1Start = isStart(be1);

                if(be1 == SVI_END && var1.type() != BND)
                    continue;

                for (final SvVarData var2 : cluster2Foldbacks)
                {
                    for (int be2 = SVI_START; be2 <= SVI_END; ++be2)
                    {
                        boolean v2Start = isStart(be2);

                        if (be2 == SVI_END && var2.type() != BND)
                            continue;

                        if (var2.inLineElement())
                            continue;

                        if (!var1.chromosome(v1Start).equals(var2.chromosome(v2Start)) || !var1.arm(v1Start).equals(var2.arm(v2Start)))
                            continue;

                        LOGGER.debug("cluster({}) SV({}) and cluster({}) SV({}) have foldbacks on same arm",
                                cluster1.getId(), var1.posId(), cluster2.getId(), var2.posId());

                        addClusterReason(cluster1, CLUSTER_REASON_FOLDBACKS, var2.id());
                        addClusterReason(cluster2, CLUSTER_REASON_FOLDBACKS, var1.id());
                        return true;
                    }
                }
            }
        }

        return false;
    }

    private boolean canMergeClustersOnCommonArms(final SvCluster cluster1, final SvCluster cluster2)
    {
        // merge if the 2 clusters have BNDs linking the same 2 arms, but not included any BNDs in the middle
        // of chains or those in line elements
        if(cluster1.getTypeCount(BND) == 0 || cluster2.getTypeCount(BND) == 0)
            return false;

        // int matchedArms = 0;
        // String linkedSvStr = "";

        // check that the BNDs in these 2 common arm links are either unlinked or the ends of chains
        final List<SvVarData> bndList1 = getUnlinkedTranslocation(cluster1);
        final List<SvVarData> bndList2 = getUnlinkedTranslocation(cluster2);

        for (final SvVarData var1 : bndList1)
        {
            for (final SvVarData var2 : bndList2)
            {
                if(haveSameChrArms(var1, var2))
                {
                    LOGGER.debug("cluster({}) and cluster({}) have common links with SV({}) and SV({})",
                            cluster1.getId(), cluster2.getId(), var1.posId(), var2.posId());

                    final String commonArms = var1.id() + "_" + var2.id();

                    addClusterReason(cluster1, CLUSTER_REASON_COMMON_ARMS, commonArms);
                    addClusterReason(cluster2, CLUSTER_REASON_COMMON_ARMS, commonArms);
                    return true;
                }
            }
        }

        return false;
    }

    private final List<SvVarData> getUnlinkedTranslocation(final SvCluster cluster)
    {
        List<SvVarData> bndList = Lists.newArrayList();

        for(final SvVarData var : cluster.getUnlinkedSVs())
        {
            if(var.type() == BND && !var.inLineElement() && !var.isReplicatedSv())
            {
                bndList.add(var);
            }
        }

        for(final SvChain chain : cluster.getChains())
        {
            if(chain.getFirstSV().type() == BND && !chain.getFirstSV().inLineElement())
            {
                bndList.add(chain.getFirstSV());
            }

            if(chain.getLastSV().type() == BND && !chain.getLastSV().inLineElement())
            {
                bndList.add(chain.getLastSV());
            }
        }

        return bndList;
    }

    private void demergeClusters(List<SvCluster> mergedClusters)
    {
        // de-merge any clusters which didn't form longer chains
        int clusterCount = mergedClusters.size();

        int index = 0;
        while(index < mergedClusters.size())
        {
            SvCluster cluster = mergedClusters.get(index);

            if (cluster.isFullyChained())
            {
                ++index;
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
                ++index;
                continue;
            }

            // add the original clusters back in
            for(final SvCluster subCluster : cluster.getSubClusters())
            {
                if(mClusters.contains(subCluster))
                {
                    continue;
                }

                mClusters.add(subCluster);
            }

            if(mainChainCount > 0)
            {
                LOGGER.debug("removed cluster({}) since maxChainCount({}) less than subclusters({}) maxSubClusterChainCount({})",
                        cluster.getId(), mainChainCount, cluster.getSubClusters().size(), maxSubClusterChainCount);
            }

            mClusters.remove(cluster);
            mergedClusters.remove(index);
        }
    }

    public void reportChainFeatures(final SvCluster cluster)
    {
        for (final SvChain chain : cluster.getChains())
        {
            if(chain.getLinkedPairs().size() <= 1)
                continue;

            findChainRepeatedSegments(cluster, chain);

            findChainTranslocationTIs(cluster, chain);
        }
    }

    private void findChainRepeatedSegments(final SvCluster cluster, final SvChain chain)
    {
        if(!chain.hasReplicatedSVs())
            return;

        List<SvVarData> replicatedSVs = Lists.newArrayList();

        final List<SvVarData> svList = chain.getSvList();

        for (int i = 0; i < svList.size(); ++i)
        {
            final SvVarData var1 = svList.get(i);

            if(replicatedSVs.contains(var1))
                continue;

            for (int j = i + 1; j < svList.size(); ++j)
            {
                final SvVarData var2 = svList.get(j);

                if (!var1.equals(var2, true))
                    continue;

                replicatedSVs.add(var1);

                // look for repeated sections forwards or backwards from this point
                List<SvVarData> forwardRepeats = getRepeatedSvSequence(svList, i, j, true);

                boolean forwardSequence = false;

                if(!forwardRepeats.isEmpty())
                {
                    forwardSequence = true;
                    replicatedSVs.addAll(forwardRepeats);
                }
                else
                {
                    forwardSequence = false;
                    forwardRepeats = getRepeatedSvSequence(svList, i, j, false);
                }

                if(!forwardRepeats.isEmpty())
                {
                    replicatedSVs.addAll(forwardRepeats);

                    forwardRepeats.set(0, var1);

                    String svIds = var1.id();
                    for(int k = 1; k < forwardRepeats.size(); ++k)
                        svIds += ";" + forwardRepeats.get(k).id();

                    LOGGER.info("sample({}) cluster({}) chain({}) {} sequence of {} SVs starting at index({}:{}) SV({})",
                            mSampleId, cluster.getId(), chain.getId(), forwardSequence ? "forward" : "reverse",
                            forwardRepeats.size(), i, j, var1.id());

                    // ClusterId,ChainId,SequenceCount,VarIds,MatchDirection
                    LOGGER.info("CF_REPEAT_SEQ: {},{},{},{},{},{}",
                            mSampleId, cluster.getId(), chain.getId(), forwardRepeats.size(), svIds, forwardSequence);

                    break;
                }

                // no sequence found
            }
        }
    }

    private List<SvVarData> getRepeatedSvSequence(final List<SvVarData> svList, int firstIndex, int secondIndex, boolean walkForwards)
    {
        // walk forward from these 2 start points comparing SVs
        List<SvVarData> sequence = Lists.newArrayList();

        int i = firstIndex;
        int j = secondIndex + 1;

        if(walkForwards)
            ++i;
        else
            --i;

        while(i < secondIndex && i >= 0 && j < svList.size())
        {
            final SvVarData var1 = svList.get(i);
            final SvVarData var2 = svList.get(j);

            if(!var1.equals(var2, true))
                break;

            sequence.add(var1);

            ++j;

            if(walkForwards)
                ++i;
            else
                --i;
        }

        return sequence;
    }

    private void findChainTranslocationTIs(final SvCluster cluster, final SvChain chain)
    {
        final List<SvLinkedPair> linkedPairs = chain.getLinkedPairs();

        int svCount = chain.getSvCount();
        int uniqueSvCount = chain.getUniqueSvCount();

        boolean hasReplication = uniqueSvCount < svCount;

        final List<SvVarData> svList = chain.getSvList();

        for(int i = 0; i < linkedPairs.size(); ++i)
        {
            final SvLinkedPair pair = linkedPairs.get(i);
            final SvVarData svBack = svList.get(i);

            if(svList.size() == i+1)
                break; // chain loops back to same SV (whether closed or not

            final SvVarData svForward = svList.get(i+1);

            // find out-and-back or out-and-on translocations showing evidence of foreign fragment insertions
            if(svForward.type() == BND && i < linkedPairs.size() - 1 && i < svList.size() - 2)
            {
                final SvLinkedPair nextPair = linkedPairs.get(i+1);
                final SvVarData nextSvForward = svList.get(i+2);

                if(nextSvForward.type() == BND)
                {
                    boolean svForwardLinkedOnStart = pair.getLinkedOnStart(svForward);
                    final String startChromosome = svForward.chromosome(svForwardLinkedOnStart);
                    // final String linkChromosome = svForward.chromosome(!svForwardLinkedOnStart);

                    boolean nextSvForwardLinkedOnStart = nextPair.getLinkedOnStart(nextSvForward);
                    final String endChromosome = nextSvForward.chromosome(!nextSvForwardLinkedOnStart);

                    boolean outAndBack = startChromosome.equals(endChromosome);

                    nextPair.setInfo("TransTI");

                    // SampleId, ClusterId,ChainId,SvId1,SvId2,IsOutAndBack,TILength,IsAssembly,ChainLinks,HasReplication
                    LOGGER.info("CF_TRANS_TI: {},{},{},{},{},{},{},{},{},{}",
                            mSampleId, cluster.getId(), chain.getId(), svForward.id(), nextSvForward.id(),
                            outAndBack, nextPair.length(), !nextPair.isInferred(), chain.getLinkCount(), hasReplication);
                }
            }
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

                checkReplicatedBreakendFoldback(breakend);
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
        final SvCluster cluster1 = var1.getCluster();

        if(cluster1.isSimpleSVs() || (cluster1.getCount() == 2 && cluster1.isConsistent()))
            return;

        SvCluster cluster2 = null;

        if(var1.equals(var2))
        {
            cluster2 = cluster1;
        }
        else
        {
            cluster2 = var2.getCluster();

            if (cluster2.isSimpleSVs() || (cluster2.getCount() == 2 && cluster2.isConsistent()))
                return;
        }

        boolean v1Start = be1.usesStart();
        boolean v2Start = be2.usesStart();

        if(!var1.equals(var2))
        {
            // must be same cluster and part of the same chain
            if(cluster1 != cluster2)
                return;

            if(var1.getReplicatedCount() != var2.getReplicatedCount())
                return;

            final SvChain chain1 = cluster1.findChain(var1);
            final SvChain chain2 = cluster2.findChain(var2);

            if(chain1 == null || chain2 == null || chain1 != chain2)
                return;

            // check if a path can be walked between these 2 breakends along the chain
            // without going back through this foldback point
            if(!chain1.breakendsAreChained(var1, !v1Start, var2, !v2Start))
            {
                return;
            }
        }

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

    private void checkReplicatedBreakendFoldback(SvBreakend be)
    {
        // a special case where one ends of an SV connects to both ends of a single other variant
        // during a replication event and in doing so forms a foldback
        final SvVarData var = be.getSV();

        if(var.type() != BND || var.getReplicatedCount() != 2)
            return;

        if(!var.getAssemblyMatchType(true).equals(ASSEMBLY_MATCH_MATCHED)
        && !var.getAssemblyMatchType(false).equals(ASSEMBLY_MATCH_MATCHED))
        {
            return;
        }

        // check if the replicated SV has the same linked pairing
        // or if the variant forms both ends of the cluster's chain
        final SvCluster cluster = var.getCluster();
        if(cluster == null)
            return;

        for(final SvChain chain : cluster.getChains())
        {
            if(chain.getFirstSV().equals(var, true)
            && chain.getLastSV().equals(var, true)
            && chain.firstLinkOpenOnStart() == chain.lastLinkOpenOnStart())
            {
                var.setFoldbackLink(true, var.id(), 0);
                var.setFoldbackLink(false, var.id(), 0);

                LOGGER.debug(String.format("cluster(%s) foldback translocation SV(%s) with self",
                        cluster.getId(), var.posId()));
            }
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

    private void cacheFinalLinkedPairs(SvCluster cluster)
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
            LOGGER.debug("cluster({}: {} count={}) has {} linked pairs",
                    cluster.getId(), cluster.getDesc(), cluster.getUniqueSvCount(), linkedPairs.size());
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
