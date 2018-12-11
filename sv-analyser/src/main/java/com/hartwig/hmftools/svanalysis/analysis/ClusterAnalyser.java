package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.MIN_TEMPLATED_INSERTION_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.NO_DB_MARKER;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areSectionBreak;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_COMMON_ARMS;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_FOLDBACKS;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_SOLO_SINGLE;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.addClusterReason;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.calcConsistency;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.markLineCluster;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_ASSEMBLY_LINK_COUNT;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_LENGTH;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.svanalysis.types.SvChain.getRepeatedSvSequence;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.CLUSTER_ANNONTATION_CT;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.CLUSTER_ANNONTATION_DM;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_COMPLEX_CHAIN;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_LINE;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_NONE;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SGL_PAIR_DEL;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SGL_PAIR_DUP;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SIMPLE_CHAIN;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SGL_PAIR_INS;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SIMPLE_SV;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.areSpecificClusters;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_SGL;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.findVariantById;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.haveSameChrArms;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ClusterAnalyser {

    final SvaConfig mConfig;
    SvClusteringMethods mClusteringMethods;

    String mSampleId;
    private List<SvVarData> mAllVariants;
    List<SvCluster> mClusters;
    private ChainFinder mChainFinder;
    private LinkFinder mLinkFinder;

    public static int SMALL_CLUSTER_SIZE = 3;
    public static int SHORT_TI_LENGTH = 1000;

    private static final Logger LOGGER = LogManager.getLogger(ClusterAnalyser.class);

    public ClusterAnalyser(final SvaConfig config, SvClusteringMethods clusteringMethods)
    {
        mConfig = config;
        mClusteringMethods = clusteringMethods;
        mClusters = Lists.newArrayList();
        mAllVariants = Lists.newArrayList();
        mSampleId = "";
        mLinkFinder = new LinkFinder();
        mChainFinder = new ChainFinder();
        mChainFinder.setLogVerbose(mConfig.LogVerbose);
    }

    // access for unit testing
    public final SvaConfig getConfig() { return mConfig; }
    public final SvClusteringMethods getClusterer() { return mClusteringMethods; }

    public void setSampleData(final String sampleId, List<SvVarData> allVariants)
    {
        mSampleId = sampleId;
        mAllVariants = allVariants;
        mClusters.clear();
    }

    public final List<SvCluster> getClusters() { return mClusters; }

    public void clusterAndAnalyse()
    {
        mClusters.clear();

        mClusteringMethods.clusterByBaseDistance(mAllVariants, mClusters);

        // mark line clusters since these are exluded from most subsequent logic
        for(SvCluster cluster : mClusters)
        {
            markLineCluster(cluster, mClusteringMethods.getProximityDistance());
        }

        findSimpleCompleteChains();

        mClusteringMethods.mergeClusters(mSampleId, mClusters);

        // log basic clustering details
        for(SvCluster cluster : mClusters)
        {
            if(cluster.getCount() > 1)
                cluster.logDetails();
        }

        findLinksAndChains();

        // INVs and other SV-pairs which make foldbacks are now used in the inconsistent clustering logic
        markFoldbacks();

        mergeClusters();

        annotateTemplatedInsertions();

        // final clean-up and analysis
        for(SvCluster cluster : mClusters)
        {
            if(!cluster.isResolved() && cluster.getResolvedType() != RESOLVED_TYPE_NONE)
            {
                // any cluster with a long DEL or DUP not merged can now be marked as resolved
                if(cluster.getCount() == 1 &&  cluster.getResolvedType() == RESOLVED_TYPE_SIMPLE_SV)
                    cluster.setResolved(true, RESOLVED_TYPE_SIMPLE_SV);

                // mark off any fully chained simple clusters
                if(cluster.getCount() == 1 &&  cluster.getResolvedType() == RESOLVED_TYPE_SIMPLE_CHAIN)
                    cluster.setResolved(true, RESOLVED_TYPE_SIMPLE_CHAIN);
            }

            reportClusterFeatures(cluster);
        }

        // validation-only: checkClusterDuplicates(mClusters);
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

            if(cluster.hasLinkingLineElements())
            {
                // find assembly links but nothing else
                mLinkFinder.findLinkedPairs(mSampleId, cluster);
                cluster.cacheLinkedPairs();
                continue;
            }

            // skip more complicated clusters for now
            if(cluster.getCount() > SMALL_CLUSTER_SIZE || !cluster.isConsistent() || cluster.hasVariedCopyNumber())
                continue;

            // first establish links between SVs (eg TIs and DBs)
            mLinkFinder.findLinkedPairs(mSampleId, cluster);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster);

            cluster.cacheLinkedPairs();

            setClusterResolvedState(cluster);

            if(cluster.isFullyChained() && cluster.isConsistent())
            {
                LOGGER.debug("sample({}) cluster({}) simple and consistent with {} SVs", mSampleId, cluster.id(), cluster.getCount());
            }
        }
    }

    public void findLinksAndChains()
    {
        for (SvCluster cluster : mClusters)
        {
            // isSpecificCluster(cluster);

            if (cluster.isResolved() || cluster.isSimpleSingleSV() || cluster.isFullyChained() || cluster.getUniqueSvCount() < 2)
            {
                // no need to chain these ones
                continue;
            }

            // first establish links between SVs (eg TIs and DBs)
            mLinkFinder.findLinkedPairs(mSampleId, cluster);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster);

            cluster.cacheLinkedPairs();

            setClusterResolvedState(cluster);
            cluster.logDetails();
        }
    }

    private void mergeClusters()
    {
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
                cluster.cacheLinkedPairs();

                setClusterResolvedState(cluster);

                cluster.logDetails();
            }
        }
    }

    private void findChains(SvCluster cluster)
    {
        if(cluster.hasLinkingLineElements()) // skipped for now
            return;

        mChainFinder.initialise(mSampleId, cluster);

        isSpecificCluster(cluster);

        boolean hasFullChain = mChainFinder.formClusterChains();

        if(!hasFullChain)
        {
            checkChainReplication(cluster);
            return;
        }

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
        int maxCopyNumber = cluster.getMaxCopyNumber();

        if(maxCopyNumber > 5 * minCopyNumber)
        {
            LOGGER.debug("cluster({}) skipping replication for large CN variation(min={} max={})",
                    cluster.id(), minCopyNumber, maxCopyNumber);
            return;
        }

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
                        cluster.id(), var.posId(), svMultiple, calcCopyNumber, minCopyNumber);

                var.setReplicatedCount(svMultiple);

                // add to the parent cluster only for now
                for(int j = 1; j < svMultiple; ++j)
                {
                    SvVarData newVar = new SvVarData(var);
                    cluster.addVariant(newVar);
                }
            }

        }
    }

    private void checkChainReplication(SvCluster cluster)
    {
        if(!cluster.hasReplicatedSVs())
            return;

        // check whether any chains are identical to others using replicated SVs
        // in which case remove the replicated SVs and chains
        List<SvChain> chains = cluster.getChains();

        int index1 = 0;
        while(index1 < chains.size())
        {
            final SvChain chain1 = chains.get(index1);

            for(int index2 = index1+1; index2 < chains.size(); ++index2)
            {
                final SvChain chain2 = chains.get(index2);

                if(chain1.identicalChain(chain2))
                {
                    boolean allReplicatedSVs = chain2.getUniqueSvCount() == 0;

                    LOGGER.debug("cluster({}) removing duplicate chain({}) vs origChain({}) all replicated({})",
                            cluster.id(), chain2.id(), chain1.id(), allReplicatedSVs);

                    // remove these replicated SVs as well as the replicated chain
                    if(allReplicatedSVs)
                    {
                        for (final SvVarData var : chain2.getSvList())
                        {
                            cluster.removeReplicatedSv(var);
                        }
                    }

                    chains.remove(index2);
                    continue;

                }

                ++index2;
            }

            ++index1;
        }
    }

    private void setClusterResolvedState(SvCluster cluster)
    {
        if(!cluster.getResolvedType().equals(RESOLVED_TYPE_NONE))
            return;

        if(cluster.hasLinkingLineElements())
        {
            // skip further classification for now
            cluster.setResolved(true, RESOLVED_TYPE_LINE);
            return;
        }

        if (cluster.isSimpleSVs())
        {
            if(cluster.getTypeCount(DEL) + cluster.getTypeCount(DUP) == 2)
            {
                if(mClusteringMethods.markDelDupPairTypes(cluster))
                    return;
            }

            boolean hasLongSVs = false;
            for(SvVarData var : cluster.getSVs())
            {
                if(var.length() >= mClusteringMethods.getDelDupCutoffLength())
                {
                    hasLongSVs = true;
                    break;
                }
            }

            cluster.setResolved(!hasLongSVs, RESOLVED_TYPE_SIMPLE_SV);
            return;
        }

        // next simple reciprocal inversions and translocations
        if (cluster.getCount() == 2 && cluster.isConsistent())
        {
            if(cluster.getTypeCount(BND) == 2)
            {
                mClusteringMethods.markBndPairTypes(cluster);
            }
            else if(cluster.getTypeCount(INV) == 2)
            {
                mClusteringMethods.markInversionPairTypes(cluster);
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

        // isSpecificCluster(cluster);

        // next clusters with which start and end on the same arm, have the same start and end orientation
        // and the same start and end copy number
        if(!cluster.getChains().isEmpty())
        {
            // set the type but don't consider long chains resolved
            if(cluster.getFoldbacks().isEmpty())
                cluster.setResolved(false, RESOLVED_TYPE_SIMPLE_CHAIN);
            else
                cluster.setResolved(false, RESOLVED_TYPE_COMPLEX_CHAIN);

            return;
        }
    }

    private List<SvCluster> mergeInconsistentClusters(List<SvCluster> existingMergedClusters)
    {
        // it's possible that to resolve arms and more complex arrangements, clusters not merged
        // by proximity of overlaps must be put together to solve inconsistencies (ie loose ends)
        List<SvCluster> mergedClusters = Lists.newArrayList();

        long longDelDupCutoffLength = mClusteringMethods.getDelDupCutoffLength();

        // merge any cluster which is itself not consistent and has breakends on the same arm as another
        int index1 = 0;
        while(index1 < mClusters.size())
        {
            SvCluster cluster1 = mClusters.get(index1);

            if(cluster1.isResolved() || cluster1.hasLinkingLineElements())
            {
                ++index1;
                continue;
            }

            boolean hasOpenSingle1 = hasOneOpenSingleVariant(cluster1);
            boolean isSoloSingle1 = cluster1.getCount() == 1 && cluster1.getSVs().get(0).type() == SGL;

            boolean cluster1Merged = false;
            SvCluster newCluster = null;

            int index2 = index1 + 1;
            while(index2 < mClusters.size())
            {
                SvCluster cluster2 = mClusters.get(index2);

                if(cluster2.isResolved() || cluster2.hasLinkingLineElements())
                {
                    ++index2;
                    continue;
                }

                boolean foundConnection = false;

                foundConnection = canMergeClustersOnFoldbacks(cluster1, cluster2);

                if(!foundConnection)
                    foundConnection = canMergeClustersOnCommonArms(cluster1, cluster2, longDelDupCutoffLength);

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

                    if(!foundConnection && isSoloSingle1 && isSingleClosestTI(cluster2, cluster1))
                    {
                        foundConnection = true;
                        addClusterReason(cluster1, CLUSTER_REASON_SOLO_SINGLE, "");
                    }

                    if(!foundConnection && isSoloSingle2 && isSingleClosestTI(cluster1, cluster2))
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
                            cluster1.id(), cluster1.getCount(), cluster2.id(), cluster2.getCount());

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
                }
                else
                {
                    cluster1Merged = true;
                    cluster2Merged = true;

                    newCluster = new SvCluster(mClusteringMethods.getNextClusterId());

                    LOGGER.debug("new cluster({}) from merge of cluster({} svs={}) and cluster({} svs={})",
                            newCluster.id(), cluster1.id(), cluster1.getCount(), cluster2.id(), cluster2.getCount());

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
        final SvVarData soloVar = soloSingleCluster.getSVs().get(0);

        final List<SvBreakend> breakendList = mClusteringMethods.getChrBreakendMap().get(soloVar.chromosome(true));

        for (int i = 1; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);

            if(breakend.getSV().equals(soloVar))
            {
                // check if the preceding or next breakend is in the cluster with the open single chain
                if((i > 0 && otherCluster.getSVs().contains(breakendList.get(i - 1).getSV()))
                || (i < breakendList.size() - 1 && otherCluster.getSVs().contains(breakendList.get(i + 1).getSV())))
                {
                    LOGGER.debug("cluster({}) and cluster({}) have adjacent SVs with solo-single({})",
                            otherCluster.id(), soloSingleCluster.id(), soloVar.id());

                    return true;
                }

                break;
            }
        }

        return false;
    }

    private boolean isSingleClosestTI(final SvCluster otherCluster, final SvCluster soloSingleCluster)
    {
        // first cluster has the closest VS to the solo single and can form a TI with it
        final SvVarData soloVar = soloSingleCluster.getSVs().get(0);

        final List<SvBreakend> breakendList = mClusteringMethods.getChrBreakendMap().get(soloVar.chromosome(true));

        for (int i = 1; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);

            if(breakend.getSV().equals(soloVar))
            {
                // check if the preceding or next breakend is in the cluster with the open single chain
                final SvBreakend prevBreakend = (i > 0) ? breakendList.get(i - 1) : null;
                final SvBreakend nextBreakend = (i < breakendList.size() - 1) ? breakendList.get(i + 1) : null;
                boolean checkPrev = prevBreakend != null && otherCluster.getSVs().contains(prevBreakend.getSV());
                boolean checkNext = nextBreakend != null && otherCluster.getSVs().contains(nextBreakend.getSV());

                if((checkPrev && canFormTIWithChainEnd(otherCluster.getChains(), prevBreakend.getSV(), soloVar))
                || (checkNext && canFormTIWithChainEnd(otherCluster.getChains(), nextBreakend.getSV(), soloVar)))
                {
                    LOGGER.debug("cluster({}) can form TI with cluster({}) and soloSingle({})",
                            otherCluster.id(), soloSingleCluster.id(), soloVar.id());
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

                        if (!var1.chromosome(v1Start).equals(var2.chromosome(v2Start)) || !var1.arm(v1Start).equals(var2.arm(v2Start)))
                            continue;

                        LOGGER.debug("cluster({}) SV({}) and cluster({}) SV({}) have foldbacks on same arm",
                                cluster1.id(), var1.posId(), cluster2.id(), var2.posId());

                        addClusterReason(cluster1, CLUSTER_REASON_FOLDBACKS, var2.id());
                        addClusterReason(cluster2, CLUSTER_REASON_FOLDBACKS, var1.id());
                        return true;
                    }
                }
            }
        }

        return false;
    }

    private boolean canMergeClustersOnCommonArms(final SvCluster cluster1, final SvCluster cluster2, long armWidthCutoff)
    {
        // merge if the 2 clusters have BNDs linking the same 2 inconsistent or long arms
        areSpecificClusters(cluster1, cluster2);

        // re-check which BNDs may link arms
        cluster1.setArmLinks();
        cluster2.setArmLinks();

        // frist find arm groups which are inconsistent in both clusters
        // BNDs only touching an arm in a short TI are ignored from the common arm check
        List<SvArmGroup> inconsistentArms1 = cluster1.getArmGroups().stream()
                .filter(x -> x.canLink(armWidthCutoff))
                .collect(Collectors.toList());

        List<SvArmGroup> inconsistentArms2 = cluster2.getArmGroups().stream()
                .filter(x -> x.canLink(armWidthCutoff))
                .collect(Collectors.toList());

        // now search for common BNDs where either end is in one of these inconsistent arms
        if(inconsistentArms1.isEmpty() || inconsistentArms2.isEmpty())
            return false;

        final List<SvVarData> crossArmList1 = cluster1.getUnlinkedRemoteSVs();
        final List<SvVarData> crossArmList2 = cluster2.getUnlinkedRemoteSVs();

        // now that the candidate arm groups have been established, just need to find a single BND
        // from each cluster that falls into the same par of arm groups

        for (final SvVarData var1 : crossArmList1)
        {
            for (final SvVarData var2 : crossArmList2)
            {
                if(!haveSameChrArms(var1, var2))
                    continue;

                for(final SvArmGroup armGroup1 : inconsistentArms1)
                {
                    if (!armGroup1.getSVs().contains(var1))
                        continue;

                    for (final SvArmGroup armGroup2 : inconsistentArms2)
                    {
                        if (!armGroup2.getSVs().contains(var2))
                            continue;

                        LOGGER.debug("cluster({}) and cluster({}) have common links with SV({}) and SV({})",
                                cluster1.id(), cluster2.id(), var1.posId(), var2.posId());

                        final String commonArms = var1.id() + "_" + var2.id();

                        addClusterReason(cluster1, CLUSTER_REASON_COMMON_ARMS, commonArms);
                        addClusterReason(cluster2, CLUSTER_REASON_COMMON_ARMS, commonArms);
                        return true;
                    }
                }
            }
        }

        return false;
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
                        cluster.id(), mainChainCount, cluster.getSubClusters().size(), maxSubClusterChainCount);
            }

            mClusters.remove(cluster);
            mergedClusters.remove(index);
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

                // pass in surrounding breakends as well to check for conflictinlinks
                final SvBreakend prePrevBreakend = (i > 1) ? breakendList.get(i - 2) : null;
                final SvBreakend nextBreakend = (i < breakendList.size() - 1) ? breakendList.get(i + 1) : null;

                checkFoldbackBreakends(breakend, prevBreakend, nextBreakend, prePrevBreakend);

                checkReplicatedBreakendFoldback(breakend);
            }
        }
    }

    private void checkFoldbackBreakends(SvBreakend be1, SvBreakend be2, SvBreakend postBe1, SvBreakend preBe2)
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

        if(cluster1.isResolved() && cluster1.getTypeCount(INV) == 0)
            return;

        SvCluster cluster2 = null;

        if(var1.equals(var2))
        {
            cluster2 = cluster1;
        }
        else
        {
            cluster2 = var2.getCluster();

            if (cluster2.isResolved())
                return;
        }

        boolean v1Start = be1.usesStart();
        boolean v2Start = be2.usesStart();

        String chainInfo = "";

        if(var1.equals(var2))
        {
            // constraint is that the ends of this INV don't link to BND taking the path off this chromosome
            final SvChain chain = cluster1.findChain(var1);

            if(chain != null)
            {
                int bndLinks = 0;
                for (final SvLinkedPair pair : chain.getLinkedPairs())
                {
                    if (pair.first().equals(var1, true) && pair.second().type() == BND)
                        ++bndLinks;
                    else if (pair.second().equals(var1, true) && pair.first().type() == BND)
                        ++bndLinks;

                    if (bndLinks == 2)
                        return;
                }
            }
        }
        else
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
            int[] chainData = chain1.breakendsAreChained(var1, !v1Start, var2, !v2Start);

            if(chainData[CHAIN_LINK_COUNT] == 0)
                return;

            chainInfo = String.format("%d;%d;%d",
                    chainData[CHAIN_LINK_COUNT], chainData[CHAIN_ASSEMBLY_LINK_COUNT], chainData[CHAIN_LENGTH]);

            // check for a conflicting deletion bridge on the backmost of the 2 breakends
            if(be2.orientation() == 1)
            {
                if(preBe2 != null && areSectionBreak(be2.getSV(), preBe2.getSV(), be2.usesStart(), preBe2.usesStart()))
                {
                    // if the DB length is short, consider this a conflicting link and skip the foldback
                    if(be2.position() - preBe2.position() <= MIN_TEMPLATED_INSERTION_LENGTH)
                        return;
                }
            }
            else
            {
                if(postBe1 != null && areSectionBreak(be1.getSV(), postBe1.getSV(), be1.usesStart(), postBe1.usesStart()))
                {
                    // if the DB length is short, consider this a conflicting link and skip the foldback
                    if(postBe1.position() - be1.position() <= MIN_TEMPLATED_INSERTION_LENGTH)
                        return;
                }
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

        var1.setFoldbackLink(v1Start, var2.id(), length, chainInfo);
        var2.setFoldbackLink(v2Start, var1.id(), length, chainInfo);

        if(var1.equals(var2))
        {
            LOGGER.debug(String.format("cluster(%s) foldback inversion SV(%s) length(%d) copyNumber(%.3f)",
                    cluster1.id(), var1.posId(), length, cn1));
        }
        else
        {
            LOGGER.debug(String.format("cluster(%s) foldback be1(%s) be2(%s) length(%d) copyNumber(%.3f)",
                    cluster1.id(), be1.toString(), be2.toString(), length, cn1));
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
                final String chainInfo = String.format("%d;%d;%d",
                        chain.getLinkCount(), chain.getAssemblyLinkCount(), chain.getLength());

                var.setFoldbackLink(true, var.id(), 0, chainInfo);
                var.setFoldbackLink(false, var.id(), 0, chainInfo);

                LOGGER.debug(String.format("cluster(%s) foldback translocation SV(%s) with self",
                        cluster.id(), var.posId()));
            }
        }
    }

    private void clearFoldbackInfo(final String varId, final String matchVarId, SvCluster cluster, boolean useStart)
    {
        SvVarData var = findVariantById(varId, cluster.getSVs());

        if(var == null)
            return;

        if(var.getFoldbackLink(true).equals(matchVarId))
            var.setFoldbackLink(true, "", -1, "");
        else
            var.setFoldbackLink(false, "", -1, "");
    }

    private void annotateTemplatedInsertions()
    {
        /* work out:
            - whether a TI has a DB on either or both sides in the same cluster
            - number of chained assembled TIs in a row
            - if local, the distance to the next SV in the cluster
         */

        final Map<String, List<SvBreakend>> chrBreakendMap = mClusteringMethods.getChrBreakendMap();

        for(final SvCluster cluster : mClusters)
        {
            if(cluster.getChains().isEmpty() || cluster.hasLinkingLineElements())
                continue;

            isSpecificCluster(cluster);

            // gather up start and end arms from each chain
            List<String> startEndArms = Lists.newArrayList();

            for(final SvChain chain : cluster.getChains())
            {
                for (int be1 = SVI_START; be1 <= SVI_END; ++be1)
                {
                    boolean isFirst = isStart(be1);

                    final SvBreakend chainEnd = chain.getChainEndSV(isFirst).getBreakend(chain.chainEndOpenOnStart(isFirst));
                    if (chainEnd == null)
                        continue;

                    if (!startEndArms.contains(chainEnd.getChrArm()))
                        startEndArms.add(chainEnd.getChrArm());
                }
            }

            for(final SvChain chain : cluster.getChains())
            {
                List<SvLinkedPair> assembledLinks = Lists.newArrayList();

                for(final SvLinkedPair pair : chain.getLinkedPairs())
                {
                    if(pair.linkType() != LINK_TYPE_TI)
                        continue;

                    if(pair.first().type() == SGL || pair.second().type() == SGL)
                        continue;

                    if(pair.isAssembled())
                    {
                        assembledLinks.add(pair);
                    }
                    else if(!assembledLinks.isEmpty())
                    {
                        for(final SvLinkedPair assembledPair : assembledLinks)
                        {
                            assembledPair.setAssembledChainCount(assembledLinks.size());
                        }

                        assembledLinks.clear();
                    }

                    // find closest SV in this cluster
                    final SvVarData first = pair.first();
                    final SvVarData second = pair.second();
                    SvBreakend firstBreakend = first.getBreakend(pair.firstLinkOnStart());
                    final List<SvBreakend> breakendList = chrBreakendMap.get(firstBreakend.chromosome());

                    int[] nextSVData = getNextClusterSVData(cluster, breakendList, pair);
                    pair.setNextSVData(nextSVData[NEXT_SV_DISTANCE], nextSVData[NEXT_SV_TRAVERSED_COUNT]);

                    pair.setTraversedSVCount(getTraversedSVs(cluster, breakendList, pair));

                    SvLinkedPair dbFirst = first.getDBLink(pair.firstLinkOnStart());
                    SvLinkedPair dbSecond = pair.second().getDBLink(pair.secondLinkOnStart());

                    pair.setDBLenFirst(dbFirst != null ? dbFirst.length() : NO_DB_MARKER);
                    pair.setDBLenSecond(dbSecond != null ? dbSecond.length() : NO_DB_MARKER);

                    if(startEndArms.contains(firstBreakend.getChrArm()))
                    {
                        pair.setOnArmOfOrigin(true);
                    }

                    // skip replicated SVs since more complicate to figure out, also those crossing centromere
                    boolean hasReplicatedSVs = (first.isReplicatedSv() || first.getReplicatedCount() > 0 || pair.second().isReplicatedSv() || pair.second().getReplicatedCount() > 0);
                    boolean sameArm = first.arm(pair.firstLinkOnStart()).equals(pair.second().arm(pair.secondLinkOnStart()));

                    if(!hasReplicatedSVs && sameArm)
                    {
                        // estimate whether this TI gives copy number gain
                        double maxTelomereCN = mClusteringMethods.chromosomeCopyNumber(firstBreakend.chromosome());

                        boolean firstHasGain = first.copyNumber(pair.firstLinkOnStart()) - first.copyNumberChange(pair.firstLinkOnStart()) >= maxTelomereCN * 0.95;
                        boolean secondHasGain = second.copyNumber(pair.secondLinkOnStart()) - second.copyNumberChange(pair.secondLinkOnStart()) >= maxTelomereCN * 0.95;

                        if(firstHasGain && secondHasGain)
                        {
                            pair.setCopyNumberGain(true);
                        }
                    }
                }
            }
        }
    }

    private static int getTraversedSVs(final SvCluster cluster, final List<SvBreakend> breakendList, final SvLinkedPair pair)
    {
        // count any non-trivial cluster's SVs crossed by this pair
        final SvBreakend firstBreakend = pair.first().getBreakend(pair.firstLinkOnStart());
        final SvBreakend secondBreakend = pair.second().getBreakend(pair.secondLinkOnStart());
        int lowerIndex = min(firstBreakend.getChrPosIndex(), secondBreakend.getChrPosIndex());
        int upperIndex = max(firstBreakend.getChrPosIndex(), secondBreakend.getChrPosIndex());

        if(lowerIndex >= upperIndex - 1)
            return 0;

        int unclusteredTraversedCount = 0;

        for (int i = lowerIndex + 1; i < upperIndex - 1; ++i)
        {
            final SvCluster otherCluster = breakendList.get(i).getSV().getCluster();

            if (otherCluster == cluster || otherCluster.isResolved())
                continue;

            ++unclusteredTraversedCount;
        }

        return unclusteredTraversedCount;
    }

    private static int NEXT_SV_DISTANCE = 0;
    private static int NEXT_SV_TRAVERSED_COUNT = 1;

    private static int[] getNextClusterSVData(final SvCluster cluster, final List<SvBreakend> breakendList, final SvLinkedPair pair)
    {
        // walk forward and backwards from this pair to the closest SV in the same cluster
        // counting the number of non-trivial SVs traversed in the process
        int[] nextSvData = {-1, 0};

        SvBreakend firstBreakend = pair.first().getBreakend(pair.firstLinkOnStart());
        SvBreakend secondBreakend = pair.second().getBreakend(pair.secondLinkOnStart());
        int lowerIndex = min(firstBreakend.getChrPosIndex(), secondBreakend.getChrPosIndex());
        int upperIndex = max(firstBreakend.getChrPosIndex(), secondBreakend.getChrPosIndex());
        int svsTraversed = 0;

        if(lowerIndex > 0)
        {
            for(int i = lowerIndex - 1; i >= 0; --i)
            {
                final SvBreakend breakend = breakendList.get(i);

                if (breakend.getSV().getCluster() == cluster)
                {
                    final SvBreakend refBreakend = breakendList.get(lowerIndex);
                    nextSvData[NEXT_SV_DISTANCE] = (int)(refBreakend.position() - breakend.position());
                    nextSvData[NEXT_SV_TRAVERSED_COUNT] = svsTraversed;
                    break;
                }
                else if(!breakend.getSV().getCluster().isResolved())
                {
                    ++svsTraversed;
                }
            }
        }

        if(upperIndex < breakendList.size() - 1)
        {
            svsTraversed = 0;

            for(int i = upperIndex + 1; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);

                if (breakend.getSV().getCluster() == cluster)
                {
                    final SvBreakend refBreakend = breakendList.get(upperIndex);
                    long distance = breakend.position() - refBreakend.position();

                    if (nextSvData[NEXT_SV_DISTANCE] == -1 || distance < nextSvData[NEXT_SV_DISTANCE])
                    {
                        nextSvData[NEXT_SV_TRAVERSED_COUNT] = svsTraversed;
                        nextSvData[NEXT_SV_DISTANCE] = (int)distance;
                    }

                    break;
                }
                else if(!breakend.getSV().getCluster().isResolved())
                {
                    ++svsTraversed;
                }
            }
        }

        return nextSvData;
    }

    public void reportClusterFeatures(final SvCluster cluster)
    {
        /*
        for (final SvChain chain : cluster.getChains())
        {
            if(chain.getLinkedPairs().size() < 4)
                continue;

            findChainRepeatedSegments(cluster, chain);
        }
        */

        reportDoubleMinutes(cluster);

        classifySimpleChainedClusters(cluster);

        // analyseFoldbacks(cluster);
    }

    private void analyseFoldbacks(final SvCluster cluster)
    {
        if(cluster.getFoldbacks().isEmpty())
            return;

        if(cluster.getResolvedType() != RESOLVED_TYPE_DEL_EXT_TI && cluster.getTypeCount(INV) != 2)
            return;

        // check whether the shorter of the 2 INVs faces the centromere
        final SvVarData var1 = cluster.getSV(0);
        final SvVarData var2 = cluster.getSV(1);

        boolean shorterInvFacesCentromere = false;

        if(var1.length() < var2.length() && var1.orientation(true) == -1)
            shorterInvFacesCentromere = true;
        else if(var2.length() < var1.length() && var2.orientation(true) == -1)
            shorterInvFacesCentromere = true;

    }

    private static int CHAIN_TI_COUNT = 0;
    private static int CHAIN_TI_CN_GAIN = 1;
    private static int CHAIN_TI_DB_COUNT = 2;
    private static int CHAIN_TI_SHORT_COUNT = 3;
    private static int CHAIN_TI_ASMB_COUNT = 4;

    private void classifySimpleChainedClusters(final SvCluster cluster)
    {
        if(cluster.isResolved())
            return;

        boolean isComplex = cluster.hasReplicatedSVs() || !cluster.getFoldbacks().isEmpty();
        boolean isIncomplete = !cluster.isFullyChained() || cluster.getTypeCount(SGL) > 0;

        // skip simple chained clusters
        if(cluster.getArmCount() == 1 && cluster.getCount() == 2)
            return;

        /* data to gather for each arm in the chain
            - number of links
            - number of short TIs without proximate deletion bridges
            - links with copy number gain
            - start and end locations
         */

        long proximityCutoff = mClusteringMethods.getProximityDistance();

        boolean allChainsConsistent = true;
        List<String> originArms = Lists.newArrayList();
        List<String> fragmentArms = Lists.newArrayList();

        for (final SvChain chain : cluster.getChains())
        {
            if(!chain.isConsistent())
            {
                allChainsConsistent = false;
                continue;
            }

            Map<String, int[]> armDataMap = new HashMap();

            final SvBreakend firstBreakend = chain.getFirstSV().getBreakend(chain.firstLinkOpenOnStart());
            final SvBreakend lastBreakend = chain.getLastSV().getBreakend(chain.lastLinkOpenOnStart());

            final String startChrArm = firstBreakend != null ? firstBreakend.getChrArm() : "";
            final String endChrArm = lastBreakend != null ? lastBreakend.getChrArm() : "";

            if(!startChrArm.isEmpty())
                armDataMap.put(startChrArm, new int[CHAIN_TI_ASMB_COUNT+1]);

            if(!endChrArm.isEmpty())
                armDataMap.put(endChrArm, new int[CHAIN_TI_ASMB_COUNT+1]);

            for(final SvLinkedPair pair : chain.getLinkedPairs())
            {
                final SvVarData first = pair.first();

                if(pair.first().type() == SGL || pair.second().type() == SGL)
                    continue;

                final String chrArm = first.getBreakend(pair.firstLinkOnStart()).getChrArm();

                int[] armData = armDataMap.get(chrArm);

                if(armData == null)
                {
                    armData = new int[CHAIN_TI_ASMB_COUNT+1];
                    armDataMap.put(chrArm, armData);
                }

                ++armData[CHAIN_TI_COUNT];

                if(pair.hasCopyNumberGain())
                {
                    ++armData[CHAIN_TI_CN_GAIN];
                }

                if((pair.getDBLenFirst() > NO_DB_MARKER && pair.getDBLenFirst() <= proximityCutoff)
                || (pair.getDBLenSecond() > NO_DB_MARKER && pair.getDBLenSecond() <= proximityCutoff))
                {
                    ++armData[CHAIN_TI_DB_COUNT];
                }

                if(pair.length() <= SHORT_TI_LENGTH)
                {
                    ++armData[CHAIN_TI_SHORT_COUNT];

                    if(pair.isAssembled())
                        ++armData[CHAIN_TI_ASMB_COUNT];
                }
            }

            String chainInfo = startChrArm + "-" + endChrArm;

            for (Map.Entry<String,int[]> entry : armDataMap.entrySet())
            {
                final String chrArm = entry.getKey();
                final int[] armData = entry.getValue();

                int linkCount = armData[CHAIN_TI_COUNT];
                int cnGain = armData[CHAIN_TI_CN_GAIN];
                int dbCount = armData[CHAIN_TI_DB_COUNT];
                int shortCount = armData[CHAIN_TI_SHORT_COUNT];
                int assembledCount = armData[CHAIN_TI_ASMB_COUNT];

                boolean isOrigin = false;
                if(chrArm.equals(startChrArm) || chrArm.equals(endChrArm))
                {
                    isOrigin = true;
                }
                else if(cnGain == linkCount && shortCount > dbCount)
                {
                    isOrigin = false;
                }

                if(isOrigin)
                {
                    // if this arm exists twice already, then more than 2 chains end up on the same arm which is invalid
                    long chrArmCount = originArms.stream().filter(x -> x.equals(chrArm)).count();

                    if(chrArmCount >= 2)
                        allChainsConsistent = false;

                    originArms.add(chrArm);

                    if(fragmentArms.contains(chrArm))
                        fragmentArms.remove(chrArm);
                }
                else if(!isOrigin && !fragmentArms.contains(chrArm))
                {
                    fragmentArms.add(chrArm);
                }

                chainInfo += String.format(" %s %s: LK=%d CG=%d DB=%d SH=%d AS=%d",
                        isOrigin ? "O" : "F", chrArm, linkCount, cnGain, dbCount, shortCount, assembledCount);
            }

            chain.setDetails(chainInfo);

            LOGGER.debug("cluster({}) chain({}) {}",
                    cluster.id(), chain.id(), chainInfo);
        }

        cluster.setArmData(originArms.size(), fragmentArms.size());

        if(allChainsConsistent && !isIncomplete && !isComplex)
        {
            cluster.addAnnotation(String.format("%s", CLUSTER_ANNONTATION_CT));
        }
    }

    private static double DOUBLE_MINUTE_COPY_NUMBER_THRESHOLD = 8;

    private void reportDoubleMinutes(final SvCluster cluster)
    {
        if(cluster.isResolved() || !cluster.isFullyChained() || cluster.getChains().get(0).getLinkCount() <= 2)
            return;

        if(!cluster.getChains().get(0).couldFormLoop())
            return;

        final SvChain chain = cluster.getChains().get(0);

        // check for high copy number within this loop
        boolean hasHighCN = false;
        double cnTotal = 0;
        double maxCN = 0;
        double ploidyTotal = 0;

        for(final SvVarData var : chain.getSvList())
        {
            ploidyTotal += var.getSvData().ploidy();
            cnTotal += (var.copyNumber(true) + var.copyNumber(false)) * 0.5;
            maxCN = max(maxCN, var.copyNumber(true));
            maxCN = max(maxCN, var.copyNumber(false));

            if(var.copyNumber(true) >= DOUBLE_MINUTE_COPY_NUMBER_THRESHOLD || var.copyNumber(false) >= DOUBLE_MINUTE_COPY_NUMBER_THRESHOLD)
            {
                hasHighCN = true;
            }
        }

        if(hasHighCN)
        {
            LOGGER.debug(String.format("sample(%s) cluster(%d) chain(%d) links(%d) closed loop, copyNumber(avg=%.1f max=%.2f) ploidy(%.1f)",
                    mSampleId, cluster.id(), chain.id(), chain.getLinkCount(),
                    cnTotal/chain.getSvList().size(), maxCN, ploidyTotal/chain.getSvList().size()));

            cluster.addAnnotation(CLUSTER_ANNONTATION_DM);
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

                    if(forwardRepeats.size() >= 4)
                    {
                        LOGGER.debug("sample({}) cluster({}) chain({}) {} sequence of {} SVs starting at index({}:{}) SV({})",
                                mSampleId, cluster.id(), chain.id(), forwardSequence ? "forward" : "reverse",
                                forwardRepeats.size(), i, j, var1.id());

                        // ClusterId,ChainId,SequenceCount,VarIds,MatchDirection
                        LOGGER.debug("CF_REPEAT_SEQ: {},{},{},{},{},{}",
                                mSampleId, cluster.id(), chain.id(), forwardRepeats.size(), svIds, forwardSequence);
                    }

                    break;
                }

                // no sequence found
            }
        }
    }

    public void logStats()
    {
        mChainFinder.getContinuousFinderPc().logStats(false);
    }

}
