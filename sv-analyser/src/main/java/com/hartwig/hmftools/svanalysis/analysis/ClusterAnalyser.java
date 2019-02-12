package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.MULTIPLE;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.MIN_TEMPLATED_INSERTION_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.NO_DB_MARKER;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.areLinkedSection;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_COMMON_ARMS;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_FOLDBACKS;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_LOOSE_OVERLAP;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.MAX_SV_REPLICATION_MULTIPLE;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.addClusterReason;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.applyCopyNumberReplication;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.markLineCluster;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.logArmClusterData;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.mergeArmClusters;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_ASSEMBLY_LINK_COUNT;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_LENGTH;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.svanalysis.types.SvChain.checkChainReplication;
import static com.hartwig.hmftools.svanalysis.types.SvChain.getRepeatedSvSequence;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.CLUSTER_ANNONTATION_CT;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.CLUSTER_ANNONTATION_DM;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_COMPLEX_CHAIN;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DEL_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_DUP_EXT_TI;
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
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.svanalysis.types.SvArmCluster;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLOH;
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

    private boolean mRunValidationChecks;

    PerformanceCounter mPcClustering;
    PerformanceCounter mPcChaining;
    PerformanceCounter mPcAnnotation;

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
        mLinkFinder.setLogVerbose(mConfig.LogVerbose);
        mRunValidationChecks = false;

        mPcClustering = new PerformanceCounter("Clustering");
        mPcChaining = new PerformanceCounter("Chaining");
        mPcAnnotation = new PerformanceCounter("Annotation");
    }

    // access for unit testing
    public final SvClusteringMethods getClusterer() { return mClusteringMethods; }
    public final ChainFinder getChainFinder() { return mChainFinder; }
    public final LinkFinder getLinkFinder() { return mLinkFinder; }

    public void setRunValidationChecks(boolean toggle) { mRunValidationChecks = toggle; }

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

        mPcClustering.start();

        mClusteringMethods.clusterByProximity(mAllVariants, mClusters);

        // mark line clusters since these are exluded from most subsequent logic
        for(SvCluster cluster : mClusters)
        {
            markLineCluster(cluster, mClusteringMethods.getProximityDistance());
        }

        if(mRunValidationChecks)
        {
            if(!mClusteringMethods.validateClustering(mClusters))
            {
                LOGGER.info("exiting with cluster-validation errors");
                return;
            }
        }

        findSimpleCompleteChains();

        mClusteringMethods.mergeClusters(mSampleId, mClusters);

        mPcClustering.pause();

        // log basic clustering details
        for(SvCluster cluster : mClusters)
        {
            if(cluster.getCount() > 1)
                cluster.logDetails();
        }

        mPcChaining.start();
        findLinksAndChains();
        mPcChaining.pause();

        // INVs and other SV-pairs which make foldbacks are now used in the inconsistent clustering logic
        markFoldbacks();

        mergeClusters();

        if(mRunValidationChecks)
        {
            if(!mClusteringMethods.validateClustering(mClusters) || !mClusteringMethods.checkClusterDuplicates(mClusters))
            {
                LOGGER.info("exiting with cluster-validation errors");
                return;
            }
        }

        mPcAnnotation.start();
        // analyseOverlappingTIs();

        annotateTemplatedInsertions();
        // checkSkippedLOHEvents();

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

            mergeArmClusters(cluster.getArmClusters());

            // if(LOGGER.isDebugEnabled())
            //     logArmClusterData(cluster);

            reportClusterFeatures(cluster);
            annotateClusterArmSegments(cluster);
        }

        mPcAnnotation.stop();

        // validation-only: checkClusterDuplicates(mClusters);
    }

    public void findSimpleCompleteChains()
    {
        // for small clusters, try to find a full chain through all SVs
        List<SvCluster> simpleClusters = Lists.newArrayList();

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

            // inferred links are used to classify simple resolved types involving 2-3 SVs
            mLinkFinder.findLinkedPairs(cluster, true);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster);

            cluster.cacheLinkedPairs();
            simpleClusters.add(cluster);
        }

        for(final SvCluster cluster : simpleClusters)
        {
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

            if (cluster.isResolved() && cluster.getResolvedType() != RESOLVED_TYPE_LINE)
                continue;

            // these are either already chained or no need to chain
            if (cluster.isSimpleSingleSV() || cluster.isFullyChained() || cluster.getUniqueSvCount() < 2)
                continue;

            cluster.dissolveLinksAndChains();

            // first establish links between SVs (eg TIs and DBs)
            mLinkFinder.findLinkedPairs(cluster, false);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster);

            cluster.cacheLinkedPairs();

            setClusterResolvedState(cluster);
            cluster.logDetails();
        }
    }

    public void mergeClusters()
    {
        // now look at merging unresolved & inconsistent clusters where they share the same chromosomal arms
        List<SvCluster> mergedClusters = Lists.newArrayList();

        mPcClustering.resume();

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

        mPcClustering.stop();

        if(!mergedClusters.isEmpty())
        {
            mPcChaining.resume();

            for(SvCluster cluster : mergedClusters)
            {
                cluster.setDesc(cluster.getClusterTypesAsString());

                // NEW LOGIC:
                cluster.dissolveLinksAndChains();
                cluster.removeReplicatedSvs();

                if(cluster.hasVariedCopyNumber())
                    applyCopyNumberReplication(cluster);

                mLinkFinder.findLinkedPairs(cluster, false);

                findChains(cluster);
            }

            mPcChaining.stop();

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
        mChainFinder.initialise(cluster);

        // isSpecificCluster(cluster);

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

    private void setClusterResolvedState(SvCluster cluster)
    {
        if(cluster.getResolvedType() != RESOLVED_TYPE_NONE)
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
                final SvVarData var1 = cluster.getSVs().get(0);
                final SvVarData var2 = cluster.getSVs().get(1);

                if(var1.orientation(true) != var2.orientation(true))
                {
                    long length = abs(var1.position(true) - var2.position(true));

                    if(length < MIN_TEMPLATED_INSERTION_LENGTH)
                    {
                        cluster.setResolved(true, RESOLVED_TYPE_SGL_PAIR_INS);
                    }
                    else
                    {
                        boolean v1First = var1.position(true) < var2.position(true);
                        boolean v1PosOrientation = var1.orientation(true) == 1;

                        if(v1First == v1PosOrientation)
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
        // second round of cluster merging on more complex criteria and inconsistencies:
        // merge on foldbacks on the same arm
        // merge on links between common arms
        // merge if one cluster has footprints which overlap unresolved complex SVs

        // any merges result in a new cluster with the original clusters made into sub-clusters
        // subsequent merging keeps the 'super' clusters and adds more sub-clusters
        // the purpose of sub-clusters was to a) allow de-merging and b) preserve chains and links
        // but this could be revisited
        List<SvCluster> mergedClusters = Lists.newArrayList();

        long longDelDupCutoffLength = mClusteringMethods.getDelDupCutoffLength();

        int index1 = 0;
        while(index1 < mClusters.size())
        {
            SvCluster cluster1 = mClusters.get(index1);

            if(cluster1.isResolved())
            {
                ++index1;
                continue;
            }

            List<SvCluster> cluster1Overlaps = getTraversedClusters(cluster1);

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

                // try each merge reason in turn
                boolean canMergeClusters = false;

                canMergeClusters = canMergeClustersOnFoldbacks(cluster1, cluster2);

                if(!canMergeClusters)
                    canMergeClusters = canMergeClustersOnCommonArms(cluster1, cluster2, longDelDupCutoffLength);

                if(!canMergeClusters)
                {
                    List<SvCluster> cluster2Overlaps = getTraversedClusters(cluster2);

                    if(cluster1Overlaps.contains(cluster2) || cluster2Overlaps.contains(cluster1))
                    {
                        addClusterReason(cluster1, CLUSTER_REASON_LOOSE_OVERLAP, "");
                        addClusterReason(cluster2, CLUSTER_REASON_LOOSE_OVERLAP, "");
                        canMergeClusters = true;
                    }
                }

                if(!canMergeClusters)
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

    private static int MAX_FOLDBACK_NEXT_CLUSTER_DISTANCE = 5000000;

    private boolean canMergeClustersOnFoldbacks(final SvCluster cluster1, final SvCluster cluster2)
    {
        final List<SvVarData> cluster1Foldbacks = cluster1.getFoldbacks();
        final List<SvVarData> cluster2Foldbacks = cluster2.getFoldbacks();

        if(cluster1Foldbacks.isEmpty() && cluster2Foldbacks.isEmpty())
            return false;

        if(!cluster1Foldbacks.isEmpty() && !cluster2Foldbacks.isEmpty())
        {
            for (final SvVarData var1 : cluster1Foldbacks)
            {
                for (int be1 = SVI_START; be1 <= SVI_END; ++be1)
                {
                    boolean v1Start = isStart(be1);

                    if (be1 == SVI_END && var1.type() != BND)
                            continue;

                    if (var1.getFoldbackLink(v1Start).isEmpty())
                        continue;

                    for (final SvVarData var2 : cluster2Foldbacks)
                    {
                        for (int be2 = SVI_START; be2 <= SVI_END; ++be2)
                        {
                            boolean v2Start = isStart(be2);

                            if (be2 == SVI_END && var2.type() != BND)
                                    continue;

                            if (var2.getFoldbackLink(v2Start).isEmpty())
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
        }

        final Map<String, List<SvBreakend>> chrBreakendMap = mClusteringMethods.getChrBreakendMap();

        // additionally check whether any of the foldbacks face an opposing SV in the other cluster
        // must be the next SV, less than 5M bases and
        for(int i = 0; i <= 1; ++i)
        {
            final List<SvVarData> foldbacks = (i == 0) ? cluster1Foldbacks : cluster2Foldbacks;
            final SvCluster foldbackCluster = (i == 0) ? cluster1 : cluster2;
            final SvCluster otherCluster = (i == 0) ? cluster2 : cluster1;

            for (final SvVarData var : foldbacks)
            {
                for (int be1 = SVI_START; be1 <= SVI_END; ++be1)
                {
                    boolean isStart = isStart(be1);

                    if (be1 == SVI_END && var.type() != BND)
                        continue;

                    if (var.getFoldbackLink(isStart).isEmpty())
                        continue;

                    final SvBreakend foldbackBreakend = var.getBreakend(isStart);

                    final List<SvBreakend> breakendList = chrBreakendMap.get(foldbackBreakend.chromosome());
                    boolean traverseUp = foldbackBreakend.orientation() == -1;
                    int nextIndex = traverseUp ? foldbackBreakend.getChrPosIndex() + 1 : foldbackBreakend.getChrPosIndex() - 1;

                    SvBreakend nextBreakend = getNextUnresolvedBreakend(breakendList, nextIndex, traverseUp);

                    if(nextBreakend == null || nextBreakend.getSV().getCluster() != otherCluster)
                        continue;

                    if(abs(nextBreakend.position() - foldbackBreakend.position()) > MAX_FOLDBACK_NEXT_CLUSTER_DISTANCE)
                        continue;

                    double fbPloidy = foldbackBreakend.getSV().getSvData().ploidy();
                    double nbPloidy = nextBreakend.getSV().getSvData().ploidy();

                    if(nbPloidy < fbPloidy && !copyNumbersEqual(nbPloidy, fbPloidy))
                        continue;

                    LOGGER.debug("cluster({}) foldback breakend({}) faces cluster({}) breakend({})",
                            foldbackCluster.id(), foldbackBreakend.toString(), otherCluster.id(), nextBreakend.toString());

                    return true;
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

        // first find arm groups which are inconsistent in both clusters
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

    private List<SvCluster> getTraversedClusters(final SvCluster cluster)
    {
        // find all clusters which are overlapped by consecutive same-orientation SVs in this cluster
        // and where the overlapped clusters contain opposing orientation BNDs or SGLs
        List<SvCluster> traversedClusters = Lists.newArrayList();

        if(cluster.isResolved() || cluster.isFullyChained())
            return traversedClusters;

        final Map<String, List<SvBreakend>> chrBreakendMap = cluster.getChrBreakendMap();

        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            List<SvBreakend> fullBreakendList = mClusteringMethods.getChrBreakendMap().get(entry.getKey());

            for (int i = 0; i < breakendList.size() - 1; ++i)
            {
                final SvBreakend lowerBreakend = breakendList.get(i);
                SvBreakend upperBreakend = breakendList.get(i + 1);

                boolean isFoldback = false;

                if (lowerBreakend.arm() != upperBreakend.arm())
                    continue;

                if (lowerBreakend.orientation() != upperBreakend.orientation())
                {
                    // allow for short DBs where the breakends remain in a foldback
                    if (i < breakendList.size() - 2)
                    {
                        final SvBreakend nextBreakend = breakendList.get(i + 2);

                        if (!lowerBreakend.getSV().getFoldbackLink(lowerBreakend.usesStart()).isEmpty()
                                && lowerBreakend.getSV().getFoldbackLink(lowerBreakend.usesStart())
                                .equals(nextBreakend.getSV().getFoldbackLink(nextBreakend.usesStart())))
                        {
                            isFoldback = true;
                            upperBreakend = nextBreakend;
                        }
                    }

                    if (!isFoldback)
                        continue;
                }

                final SvBreakend frontBE = lowerBreakend.orientation() == 1 ? lowerBreakend : upperBreakend;

                if (frontBE.getSV().getDBLink(frontBE.usesStart()) != null)
                {
                    // check for an overlapping short DB which would invalidate these consecutive breakends
                    if (frontBE.getSV().getDBLink(frontBE.usesStart()).length() < 0)
                        continue;
                }

                for (int j = lowerBreakend.getChrPosIndex() + 1; j <= upperBreakend.getChrPosIndex() - 1; ++j)
                {
                    final SvBreakend breakend = fullBreakendList.get(j);

                    if(breakend.getSV().type() != BND && breakend.getSV().type() != SGL)
                        continue;

                    if(breakend.orientation() == lowerBreakend.orientation())
                        continue;

                    final SvCluster otherCluster = breakend.getSV().getCluster();

                    if (otherCluster == cluster || otherCluster.isResolved())
                        continue;

                    if(!traversedClusters.contains(otherCluster))
                    {
                        LOGGER.debug("cluster({}) breakends({} & {}) overlap cluster({}) breakend({})",
                                cluster.id(), lowerBreakend.toString(), upperBreakend.toString(), otherCluster.id(), breakend.toString());
                        traversedClusters.add(otherCluster);
                    }
                }
            }
        }

        return traversedClusters;
    }

    public void markFoldbacks()
    {
        for(final Map.Entry<String, List<SvBreakend>> entry : mClusteringMethods.getChrBreakendMap().entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            for(int i = 0; i < breakendList.size() - 1; ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvBreakend nextBreakend = breakendList.get(i + 1);

                SvBreakend beFront = null; // the lower position for orientation +1 and vice versa
                SvBreakend beBack = null;

                // isSpecificSV(breakend.getSV().id());

                if(breakend.orientation() == nextBreakend.orientation())
                {
                    beFront = breakend.orientation() == 1 ? breakend : nextBreakend;
                    beBack = breakend.orientation() == 1 ? nextBreakend : breakend;
                }
                else if(i < breakendList.size() - 2)
                {
                    SvBreakend postNextBreakend = breakendList.get(i + 2);

                    // check for a short overlapping deletion bridge on the either breakend
                    // that would otherwise mask this foldback
                    if(breakend.orientation() == postNextBreakend.orientation())
                    {
                        if(postNextBreakend.orientation() == 1
                        && postNextBreakend.position() - nextBreakend.position() < MIN_TEMPLATED_INSERTION_LENGTH)
                        {
                            beFront = breakend;
                            beBack = postNextBreakend;
                        }
                        else if(postNextBreakend.orientation() == -1
                        && nextBreakend.position() - breakend.position() < MIN_TEMPLATED_INSERTION_LENGTH)
                        {
                            beBack = breakend;
                            beFront = postNextBreakend;
                        }
                    }
                }

                if(beFront != null && beBack != null)
                {
                    // the foldback is invalid if it has a deletion bridge (including < 30b overhang) on the front-facing breakend
                    final SvLinkedPair dbLink = beFront.getSV().getDBLink(beFront.usesStart());
                    if(dbLink == null || dbLink.length() >= MIN_TEMPLATED_INSERTION_LENGTH)
                    {
                        checkFoldbackBreakends(beFront, beBack);
                    }
                }

                checkReplicatedBreakendFoldback(breakend);
            }
        }
    }

    public static int MAX_FOLDBACK_CHAIN_LENGTH = 5000;

    private void checkFoldbackBreakends(SvBreakend beStart, SvBreakend beEnd)
    {
        // SVs are marked as being in a foldback if they are consecutive breakends,
        // have the same orientation, and are either an INV or part of a chain

        // beStart is the one with the lower position

        SvVarData varEnd = beEnd.getSV();
        SvVarData varStart = beStart.getSV();

        if(varEnd.type() == INS || varStart.type() == INS)
            return;

        // skip unclustered DELs & DUPs, reciprocal INV or reciprocal BNDs
        final SvCluster cluster1 = varEnd.getCluster();

        if(cluster1.isResolved() && cluster1.getTypeCount(INV) == 0)
            return;

        SvCluster cluster2 = null;

        if(varEnd.equals(varStart))
        {
            cluster2 = cluster1;
        }
        else
        {
            cluster2 = varStart.getCluster();

            // must be same cluster
            if(cluster1 != cluster2)
                return;
        }

        boolean v1Start = beEnd.usesStart();
        boolean v2Start = beStart.usesStart();

        String chainInfo = "";

        if(varEnd.equals(varStart))
        {
            // constraint is that the ends of this INV don't link to BND taking the path off this chromosome
            final SvChain chain = cluster1.findChain(varEnd);

            if(chain != null)
            {
                int bndLinks = 0;
                for (final SvLinkedPair pair : chain.getLinkedPairs())
                {
                    if (pair.first().equals(varEnd, true) && pair.second().type() == BND)
                        ++bndLinks;
                    else if (pair.second().equals(varEnd, true) && pair.first().type() == BND)
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

            if(varEnd.getReplicatedCount() != varStart.getReplicatedCount())
                return;

            final SvChain chain1 = cluster1.findChain(varEnd);
            final SvChain chain2 = cluster2.findChain(varStart);

            if(chain1 == null || chain2 == null || chain1 != chain2)
                return;

            // check if a path can be walked between these 2 breakends along the chain
            // without going back through this foldback point
            int[] chainData = chain1.breakendsAreChained(varEnd, !v1Start, varStart, !v2Start);

            if(chainData[CHAIN_LINK_COUNT] == 0)
                return;

            int chainLength = chainData[CHAIN_LENGTH];

            if(chainLength > MAX_FOLDBACK_CHAIN_LENGTH)
                return;

            chainInfo = String.format("%d;%d;%d",
                    chainData[CHAIN_LINK_COUNT], chainData[CHAIN_ASSEMBLY_LINK_COUNT], chainLength);
        }

        // check copy numbers match
        double cn1 = 0;
        double cn2 = 0;
        if((beEnd.orientation() == 1 && beEnd.position() < beStart.position()) || (beEnd.orientation() == -1 && beEnd.position() > beStart.position()))
        {
            // be1 is facing away from be2, so need to take its copy number on the break side
            cn1 = varEnd.copyNumber(v1Start) - varEnd.copyNumberChange(v1Start);
            cn2 = varStart.copyNumber(v2Start);
        }
        else
        {
            cn2 = varStart.copyNumber(v2Start) - varStart.copyNumberChange(v2Start);
            cn1 = varEnd.copyNumber(v1Start);
        }

        if(!copyNumbersEqual(cn1, cn2))
            return;

        int length = (int)abs(beEnd.position() - beStart.position());

        // if either variant already has foldback info set, favour
        // a) simple inversions then
        // b) shortest length

        boolean skipFoldback = false;
        if(!varEnd.getFoldbackLink(v1Start).isEmpty() && !varEnd.equals(varStart) && varEnd.getFoldbackLen(v1Start) < length)
        {
            skipFoldback = true;
        }
        else if(!varStart.getFoldbackLink(v2Start).isEmpty() && !varEnd.equals(varStart) && varStart.getFoldbackLen(v2Start) < length)
        {
            skipFoldback = true;
        }

        if(skipFoldback)
            return;

        if(!varEnd.getFoldbackLink(v1Start).isEmpty())
            clearFoldbackInfo(varEnd.getFoldbackLink(v1Start), varEnd.id(), cluster1);

        if(!varStart.getFoldbackLink(v2Start).isEmpty())
            clearFoldbackInfo(varStart.getFoldbackLink(v2Start), varStart.id(), cluster2);

        varEnd.setFoldbackLink(v1Start, varStart.id(), length, chainInfo);
        varStart.setFoldbackLink(v2Start, varEnd.id(), length, chainInfo);

        if(varEnd.equals(varStart))
        {
            LOGGER.debug(String.format("cluster(%s) foldback inversion SV(%s) length(%d) copyNumber(%.3f)",
                    cluster1.id(), varEnd.posId(), length, cn1));
        }
        else
        {
            LOGGER.debug(String.format("cluster(%s) foldback be1(%s) be2(%s) length(%d) copyNumber(%.3f)",
                    cluster1.id(), beEnd.toString(), beStart.toString(), length, cn1));
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

                boolean foldbackIsStart = chain.firstLinkOpenOnStart();

                var.setFoldbackLink(foldbackIsStart, var.id(), 0, chainInfo);

                LOGGER.debug("cluster({}) foldback translocation SV({}) with self on {}",
                        cluster.id(), var.posId(), foldbackIsStart ? "start" : "end");
            }
        }
    }

    private void clearFoldbackInfo(final String varId, final String matchVarId, SvCluster cluster)
    {
        SvVarData var = findVariantById(varId, cluster.getSVs());

        if(var == null)
            return;

        if(var.getFoldbackLink(true).equals(matchVarId))
            var.setFoldbackLink(true, "", -1, "");
        else
            var.setFoldbackLink(false, "", -1, "");
    }

    private void annotateClusterArmSegments(final SvCluster cluster)
    {
        if (cluster.isResolved())
            return;

        if (cluster.getCount() < 6)
            return;

        // isSpecificCluster(cluster);

        final Map<String, List<SvBreakend>> chrBreakendMap = cluster.getChrBreakendMap();

        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            if (breakendList.size() < 6)
                continue;

            // first establish the lowest copy number segments
            int bndCount = 0;
            boolean hasConsecutiveOrientations = false;
            double lowestCopyNumber = -1;

            List<SvArmCluster> localSegments = Lists.newArrayList();
            SvArmCluster armCuster = null;

            for (int i = 0; i < breakendList.size() - 1; ++i)
            {
                final SvBreakend breakend = breakendList.get(i);

                if(!hasConsecutiveOrientations && i < breakendList.size() - 1)
                {
                    hasConsecutiveOrientations = (breakend.orientation() != breakendList.get(i+1).orientation());
                }

                double lowCopyNumber = breakend.getCopyNumber(false);

                if (breakend.orientation() == 1 && (lowestCopyNumber == -1 || lowCopyNumber < lowestCopyNumber))
                {
                    // a new lowest segment, which will be the first original DSB
                    // all previous breakends will be put into the first arm cluster
                    lowestCopyNumber = lowCopyNumber;

                    localSegments.clear();

                    // add all previous breakends to the first arm cluster
                    armCuster = new SvArmCluster(localSegments.size(), cluster, breakend.chromosome(), breakend.arm());

                    for (int j = 0; j <= i; ++j)
                        armCuster.addBreakend(breakendList.get(j));

                    localSegments.add(armCuster);
                    armCuster = null;
                }
                else if (breakend.orientation() == 1 && copyNumbersEqual(lowestCopyNumber, lowCopyNumber))
                {
                    // end of a segment since next CN equals the low - add the last breakend
                    if(armCuster == null)
                    {
                        // probably an indication of failed clustering, so work around it
                        armCuster = new SvArmCluster(localSegments.size(), cluster, breakend.chromosome(), breakend.arm());
                        localSegments.add(armCuster);
                    }

                    armCuster.addBreakend(breakend);
                    armCuster = null;
                }
                else
                {
                    // either continuation of an existing segment or start of a new one
                    if (armCuster == null)
                    {
                        armCuster = new SvArmCluster(localSegments.size(), cluster, breakend.chromosome(), breakend.arm());
                        localSegments.add(armCuster);
                    }

                    armCuster.addBreakend(breakend);
                }

                if (breakend.getSV().type() == BND)
                    ++bndCount;
            }

            if (bndCount > 0 || !hasConsecutiveOrientations)
                continue;

            for (final SvArmCluster armCluster : localSegments)
            {
                List<SvBreakend> acBreakendList = armCluster.getBreakends();

                if(acBreakendList.size() < 4)
                    continue;

                int backwardCount = (int)acBreakendList.stream().filter(x -> x.orientation() == 1).count();
                int forwardCount = (int)acBreakendList.stream().filter(x -> x.orientation() == -1).count();

                // form into mutually exclusive pairs based on copy number change match
                List<SvBreakend> unlinkedBreakends = Lists.newArrayList(acBreakendList);
                List<SvLinkedPair> pairs = formPossibleLinkedPairsByShortest(unlinkedBreakends);

                LOGGER.debug("cluster({}) armCluster({} : {}_{}) count({} fwd={} bak={}) pairs({}) unlinked({})",
                        cluster.id(), armCluster.id(), armCluster.chromosome(), armCluster.arm(),
                        acBreakendList.size(), forwardCount, backwardCount, pairs.size(), unlinkedBreakends.size());
            }
        }
    }

    private List<SvLinkedPair> formPossibleLinkedPairsByShortest(List<SvBreakend> unlinkedBreakends)
    {
        List<SvLinkedPair> pairs = Lists.newArrayList();

        while(!unlinkedBreakends.isEmpty())
        {
            boolean pairFound = false;

            for(int i = 0; i < unlinkedBreakends.size() - 2; ++i)
            {
                final SvBreakend be1 = unlinkedBreakends.get(i);
                final SvBreakend be2 = unlinkedBreakends.get(i+1);

                if(be1.orientation() != -1 || be2.orientation() != 1)
                    continue;

                // pair off these if the CN matches
                if (!copyNumbersEqual(be1.copyNumberChange(), be2.copyNumberChange()))
                    continue;

                SvLinkedPair pair = SvLinkedPair.from(be1, be2, LINK_TYPE_TI);
                pairs.add(pair);

                pairFound = true;
                unlinkedBreakends.remove(i+1); // higher index removed first
                unlinkedBreakends.remove(i);
                break;
            }

            if(!pairFound)
                break;
        }

        return pairs;
    }

    private List<SvLinkedPair> formPossibleLinkedPairsConsecutively(List<SvBreakend> unlinkedBreakends)
    {
        List<SvLinkedPair> pairs = Lists.newArrayList();

        int j = 0;
        while(j < unlinkedBreakends.size())
        {
            final SvBreakend be1 = unlinkedBreakends.get(j);

            if(be1.orientation() == 1)
            {
                ++j;
                continue;
            }

            boolean linkFound = false;
            for(int k = j+1; k < unlinkedBreakends.size(); ++k)
            {
                final SvBreakend be2 = unlinkedBreakends.get(k);

                if (be2.orientation() == -1)
                    continue;

                // pair off these if the CN matches
                if (!copyNumbersEqual(be1.copyNumberChange(), be2.copyNumberChange()))
                    continue;

                SvLinkedPair pair = SvLinkedPair.from(be1, be2, LINK_TYPE_TI);
                pairs.add(pair);

                linkFound = true;
                unlinkedBreakends.remove(k); // higher index removed first
                unlinkedBreakends.remove(j);
                break;

            }

            if(linkFound)
                continue;
            else
                ++j;
        }

        return pairs;
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
            if(cluster.getChains().isEmpty())
                continue;

            // isSpecificCluster(cluster);

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

                    pair.setTraversedSVCount(getTraversedSvCount(cluster, breakendList,
                            pair.getBreakend(true).getChrPosIndex(), pair.getBreakend(false).getChrPosIndex()));

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

    public void checkSkippedLOHEvents()
    {
        List<SvLOH> lohList = mClusteringMethods.getSampleLohData().get(mSampleId);
        List<SvLOH> unmatchedLohList = Lists.newArrayList();

        if(lohList != null)
            unmatchedLohList.addAll(lohList.stream().filter(x -> x.Skipped).collect(Collectors.toList()));

        int matchedLohCount = 0;

        // check if an LOH was a skipped for being a potential TI or DB
        int index = 0;
        while(index < unmatchedLohList.size())
        {
            final SvLOH lohEvent = unmatchedLohList.get(index);

            boolean matched = false;
            long lohLength = lohEvent.PosEnd - lohEvent.PosStart;

            final List<SvBreakend> breakendList = mClusteringMethods.getChrBreakendMap().get(lohEvent.Chromosome);

            for(int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvVarData var = breakend.getSV();

                if(!lohEvent.StartSV.equals(var.id()))
                    continue;

                if(lohEvent.StartSV.equals(var.id()) && lohEvent.EndSV.equals(var.id()))
                {
                    LOGGER.debug("var({} {}) matches skipped LOH: chr({}) breaks({} -> {}, len={})",
                            var.id(), var.type(), lohEvent.Chromosome, lohEvent.PosStart, lohEvent.PosEnd, lohLength);

                    if(var.type() == INV || var.type() == DUP)
                        matched = true;

                    break;
                }

                for (int be1 = SVI_START; be1 <= SVI_END; ++be1)
                {
                    boolean v1Start = isStart(be1);

                    final SvLinkedPair dbPair = var.getDBLink(v1Start);

                    if (dbPair != null && dbPair.getOtherSV(var).id().equals(lohEvent.EndSV)
                    && dbPair.getBreakend(true).position() == lohEvent.PosStart
                    && dbPair.getBreakend(false).position() == lohEvent.PosEnd - 1)
                    {
                        LOGGER.debug("deletionBridge({}) matches skipped LOH: chr({}) breaks({} -> {}, len={})",
                                var.getDBLink(v1Start).toString(), lohEvent.Chromosome,
                                lohEvent.PosStart, lohEvent.PosEnd, lohLength);
                        matched = true;
                        break;
                    }

                    final SvLinkedPair tiPair = var.getLinkedPair(v1Start);

                    if (tiPair != null && tiPair.getOtherSV(var).id().equals(lohEvent.EndSV)
                    && tiPair.getBreakend(true).position() == lohEvent.PosStart
                    && tiPair.getBreakend(false).position() == lohEvent.PosEnd - 1)
                    {
                        LOGGER.debug("templatedInsertion({}) matches skipped LOH: chr({}) breaks({} -> {}, len={})",
                                var.getLinkedPair(v1Start).toString(), lohEvent.Chromosome,
                                lohEvent.PosStart, lohEvent.PosEnd, lohLength);
                        matched = true;
                        break;
                    }
                }

                if(!matched)
                {
                    // check for line and SGLs which may not have formed TIs
                    SvVarData varEnd = null;

                    if (i < breakendList.size() - 1 && breakendList.get(i + 1).getSV().id().equals(lohEvent.EndSV))
                    {
                        // should be the next SV
                        varEnd = breakendList.get(i + 1).getSV();
                    }

                    if (var.inLineElement() || (varEnd != null && varEnd.inLineElement()))
                    {
                        LOGGER.debug("line SVs({} and {}) match skipped LOH: chr({}) breaks({} -> {}, len={})",
                                var.id(), varEnd != null ? varEnd.id() : "null",
                                lohEvent.Chromosome, lohEvent.PosStart, lohEvent.PosEnd, lohLength);
                        matched = true;
                    }
                    else if (var.type() == SGL || (varEnd != null && varEnd.type() == SGL))
                    {
                        matched = true;
                    }
                }

                break;
            }

            if(matched || lohEvent.matchesSegment(MULTIPLE, true) || lohEvent.matchesSegment(MULTIPLE, false))
            {
                unmatchedLohList.remove(index);
                ++matchedLohCount;
            }
            else
            {
                ++index;
            }
        }

        if(!unmatchedLohList.isEmpty())
        {
            LOGGER.info("sample({}) has matched({}) unmatched({}) skipped LOH events",
                    mSampleId, matchedLohCount, unmatchedLohList.size());

            for(final SvLOH lohEvent : unmatchedLohList)
            {
                LOGGER.info("unmatched LOH: chr({}) breaks({} -> {}, len={}) SV start({} {}) end({} {}) {} SV",
                        lohEvent.Chromosome, lohEvent.PosStart, lohEvent.PosEnd, lohEvent.PosEnd - lohEvent.PosStart,
                        lohEvent.StartSV, lohEvent.SegStart, lohEvent.EndSV, lohEvent.SegEnd,
                        lohEvent.StartSV == lohEvent.EndSV ? "same" : "diff");
            }
        }

    }

    private static int getTraversedSvCount(final SvCluster cluster, final List<SvBreakend> breakendList, int lowerIndex, int upperIndex)
    {
        String traversedInfo = getTraversedSvData(cluster, breakendList, lowerIndex, upperIndex);

        if(traversedInfo.isEmpty())
            return 0;

        String[] items = traversedInfo.split(";");
        return items.length;
    }

    private static String getTraversedSvData(final SvCluster cluster, final List<SvBreakend> breakendList, int lowerIndex, int upperIndex)
    {
        if(lowerIndex >= upperIndex - 1)
            return "";

        String traversedInfo = "";

        for (int i = lowerIndex + 1; i <= upperIndex - 1; ++i)
        {
            final SvBreakend breakend = breakendList.get(i);
            final SvCluster otherCluster = breakend.getSV().getCluster();

            if (otherCluster == cluster || otherCluster.isResolved())
                continue;

            if(!traversedInfo.isEmpty())
                traversedInfo += ";";

            traversedInfo += String.format("%d %.2f", breakend.orientation(), breakend.getSV().copyNumberChange(breakend.usesStart()));
        }

        return traversedInfo;
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

        checkLooseFoldbacks(cluster);

        classifyChainedClusters(cluster);
    }

    private final SvBreakend getNextUnresolvedBreakend(final List<SvBreakend> breakendList, int startIndex, boolean traverseUp)
    {
        int index = startIndex;

        while(index >= 0 && index < breakendList.size())
        {
            final SvBreakend breakend = breakendList.get(index);
            final SvCluster cluster = breakend.getSV().getCluster();

            if (!cluster.isResolved())
                return breakend;

            if(traverseUp)
                ++index;
            else
                --index;
        }

        return null;
    }

    private void checkLooseFoldbacks(final SvCluster cluster)
    {
        if (cluster.isResolved() || cluster.isFullyChained())
            return;

        final Map<String, List<SvBreakend>> chrBreakendMap = cluster.getChrBreakendMap();

        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            for (int i = 0; i < breakendList.size() - 1; ++i)
            {
                final SvBreakend lowerBreakend = breakendList.get(i);
                SvBreakend upperBreakend = breakendList.get(i + 1);

                boolean isFoldback = false;

                if (lowerBreakend.arm() != upperBreakend.arm())
                    continue;

                if (lowerBreakend.orientation() != upperBreakend.orientation())
                {
                    // allow for short DBs where the breakends remain in a foldback
                    if (i < breakendList.size() - 2)
                    {
                        final SvBreakend nextBreakend = breakendList.get(i + 2);

                        if (!lowerBreakend.getSV().getFoldbackLink(lowerBreakend.usesStart()).isEmpty()
                                && lowerBreakend.getSV().getFoldbackLink(lowerBreakend.usesStart())
                                .equals(nextBreakend.getSV().getFoldbackLink(nextBreakend.usesStart())))
                        {
                            isFoldback = true;
                            upperBreakend = nextBreakend;
                        }
                    }

                    if (!isFoldback)
                        continue;
                }

                final SvBreakend frontBE = lowerBreakend.orientation() == 1 ? lowerBreakend : upperBreakend;
                final SvBreakend backBE = lowerBreakend.orientation() == 1 ? upperBreakend : lowerBreakend;

                if (frontBE.getSV().getDBLink(frontBE.usesStart()) != null)
                {
                    // check for an overlapping short DB which would invalidate these consecutive breakends
                    if (frontBE.getSV().getDBLink(frontBE.usesStart()).length() < 0)
                        continue;
                }

                final SvVarData frontSv = frontBE.getSV();
                final SvVarData backSv = backBE.getSV();

                frontSv.setConsecBEStart(backSv.origId(), frontBE.usesStart());
                backSv.setConsecBEStart(frontSv.origId(), backBE.usesStart());
            }
        }
    }

    private void analyseOverlappingTIs()
    {
        // look for any TIs where one falls within the bounds of the other, as a sign of replication after shattering
        for(final Map.Entry<String, List<SvBreakend>> entry : mClusteringMethods.getChrBreakendMap().entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            for(int i = 1; i < breakendList.size() - 3; ++i)
            {
                final SvBreakend beStart1 = breakendList.get(i);
                final SvBreakend beStart2 = breakendList.get(i + 1);

                if(beStart1.orientation() == 1 || beStart2.orientation() == 1)
                    continue;

                final SvBreakend beEnd1 = breakendList.get(i + 2);
                final SvBreakend beEnd2 = breakendList.get(i + 3);

                if(beEnd1.orientation() == -1 || beEnd2.orientation() == -1)
                    continue;

                // now the orientations are correct, check if these 2 sets of breakends form one TI enclosed by the other

                SvLinkedPair pair1 = beStart1.getSV().getLinkedPair(beStart1.usesStart());
                SvLinkedPair pair2 = beStart2.getSV().getLinkedPair(beStart2.usesStart());

                boolean pair1Matched = pair1 != null && pair1.hasBreakend(beEnd2);
                boolean pair2Matched = pair2 != null && pair2.hasBreakend(beEnd1);

                boolean sameCluster = beStart1.getSV().getCluster() == beStart2.getSV().getCluster()
                        && beStart1.getSV().getCluster() == beEnd1.getSV().getCluster()
                        && beStart1.getSV().getCluster() == beEnd2.getSV().getCluster();

                long startGap = beStart2.position() - beStart1.position();

                if((pair1Matched && pair2Matched) || (pair1 != null && pair2 != null))
                {
                    LOGGER.info("sample({}) cluster({}) matched TIs: outer({}) inner({}) gap(start={} middle={} end={})",
                            mSampleId, sameCluster ? beStart1.getSV().getCluster().id() : "multiple",
                            pair1.toString(), pair2.toString(), beStart2.position() - beStart1.position(),
                            beEnd2.position() - beStart2.position(), beEnd2.position() - beEnd1.position());

                    long endGap = pair1.length() - pair2.length() - startGap;

                    if(!pair1Matched || !pair2Matched)
                    {
                        if (endGap < 0)
                            continue;
                    }

                    LOGGER.info("sample({}) cluster({}) enclosed TIs: outer({}) inner({}) gap(start={} middle={} end={})",
                            mSampleId, sameCluster ? beStart1.getSV().getCluster().id() : "multiple",
                            pair1.toString(), pair2.toString(), startGap, pair2.length(), endGap);

                    LOGGER.info("ENCLOSED_TI: {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                            mSampleId, sameCluster ? beStart1.getSV().getCluster().id() : "multiple",
                            pair1.first().chromosome(pair1.firstLinkOnStart()), pair1Matched, pair1Matched,
                            pair1.first().id(), pair1.first().position(pair1.firstLinkOnStart()),
                            pair1.second().id(), pair1.second().position(pair1.secondLinkOnStart()),
                            pair2.first().id(), pair2.first().position(pair2.firstLinkOnStart()),
                            pair2.second().id(), pair2.second().position(pair2.secondLinkOnStart()),
                            pair1.length(), pair2.length(), startGap, endGap);

                }
                else
                {
                    // breakend form 2 TIs with one enclosed by pairings don't match for some reason
                    LOGGER.info("ENCLOSED_TI: {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                            mSampleId, sameCluster ? beStart1.getSV().getCluster().id() : "multiple",
                            beStart1.chromosome(), pair1Matched, pair1Matched,
                            beStart1.getSV().id(), beStart1.position(), beStart2.getSV().id(), beStart2.position(),
                            beEnd1.getSV().id(), beEnd1.position(), beEnd2.getSV().id(), beEnd2.position(),
                            beEnd2.position() - beStart1.position(), beEnd1.position() - beStart2.position(),
                            startGap, beEnd2.position() - beEnd1.position());
                }
            }
        }
    }

    private static int CHAIN_TI_COUNT = 0;
    private static int CHAIN_TI_CN_GAIN = 1;
    private static int CHAIN_TI_DB_COUNT = 2;
    private static int CHAIN_TI_SHORT_COUNT = 3;
    private static int CHAIN_TI_ASMB_COUNT = 4;

    private void classifyChainedClusters(final SvCluster cluster)
    {
        if(cluster.isResolved())
            return;

        boolean isComplex = cluster.hasReplicatedSVs() || !cluster.getFoldbacks().isEmpty();
        boolean isIncomplete = !cluster.isFullyChained() || cluster.getTypeCount(SGL) > 0;

        // skip simple chained clusters
        if(cluster.getArmCount() == 1 && cluster.getCount() == 2)
            return;

        // isSpecificCluster(cluster);

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

            int shortTICount = 0;
            long chainLinkLength = 0;

            for(final SvLinkedPair pair : chain.getLinkedPairs())
            {
                final SvVarData first = pair.first();

                if(pair.first().type() == SGL || pair.second().type() == SGL)
                    continue;

                chainLinkLength += pair.length();

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
                    ++shortTICount;

                    if(pair.isAssembled())
                        ++armData[CHAIN_TI_ASMB_COUNT];
                }
            }

            // check for synthetic DELs and DUPs from longer chains
            if(allChainsConsistent && !isIncomplete && !isComplex
            && cluster.getChains().size() == 1 && chain.getLinkCount() == shortTICount
            && firstBreakend.getChrArm().equals(lastBreakend.getChrArm())
            && firstBreakend.orientation() != lastBreakend.orientation())
            {
                long syntheticLength = abs(lastBreakend.position() - firstBreakend.position());
                long avgLinkLength = round(chainLinkLength/chain.getLinkCount());

                cluster.setSynDelDupData(syntheticLength, avgLinkLength);

                if((firstBreakend.position() < lastBreakend.position() && firstBreakend.orientation() == 1)
                || (lastBreakend.position() < firstBreakend.position() && lastBreakend.orientation() == 1))
                {
                    cluster.setResolved(false, RESOLVED_TYPE_DEL_EXT_TI);
                }
                else
                {
                    cluster.setResolved(false, RESOLVED_TYPE_DUP_EXT_TI);
                }

                LOGGER.debug("cluster({}) chainLinks({}) synLen({}) avgTILen({}) marked as {}",
                        cluster.id(), chain.getLinkCount(), syntheticLength, avgLinkLength, cluster.getResolvedType());
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

            // LOGGER.debug("cluster({}) chain({}) {}", cluster.id(), chain.id(), chainInfo);
        }

        cluster.setArmData(originArms.size(), fragmentArms.size());

        if(allChainsConsistent && !isIncomplete && !isComplex)
        {
            cluster.addAnnotation(String.format("%s", CLUSTER_ANNONTATION_CT));
        }
    }

    private void reportClusterNeoChromosomes(final SvCluster cluster)
    {
        // count up number of arms from the breakends, excluding those in a TI
        int unlinkedSvCount = cluster.getUnlinkedSVs().size();
        int inconsistentChains = 0;
        int neoChromosomes = 0;

        for (final SvChain chain : cluster.getChains())
        {
            if (!chain.isConsistent())
            {
                ++inconsistentChains;
                continue;
            }

            ++neoChromosomes;
        }

        int armGroupCount = cluster.getArmGroups().size();

        // skip simple chained clusters
        if (cluster.getArmCount() == 1 && cluster.getCount() == 2)
            return;

        // isSpecificCluster(cluster);

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
            if (!chain.isConsistent())
            {
                allChainsConsistent = false;
                continue;
            }

            Map<String, int[]> armDataMap = new HashMap();

            final SvBreakend firstBreakend = chain.getFirstSV().getBreakend(chain.firstLinkOpenOnStart());
            final SvBreakend lastBreakend = chain.getLastSV().getBreakend(chain.lastLinkOpenOnStart());

            final String startChrArm = firstBreakend != null ? firstBreakend.getChrArm() : "";
            final String endChrArm = lastBreakend != null ? lastBreakend.getChrArm() : "";

            if (!startChrArm.isEmpty())
                armDataMap.put(startChrArm, new int[CHAIN_TI_ASMB_COUNT + 1]);

            if (!endChrArm.isEmpty())
                armDataMap.put(endChrArm, new int[CHAIN_TI_ASMB_COUNT + 1]);
        }
    }

            private static double DOUBLE_MINUTE_PLOIDY_THRESHOLD = 8;
    private static double DOUBLE_MINUTE_PLOIDY_GAP_RATIO = 3;

    private void reportDoubleMinutes(final SvCluster cluster)
    {
        // order SVs in descending ploidy order
        List<Double> ploidyList = Lists.newArrayList();
        List<SvVarData> indexSvList = Lists.newArrayList();
        boolean hasHighPloidy = false;

        for(final SvVarData var : cluster.getSVs())
        {
            if(var.isReplicatedSv())
                continue;

            double ploidy = var.getSvData().ploidy();
            int i = 0;
            for(; i < ploidyList.size(); ++i)
            {
                Double otherPloidy = ploidyList.get(i);
                if(ploidy > otherPloidy)
                    break;
            }

            ploidyList.add(i, ploidy);
            indexSvList.add(i, var);

            if(ploidy >= DOUBLE_MINUTE_PLOIDY_THRESHOLD)
                hasHighPloidy = true;
        }

        if(!hasHighPloidy)
            return;

        boolean isPotentialDM = false;

        // keep track of any chromosome with high ploidy as an indication of false DMs and/or under-clustering
        List<String> highPloidyChromosomes = Lists.newArrayList();

        int svsAboveThreshold = 0;
        double minPloidyAboveThreshold = 0;
        for(int i = 0; i < ploidyList.size(); ++i)
        {
            Double ploidy = ploidyList.get(i);

            if(ploidy < DOUBLE_MINUTE_PLOIDY_THRESHOLD)
                break;

            ++svsAboveThreshold;
            minPloidyAboveThreshold = ploidy;

            final SvVarData var = indexSvList.get(i);

            for(int be = SVI_START; be <= SVI_END; ++be)
            {
                boolean useStart = isStart(be);

                if(var.isNullBreakend() && !useStart)
                    continue;

                if(!highPloidyChromosomes.contains(var.chromosome(useStart)))
                    highPloidyChromosomes.add(var.chromosome(useStart));
            }

            // check vs next
            if(i == ploidyList.size() - 1)
            {
                LOGGER.debug(String.format("cluster(%s count=%d) DM highPloidyCount(%d chr=%d) currentSV(%s) ploidy(%.2f) with no others",
                        cluster.id(), cluster.getUniqueSvCount(), svsAboveThreshold, highPloidyChromosomes.size(), var.posId(), ploidy));
                isPotentialDM = true;
                break;
            }
            else
            {
                double nextPloidy = ploidyList.get(i+1);

                if(nextPloidy * DOUBLE_MINUTE_PLOIDY_GAP_RATIO < ploidy)
                {
                    LOGGER.debug(String.format("cluster(%s count=%d) DM highPloidyCount(%d chr=%d) currentSV(%s) ploidy(%.2f) vs next(%.3f)",
                            cluster.id(), cluster.getUniqueSvCount(), svsAboveThreshold, highPloidyChromosomes.size(),
                            indexSvList.get(i).posId(), ploidy, nextPloidy));
                    isPotentialDM = true;
                    break;
                }
            }
        }

        if(isPotentialDM)
        {
            // check for high ploidy in other variants on the relevant chromosomes
            boolean otherClustersHaveHighPloidy = false;

            for(final String chromsome : highPloidyChromosomes)
            {
                final List<SvBreakend> breakendList = mClusteringMethods.getChrBreakendMap().get(chromsome);

                if(breakendList == null)
                    continue;

                for(final SvBreakend breakend : breakendList)
                {
                    if(breakend.getSV().getCluster() == cluster)
                        continue;

                    if(breakend.getSV().getSvData().ploidy() * DOUBLE_MINUTE_PLOIDY_GAP_RATIO >= minPloidyAboveThreshold)
                    {
                        otherClustersHaveHighPloidy = true;
                        break;
                    }
                }

                if(otherClustersHaveHighPloidy)
                    break;
            }

            final String dmAnnotation = otherClustersHaveHighPloidy ? CLUSTER_ANNONTATION_DM + "_Unclear" : CLUSTER_ANNONTATION_DM;
            cluster.addAnnotation(dmAnnotation);
        }

        /*
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
        */
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
        mPcClustering.logStats(false);
        mPcChaining.logStats(false);
        mPcAnnotation.logStats(false);
    }

}
