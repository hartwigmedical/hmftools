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
import static com.hartwig.hmftools.svanalysis.analysis.CNAnalyser.CN_SEG_DATA_MAP_BEFORE;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.DOUBLE_MINUTES;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.FOLDBACK_MATCHES;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.REPLICATION_REPAIR;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.annotateChainedClusters;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.annotateFoldbacks;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.annotateTemplatedInsertions;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.calcNetCopyNumberChangeAcrossCluster;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.checkLooseFoldbacks;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.findIncompleteFoldbackCandidates;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.findPotentialDoubleMinuteClusters;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.reportClusterRepRepairSegments;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.runAnnotation;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.MIN_TEMPLATED_INSERTION_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.haveLinkedAssemblies;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_COMMON_ARMS;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_FOLDBACKS;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_LOH_CHAIN;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_LOOSE_OVERLAP;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_NET_ARM_END_PLOIDY;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.CLUSTER_REASON_BE_PLOIDY_DROP;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.applyCopyNumberReplication;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.findCentromereBreakendIndex;
import static com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator.markLineCluster;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.mergeArmClusters;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_ASSEMBLY_LINK_COUNT;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_LENGTH;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_COMPLEX_CHAIN;
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
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.haveSameChrArms;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLOH;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvaConfig;
import com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ClusterAnalyser {

    private final SvaConfig mConfig;
    private SvClusteringMethods mClusteringMethods;
    private CNAnalyser mCopyNumberAnalyser;
    private SvGeneTranscriptCollection mGeneCollection;

    String mSampleId;
    private List<SvVarData> mAllVariants;
    List<SvCluster> mClusters;
    private ChainFinder mChainFinder;
    private LinkFinder mLinkFinder;

    private boolean mUseAllelePloidies;
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
        mCopyNumberAnalyser = null;
        mGeneCollection = null;
        mClusters = Lists.newArrayList();
        mAllVariants = Lists.newArrayList();
        mSampleId = "";
        mLinkFinder = new LinkFinder();
        mChainFinder = new ChainFinder();
        mChainFinder.setLogVerbose(mConfig.LogVerbose);
        mLinkFinder.setLogVerbose(mConfig.LogVerbose);
        mRunValidationChecks = false;
        mUseAllelePloidies = false;

        mPcClustering = new PerformanceCounter("Clustering");
        mPcChaining = new PerformanceCounter("Chaining");
        mPcAnnotation = new PerformanceCounter("Annotation");
    }

    public void setCopyNumberAnalyser(CNAnalyser cnAnalyser)
    {
        mCopyNumberAnalyser = cnAnalyser;
    }
    public void setGeneCollection(final SvGeneTranscriptCollection geneCollection) { mGeneCollection = geneCollection; }
    public void setUseAllelePloidies(boolean toggle)
    {
        mChainFinder.setUseAllelePloidies(toggle);
        mUseAllelePloidies = toggle;
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
        mPcClustering.pause();

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

        mPcChaining.start();
        findLimitedChains();
        mPcChaining.pause();

        mPcClustering.resume();
        mClusteringMethods.mergeClusters(mSampleId, mClusters);
        mPcClustering.pause();

        // log basic clustering details
        for(SvCluster cluster : mClusters)
        {
            if(cluster.getSvCount() > 1)
                cluster.logDetails();
        }

        // INVs and other SV-pairs which make foldbacks are now used in the inconsistent clustering logic
        markFoldbacks();

        mPcClustering.resume();
        mergeClusters();
        mPcClustering.stop();

        mPcChaining.resume();
        findLinksAndChains();
        mPcChaining.stop();

        if(mRunValidationChecks)
        {
            if(!mClusteringMethods.validateClustering(mClusters) || !mClusteringMethods.checkClusterDuplicates(mClusters))
            {
                LOGGER.info("exiting with cluster-validation errors");
                return;
            }
        }

        mPcAnnotation.start();

        // final clean-up and analysis
        for(SvCluster cluster : mClusters)
        {
            if(!cluster.isResolved() && cluster.getResolvedType() != RESOLVED_TYPE_NONE)
            {
                // any cluster with a long DEL or DUP not merged can now be marked as resolved
                if(cluster.getSvCount() == 1 &&  cluster.getResolvedType() == RESOLVED_TYPE_SIMPLE_SV)
                    cluster.setResolved(true, RESOLVED_TYPE_SIMPLE_SV);

                // mark off any fully chained simple clusters
                if(cluster.getSvCount() == 1 &&  cluster.getResolvedType() == RESOLVED_TYPE_SIMPLE_CHAIN)
                    cluster.setResolved(true, RESOLVED_TYPE_SIMPLE_CHAIN);
            }

            mergeArmClusters(cluster.getArmClusters());

            // if(LOGGER.isDebugEnabled())
            //     logArmClusterData(cluster);

            reportClusterFeatures(cluster);
        }

        reportOtherFeatures();

        mPcAnnotation.stop();

        // validation-only: checkClusterDuplicates(mClusters);
    }

    public void findLimitedChains()
    {
        // chain small clusters and only assembled links in larger ones
        for(SvCluster cluster : mClusters)
        {
            if(cluster.isSimpleSingleSV())
            {
                setClusterResolvedState(cluster);
                continue;
            }

            // more complicated clusters for now
            boolean isSimple = cluster.getSvCount() <= SMALL_CLUSTER_SIZE && cluster.isConsistent() && !cluster.hasVariedCopyNumber();

            // inferred links are used to classify simple resolved types involving 2-3 SVs
            mLinkFinder.findLinkedPairs(cluster, isSimple);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster, !isSimple);

            cluster.cacheLinkedPairs();

            if(isSimple)
            {
                setClusterResolvedState(cluster);

                if(cluster.isFullyChained() && cluster.isConsistent())
                {
                    LOGGER.debug("cluster({}) simple and consistent with {} SVs", cluster.id(), cluster.getSvCount());
                }
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
            if (cluster.isSimpleSingleSV() || cluster.isFullyChained() || cluster.getSvCount() < 2)
                continue;

            cluster.dissolveLinksAndChains();

            cluster.removeReplicatedSvs();
            applyCopyNumberReplication(cluster);

            // first establish links between SVs (eg TIs and DBs)
            mLinkFinder.findLinkedPairs(cluster, false);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster, false);

            cluster.cacheLinkedPairs();
            setClusterResolvedState(cluster);
            cluster.logDetails();
        }
    }

    public void mergeClusters()
    {
        boolean foundMerges = mergeInconsistentClusters();

        int iterations = 1;

        while(foundMerges)
        {
            foundMerges = mergeInconsistentClusters();
            ++iterations;

            if(iterations > 20)
            {
                LOGGER.warn("reached {} iterations of clustering merging", iterations);
                break;
            }
        }
    }

    private void findChains(SvCluster cluster, boolean assembledLinksOnly)
    {
        mChainFinder.initialise(cluster);

        // isSpecificCluster(cluster);

        if(mConfig.MaxClusterSize > 0 && (cluster.getSvCount() > mConfig.MaxClusterSize || cluster.getSvCount(true) > mConfig.MaxClusterSize * 5))
        {
            LOGGER.info("cluster({}) skipping large cluster: unique({}) replicated({})",
                    cluster.id(), cluster.getSvCount(), cluster.getSvCount(true));
            return;
        }

        mChainFinder.formClusterChains(assembledLinksOnly);
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
        if (cluster.getSvCount() == 2 && cluster.isConsistent())
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
                SvBreakend breakend1 = cluster.getSV(0).getBreakend(true);
                SvBreakend breakend2 = cluster.getSV(1).getBreakend(true);

                if(breakend1.orientation() != breakend2.orientation())
                {
                    long length = abs(breakend1.position() - breakend2.position());
                    int minTiLength = getMinTemplatedInsertionLength(breakend1, breakend2);

                    if(length < minTiLength)
                    {
                        cluster.setResolved(true, RESOLVED_TYPE_SGL_PAIR_INS);
                    }
                    else
                    {
                        boolean v1First = breakend1.position() < breakend2.position();
                        boolean v1PosOrientation = (breakend1.orientation() == 1);

                        if(v1First == v1PosOrientation)
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

    private boolean mergeInconsistentClusters()
    {
        // second round of cluster merging on more complex criteria and inconsistencies:
        // merge on foldbacks on the same arm
        // merge on links between common arms
        // merge if one cluster has footprints which overlap unresolved complex SVs
        // merge clusters which resolve another's LOH DUP
        boolean foundClustersToMerge = false;

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
                    ++index2;
                    continue;
                }

                foundClustersToMerge = true;
                cluster1.mergeOtherCluster(cluster2, false);

                mClusters.remove(index2);
            }

            ++index1;
        }

        if(mergeTraversingClusters())
            foundClustersToMerge = true;

        if(mergeLOHResolvingClusters())
            foundClustersToMerge = true;

        if(mergePloidyResolvingClusters())
            foundClustersToMerge = true;

        if(mergeNetCopyNumberResolvingClusters())
            foundClustersToMerge = true;

        return foundClustersToMerge;
    }

    private boolean mergePloidyResolvingClusters()
    {
        if(!mUseAllelePloidies)
            return false;

        List<SvCluster> mergedClusters = Lists.newArrayList();

        int clusterIndex = 0;
        while(clusterIndex < mClusters.size())
        {
            SvCluster cluster = mClusters.get(clusterIndex);

            if(mergedClusters.contains(cluster) || cluster.isResolved())
            {
                ++clusterIndex;
                continue;
            }

            boolean mergedOtherClusters = false;

            for (final Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
            {
                List<SvBreakend> breakendList = entry.getValue();
                long lowerArmBoundary = breakendList.get(0).position();
                long upperArmBoundary = breakendList.get(breakendList.size() - 1).position();

                List<SvBreakend> fullBreakendList = mClusteringMethods.getChrBreakendMap().get(entry.getKey());

                for (SvBreakend breakend : breakendList)
                {
                    double breakendPloidy = breakend.ploidy();
                    SvCluster resolvingCluster = null;
                    SvBreakend resolvingBreakend = null;

                    boolean traverseUp = breakend.orientation() == -1;
                    int index = breakend.getChrPosIndex();

                    while(true)
                    {
                        index += traverseUp ? 1 : -1;

                        if(index < 0 || index >= fullBreakendList.size())
                            break;

                        SvBreakend nextBreakend = fullBreakendList.get(index);

                        if(nextBreakend.arm() != breakend.arm())
                            break;

                        // only cluster with variants within this cluster's boundaries
                        if(nextBreakend.position() < lowerArmBoundary || nextBreakend.position() > upperArmBoundary)
                            break;

                        if(nextBreakend.orientation() == breakend.orientation())
                            continue;

                        SvCluster nextCluster = nextBreakend.getSV().getCluster();
                        if(nextCluster == cluster)
                            break;

                        if(nextCluster.isResolved() || mergedClusters.contains(nextCluster))
                            continue;

                        resolvingCluster = nextCluster;
                        resolvingBreakend = nextBreakend;

                        double majorAP = nextBreakend.majorAllelePloidy(!traverseUp);

                        if(majorAP < breakendPloidy - 1 && !copyNumbersEqual(majorAP, breakendPloidy))
                        {
                            LOGGER.debug("cluster({}) SV({}) requires cluster({}) breakend({}) prior to MAP drop({})",
                                    cluster.id(), breakend.getSV().posId(), resolvingCluster.id(), resolvingBreakend.toString(),
                                    String.format("%.2f -> %.2f", breakendPloidy, majorAP));

                            breakend.getSV().addClusterReason(CLUSTER_REASON_BE_PLOIDY_DROP, nextBreakend.getSV().id());
                            nextBreakend.getSV().addClusterReason(CLUSTER_REASON_BE_PLOIDY_DROP, breakend.getSV().id());

                            resolvingCluster.addClusterReason(CLUSTER_REASON_BE_PLOIDY_DROP);
                            cluster.addClusterReason(CLUSTER_REASON_BE_PLOIDY_DROP);

                            cluster.mergeOtherCluster(resolvingCluster);

                            /*
                            // Time,SampleId,ClusterId,SvId1,SvId2,BreakendPloidy,MajorAP
                            LOGGER.info(String.format("CR_BE_PLOIDY: %s,%d,%s,%s,%.2f,%.2f",
                                    mSampleId, cluster.id(), breakend.getSV().id(), resolvingBreakend.getSV().id(),
                                    breakendPloidy, majorAP));
                            */

                            mergedClusters.add(resolvingCluster);

                            mergedOtherClusters = true;
                            break;
                        }
                    }

                    if(mergedOtherClusters)
                        break;
                }

                if(mergedOtherClusters)
                    break;
            }

            if(mergedOtherClusters)
            {
                // repeat this cluster
            }
            else
            {
                ++clusterIndex;
            }
        }

        if(mergedClusters.isEmpty())
            return false;

        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }

    private boolean mergeLOHResolvingClusters()
    {
        // merge clusters if one resolves another's LOH event with a DUP on one side
        List<SvCluster> clustersWithLohEvents = mClusters.stream()
                .filter(x -> !x.getLohEvents().isEmpty())
                .filter(x -> !x.hasLinkingLineElements())
                .collect(Collectors.toList());

        List<SvCluster> mergedClusters = Lists.newArrayList();

        for(SvCluster lohCluster : clustersWithLohEvents)
        {
            if(mergedClusters.contains(lohCluster)) // if this has LOH events they will have been added to the parent cluster
                continue;

            List<SvLOH> lohEvents = lohCluster.getLohEvents();

            int lohIndex = 0; // used since merging another cluster can add more LOH events
            while(lohIndex < lohEvents.size())
            {
                SvLOH lohEvent = lohEvents.get(lohIndex);

                if(!lohEvent.IsValid)
                {
                    ++lohIndex;
                    continue;
                }

                for(int be = SVI_START; be <= SVI_END; ++be)
                {
                    SvBreakend lohBreakend = lohEvent.getBreakend(isStart(be));

                    if(lohBreakend == null || lohBreakend.getSV().type() != DUP)
                        continue;

                    // it's possible that the breakends for this LOH are not clustered, eg if one is LINE
                    if(lohBreakend.getSV().getCluster() != lohCluster)
                        continue;

                    // walk towards the LOH from the other end of this DUP to see if it can find a resolving event within the cluster
                    List<SvBreakend> fullBreakendList = mClusteringMethods.getChrBreakendMap().get(lohBreakend.chromosome());

                    SvBreakend otherBreakend = lohBreakend.getOtherBreakend();
                    int index = otherBreakend.getChrPosIndex();
                    boolean traverseUp = otherBreakend.orientation() == -1;
                    SvCluster resolvingCluster = null;
                    SvBreakend resolvingBreakend = null;

                    while(true)
                    {
                        index += traverseUp ? 1 : -1;

                        if(index < 0 || index >= fullBreakendList.size())
                            break;

                        SvBreakend nextBreakend = fullBreakendList.get(index);

                        if(nextBreakend == lohBreakend)
                        {
                            // the LOH was reached without finding an offsetting SV
                            if(resolvingCluster == null)
                                break;

                            LOGGER.debug("cluster({}) SV({}) resolved prior to LOH by other cluster({}) breakend({})",
                                    lohCluster.id(), lohBreakend.getSV().posId(), resolvingCluster.id(), resolvingBreakend.toString());

                            lohBreakend.getSV().addClusterReason(CLUSTER_REASON_LOH_CHAIN, resolvingBreakend.getSV().id());
                            resolvingBreakend.getSV().addClusterReason(CLUSTER_REASON_LOH_CHAIN, lohBreakend.getSV().id());

                            resolvingCluster.addClusterReason(CLUSTER_REASON_LOH_CHAIN);
                            lohCluster.addClusterReason(CLUSTER_REASON_LOH_CHAIN);

                            lohCluster.mergeOtherCluster(resolvingCluster);

                            mergedClusters.add(resolvingCluster);
                            break;
                        }

                        if(nextBreakend.orientation() == otherBreakend.orientation())
                            continue;

                        SvCluster otherCluster = nextBreakend.getSV().getCluster();

                        if(otherCluster == lohCluster)
                            break; // own cluster resolves this LOH breakend

                        if(mergedClusters.contains(otherCluster) || otherCluster.isResolved())
                            continue;

                        if(resolvingCluster != null)
                            break; // cannot apply this rule if more than 1 cluster meet the conditions

                        // found an option, but continue on to see if any other clusters also satisfy the same conditions
                        resolvingBreakend = nextBreakend;
                        resolvingCluster = otherCluster;
                    }
                }

                ++lohIndex;
            }
        }

        if(mergedClusters.isEmpty())
            return false;

        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }

    private boolean mergeNetCopyNumberResolvingClusters()
    {
        if(!mUseAllelePloidies)
            return false;

        // calculate the net CN change across all breakends on an arm
        // if it needs to be resolved prior to the telomere or centromere and another cluster
        // can help do that, then merge in that cluster

        /* for the first and last uninterrupted footprint of each cluster on each chromosome, calculate the minimal number
        of telomeric/centromeric facing ploidy that is required to explain the orientation and ploidy of the breakends limited
        to the major allele ploidy immediately flanking the cluster

          If this exceeds the telomeric / centromeric major allele ploidy then search for other (non-resolved) clusters that
          could explain the drop in ploidy. If there is only one cluster that can explain the full change in major allele
          ploidy cluster with that, else choose the nearest cluster.
        */

        List<SvCluster> mergedClusters = Lists.newArrayList();

        int clusterIndex = 0;
        while(clusterIndex < mClusters.size())
        {
            SvCluster cluster = mClusters.get(clusterIndex);

            if(mergedClusters.contains(cluster) || cluster.isResolved())
            {
                ++clusterIndex;
                continue;
            }

            boolean mergedOtherClusters = false;

            for (final Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
            {
                final String chromosome = entry.getKey();
                List<SvBreakend> breakendList = entry.getValue();
                List<SvCNData> cnDataList = mCopyNumberAnalyser.getChrCnDataMap().get(chromosome);

                if(cnDataList == null || cnDataList.isEmpty())
                    continue;

                List<SvBreakend> fullBreakendList = mClusteringMethods.getChrBreakendMap().get(entry.getKey());

                for(int armIndex = 0; armIndex <= 1; ++armIndex)
                {
                    String arm = (armIndex == 0) ? CHROMOSOME_ARM_P : CHROMOSOME_ARM_Q;

                    int centomereBreakendIndex = findCentromereBreakendIndex(breakendList, arm);

                    if(centomereBreakendIndex == -1)
                        continue; // no breakends on this arm

                    SvBreakend lowerBreakend = arm == CHROMOSOME_ARM_P ? breakendList.get(0) : breakendList.get(centomereBreakendIndex);

                    SvBreakend upperBreakend = arm == CHROMOSOME_ARM_P ? breakendList.get(centomereBreakendIndex)
                            : breakendList.get(breakendList.size() - 1);

                    double[] boundaryCNData = calcNetCopyNumberChangeAcrossCluster(breakendList, arm, true);
                    double telomereMinFacingPloidy = boundaryCNData[0];
                    double centromereMinFacingPloidy = boundaryCNData[1];

                    double[] centromereCNData = mCopyNumberAnalyser.getCentromereCopyNumberData(chromosome, arm.equals(CHROMOSOME_ARM_P));

                    SvCNData telemoreData = arm == CHROMOSOME_ARM_P ? cnDataList.get(0) : cnDataList.get(cnDataList.size() - 1);
                    double telomereMAP = telemoreData.majorAllelePloidy();

                    // now look towards the telomere and centromere and work out there is likely a breakend missing
                    // which would explain the net CN change
                    for(int directionIndex = 0; directionIndex <= 1; ++directionIndex)
                    {
                        boolean traverseUp = (directionIndex == 0);
                        SvBreakend clusterBreakend = traverseUp ? upperBreakend : lowerBreakend;
                        boolean facingCentromere = traverseUp == (arm == CHROMOSOME_ARM_P);

                        double tOrCMinFacingPloidy = facingCentromere ? centromereMinFacingPloidy : telomereMinFacingPloidy;
                        double clusterBoundaryMAP = clusterBreakend.majorAllelePloidy(!traverseUp);
                        double clusterBoundaryMinPloidy = min(tOrCMinFacingPloidy, clusterBoundaryMAP);

                        double armEndMAP = facingCentromere ? centromereCNData[CN_SEG_DATA_MAP_BEFORE] : telomereMAP;

                        if(clusterBoundaryMinPloidy < armEndMAP || copyNumbersEqual(armEndMAP, clusterBoundaryMinPloidy))
                            continue;

                        int index = clusterBreakend.getChrPosIndex();
                        while (true)
                        {
                            index += traverseUp ? 1 : -1;

                            if (index < 0 || index >= fullBreakendList.size())
                                break;

                            SvBreakend nextBreakend = fullBreakendList.get(index);

                            if (nextBreakend.arm() != arm)
                                break;

                            SvCluster otherCluster = nextBreakend.getSV().getCluster();

                            if (otherCluster.isResolved() || otherCluster == cluster || mergedClusters.contains(otherCluster))
                                continue;

                            if (nextBreakend.orientation() != clusterBreakend.orientation())
                            {
                                // candidate found
                                LOGGER.debug("cluster({}) arm({} facing {}) outerBE({}) requires cluster({}) breakend({}) prior to CN drop({})",
                                        cluster.id(), arm, facingCentromere ? "centro" : "telo",
                                        clusterBreakend.toString(), otherCluster.id(), nextBreakend.toString(),
                                        String.format("cluster minPloidy=%.2f map=%.2f arm map=%.2f",
                                                clusterBoundaryMinPloidy, clusterBoundaryMAP, armEndMAP));

                                clusterBreakend.getSV().addClusterReason(CLUSTER_REASON_NET_ARM_END_PLOIDY, nextBreakend.getSV().id());
                                nextBreakend.getSV().addClusterReason(CLUSTER_REASON_NET_ARM_END_PLOIDY, clusterBreakend.getSV().id());

                                otherCluster.addClusterReason(CLUSTER_REASON_NET_ARM_END_PLOIDY);
                                cluster.addClusterReason(CLUSTER_REASON_NET_ARM_END_PLOIDY);

                                /*
                                // TEMP logging for further analysis
                                // Time,SampleId,ClusterId,SvId1,SvId2,ClusterBoundaryMinPloidy,ClusterBoundaryMAP,ArmEndMAP,FacingCentromere
                                LOGGER.info(String.format("CR_ARM_END_PLOIDY: %s,%d,%s,%s,%.2f,%.2f,%.2f,%s",
                                        mSampleId, cluster.id(), clusterBreakend.getSV().id(), nextBreakend.getSV().id(),
                                        clusterBoundaryMinPloidy, clusterBoundaryMAP, armEndMAP, facingCentromere));
                                */

                                cluster.mergeOtherCluster(otherCluster);

                                mergedClusters.add(otherCluster);

                                mergedOtherClusters = true;
                                break;
                            }
                        }

                        if(mergedOtherClusters)
                            break;

                    } // end each direction within an arm

                    if(mergedOtherClusters)
                        break;

                } // end each arm

                if(mergedOtherClusters)
                    break;

            } // end each chromosome

            if(mergedOtherClusters)
            {
                // repeat this cluster
            }
            else
            {
                ++clusterIndex;
            }
        }

        if(mergedClusters.isEmpty())
            return false;

        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }

    private static int MAX_FOLDBACK_NEXT_CLUSTER_DISTANCE = 5000000;

    private boolean canMergeClustersOnFoldbacks(final SvCluster cluster1, final SvCluster cluster2)
    {
        // merge any clusters with foldbacks on the same arm
        final List<SvVarData> cluster1Foldbacks = cluster1.getFoldbacks();
        final List<SvVarData> cluster2Foldbacks = cluster2.getFoldbacks();

        if (cluster1Foldbacks.isEmpty() && cluster2Foldbacks.isEmpty())
            return false;

        if (!cluster1Foldbacks.isEmpty() && !cluster2Foldbacks.isEmpty())
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

                            var1.addClusterReason(CLUSTER_REASON_FOLDBACKS, var2.id());
                            var2.addClusterReason(CLUSTER_REASON_FOLDBACKS, var1.id());

                            cluster1.addClusterReason(CLUSTER_REASON_FOLDBACKS);
                            cluster2.addClusterReason(CLUSTER_REASON_FOLDBACKS);
                            return true;
                        }
                    }
                }
            }
        }

        final Map<String, List<SvBreakend>> chrBreakendMap = mClusteringMethods.getChrBreakendMap();

        // additionally check whether any of the foldbacks face an opposing SV in the other cluster,
        // and they must be the next SV and within 5M bases away
        for (int i = 0; i <= 1; ++i)
        {
            final List<SvVarData> foldbacks = (i == 0) ? cluster1Foldbacks : cluster2Foldbacks;
            final SvCluster foldbackCluster = (i == 0) ? cluster1 : cluster2;
            final SvCluster otherCluster = (i == 0) ? cluster2 : cluster1;

            for (final SvVarData var : foldbacks)
            {
                // get the inner-most breakend
                SvBreakend foldbackBreakend = null;

                if (!var.isChainedFoldback())
                {
                    foldbackBreakend = var.orientation(true) == 1 ? var.getBreakend(true) : var.getBreakend(false);
                }
                else
                {
                    SvBreakend be1 = var.getFoldbackBreakend(true) != null ? var.getBreakend(true) : var.getBreakend(false);
                    SvBreakend be2 = var.getChainedFoldbackBreakend();
                    foldbackBreakend = (var.orientation(true) == 1) == (be1.position() < be2.position()) ? be1 : be2;
                }

                final List<SvBreakend> breakendList = chrBreakendMap.get(foldbackBreakend.chromosome());

                SvBreakend nextBreakend = getNextUnresolvedBreakend(foldbackBreakend, breakendList);

                if (nextBreakend == null || nextBreakend.orientation() == foldbackBreakend.orientation()
                || nextBreakend.getSV().getCluster() != otherCluster)
                    continue;

                if (abs(nextBreakend.position() - foldbackBreakend.position()) > MAX_FOLDBACK_NEXT_CLUSTER_DISTANCE)
                    continue;

                double fbPloidy = foldbackBreakend.ploidy();
                double nbPloidy = nextBreakend.ploidy();

                if (nbPloidy < fbPloidy && !copyNumbersEqual(nbPloidy, fbPloidy))
                    continue;

                LOGGER.debug("cluster({}) foldback breakend({}) faces cluster({}) breakend({})",
                        foldbackCluster.id(), foldbackBreakend.toString(), otherCluster.id(), nextBreakend.toString());

                var.addClusterReason(CLUSTER_REASON_FOLDBACKS, nextBreakend.getSV().id());
                nextBreakend.getSV().addClusterReason(CLUSTER_REASON_FOLDBACKS, var.id());

                foldbackCluster.addClusterReason(CLUSTER_REASON_FOLDBACKS);
                otherCluster.addClusterReason(CLUSTER_REASON_FOLDBACKS);

                return true;
            }
        }

        return false;
    }

    private final SvBreakend getNextUnresolvedBreakend(final SvBreakend foldbackBreakend, final List<SvBreakend> breakendList)
    {
        // select the next breakend after this foldback if it's in a different, unresolved cluster
        boolean traverseUp = foldbackBreakend.orientation() == -1;
        int startIndex = traverseUp ? foldbackBreakend.getChrPosIndex() + 1 : foldbackBreakend.getChrPosIndex() - 1;
        final SvCluster fbCluster = foldbackBreakend.getSV().getCluster();

        int index = startIndex;

        while(index >= 0 && index < breakendList.size())
        {
            final SvBreakend breakend = breakendList.get(index);
            final SvCluster cluster = breakend.getSV().getCluster();

            if(!cluster.isResolved())
            {
                if(cluster != fbCluster)
                    return breakend;
                else
                    return null;
            }

            if(traverseUp)
                ++index;
            else
                --index;
        }

        return null;
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

                        // final String commonArms = var1.id() + "_" + var2.id();

                        var1.addClusterReason(CLUSTER_REASON_COMMON_ARMS, var2.id());
                        var2.addClusterReason(CLUSTER_REASON_COMMON_ARMS, var1.id());

                        cluster1.addClusterReason(CLUSTER_REASON_COMMON_ARMS);
                        cluster2.addClusterReason(CLUSTER_REASON_COMMON_ARMS);
                        return true;
                    }
                }
            }
        }

        return false;
    }

    private boolean mergeTraversingClusters()
    {
        // merge clusters if one is straddled by opposing breakends in the other
        List<SvCluster> mergedClusters = Lists.newArrayList();

        int clusterIndex = 0;
        while(clusterIndex < mClusters.size())
        {
            SvCluster cluster = mClusters.get(clusterIndex);

            if(mergedClusters.contains(cluster) || cluster.isResolved())
            {
                ++clusterIndex;
                continue;
            }

            boolean mergedOtherClusters = false;

            for (final Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
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
                        final SvBreakend otherBreakend = fullBreakendList.get(j);

                        if(otherBreakend.getSV().type() != BND && otherBreakend.getSV().type() != SGL)
                            continue;

                        if(otherBreakend.orientation() == lowerBreakend.orientation())
                            continue;

                        final SvCluster otherCluster = otherBreakend.getSV().getCluster();

                        if (otherCluster == cluster || otherCluster.isResolved() || mergedClusters.contains(otherCluster))
                            continue;

                        LOGGER.debug("cluster({}) breakends({} & {}) overlap cluster({}) breakend({})",
                                cluster.id(), lowerBreakend.toString(), upperBreakend.toString(), otherCluster.id(), otherBreakend.toString());

                        otherBreakend.getSV().addClusterReason(CLUSTER_REASON_LOOSE_OVERLAP, lowerBreakend.getSV().id());
                        lowerBreakend.getSV().addClusterReason(CLUSTER_REASON_LOOSE_OVERLAP, otherBreakend.getSV().id());

                        otherCluster.addClusterReason(CLUSTER_REASON_LOOSE_OVERLAP);
                        cluster.addClusterReason(CLUSTER_REASON_LOOSE_OVERLAP);

                        cluster.mergeOtherCluster(otherCluster);

                        mergedClusters.add(otherCluster);

                        mergedOtherClusters = true;
                        break;
                    }

                    if(mergedOtherClusters)
                        break;
                }

                if(mergedOtherClusters)
                    break;
            }

            if(mergedOtherClusters)
            {
                // repeat this cluster
            }
            else
            {
                ++clusterIndex;
            }
        }

        if(mergedClusters.isEmpty())
            return false;

        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }

    public void markFoldbacks()
    {
        // find all valid consective breakends formed either from a single SV or a chained set
        for(final Map.Entry<String, List<SvBreakend>> entry : mClusteringMethods.getChrBreakendMap().entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            for(int i = 0; i < breakendList.size() - 1; ++i)
            {
                SvBreakend breakend = breakendList.get(i);

                if(breakend.isAssembledLink())
                    continue;

                SvBreakend nextBreakend = breakendList.get(i + 1);

                SvBreakend beFront = null; // the lower position for orientation +1 and vice versa
                SvBreakend beBack = null;

                int skippedAssemblies = 0;

                int j = i + 1;
                while(j < breakendList.size())
                {
                    nextBreakend = breakendList.get(j);

                    // first skip over an breakends in a DB with the initial breakend
                    if(j == i + 1 && breakend.orientation() == -1 && nextBreakend.orientation() == 1
                    && nextBreakend.position() - breakend.position() < getMinTemplatedInsertionLength(nextBreakend, breakend))
                    {
                        ++j;
                        continue;
                    }

                    // check check for any assembled links in between the potential foldback breakends
                    if(j + 1 < breakendList.size() && nextBreakend.isAssembledLink()
                    && nextBreakend.getSV().getLinkedPair(nextBreakend.usesStart()) != null)
                    {
                        SvLinkedPair asmbLink = nextBreakend.getSV().getLinkedPair(nextBreakend.usesStart());
                        SvBreakend nextNextBreakend = breakendList.get(j + 1);
                        if(asmbLink.getOtherBreakend(nextBreakend) == nextNextBreakend)
                        {
                            // skip over both these assembled links
                            j += 2;
                            ++skippedAssemblies;
                            continue;
                        }
                    }

                    // check again for an overlapping DB at the outer (potential) foldback breakend
                    if(breakend.orientation() == 1 && nextBreakend.orientation() == -1 && j < breakendList.size() - 1)
                    {
                        SvBreakend nextNextBreakend = breakendList.get(j + 1);

                        if(nextNextBreakend.orientation() == breakend.orientation()
                        && nextNextBreakend.position() - nextBreakend.position() < getMinTemplatedInsertionLength(nextBreakend, nextNextBreakend))
                        {
                            nextBreakend = nextNextBreakend;
                        }
                    }

                    // now check for opposite orientation for the potential foldback
                    if(nextBreakend.orientation() == breakend.orientation() && !nextBreakend.isAssembledLink())
                    {
                        beFront = breakend.orientation() == 1 ? breakend : nextBreakend;
                        beBack = breakend.orientation() == 1 ? nextBreakend : breakend;
                    }

                    break;
                }

                boolean foldbackFound = false;

                if(beFront != null && beBack != null)
                {
                    // the foldback is invalid if it has a deletion bridge with overhang on the front-facing breakend
                    final SvLinkedPair dbLink = beFront.getSV().getDBLink(beFront.usesStart());

                    if(dbLink == null || dbLink.length() > 0)
                    {
                        foldbackFound = checkFoldbackBreakends(beFront, beBack);
                    }
                }

                if(!foldbackFound)
                {
                    checkReplicatedBreakendFoldback(breakend);
                }
            }
        }
    }

    public static int MAX_FOLDBACK_CHAIN_LENGTH = 5000;

    private boolean checkFoldbackBreakends(SvBreakend beStart, SvBreakend beEnd)
    {
        // SVs are marked as being in a foldback if they are consecutive breakends,
        // have the same orientation, and are either an INV or part of a chain

        // beStart is the one with the lower position

        SvVarData varEnd = beEnd.getSV();
        SvVarData varStart = beStart.getSV();

        if(varEnd.type() == INS || varStart.type() == INS)
            return false;

        // isSpecificSV(varStart);

        // skip unclustered DELs & DUPs, reciprocal INV or reciprocal BNDs
        final SvCluster cluster1 = varEnd.getCluster();

        if(cluster1.isResolved() && cluster1.getTypeCount(INV) == 0)
            return false;

        boolean singleSV = varEnd.equals(varStart);

        SvCluster cluster2 = null;

        if(singleSV)
        {
            cluster2 = cluster1;
        }
        else
        {
            cluster2 = varStart.getCluster();

            // must be same cluster
            if(cluster1 != cluster2)
                return false;
        }

        boolean beEndUsesStart = beEnd.usesStart();
        boolean beStartUsesStart = beStart.usesStart();

        String chainInfo = "0;0;0";

        if(singleSV)
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
                        return false;
                }
            }
        }
        else
        {
            // must be same cluster and part of the same chain
            if(cluster1 != cluster2)
                return false;

            if(varEnd.getReplicatedCount() != varStart.getReplicatedCount())
                return false;

            final SvChain chain1 = cluster1.findChain(varEnd);
            final SvChain chain2 = cluster2.findChain(varStart);

            if(chain1 == null || chain2 == null || chain1 != chain2)
                return false;

            // check if a path can be walked between these 2 breakends along the chain
            // without going back through this foldback point
            int[] chainData = chain1.breakendsAreChained(varEnd, !beEndUsesStart, varStart, !beStartUsesStart);

            if(chainData[CHAIN_LINK_COUNT] == 0 ) // || chainData[CHAIN_LINK_COUNT] != chainData[CHAIN_ASSEMBLY_LINK_COUNT]
                return false;

            int chainLength = chainData[CHAIN_LENGTH];

            if(chainLength > MAX_FOLDBACK_CHAIN_LENGTH)
            {
                /*
                LOGGER.info("sample({}) chained foldback breakends({} and {}) have long length({} links={} asmb={})",
                        mSampleId, beEnd.toString(), beStart.toString(),
                        chainLength, chainData[CHAIN_LINK_COUNT], chainData[CHAIN_ASSEMBLY_LINK_COUNT]);
                */
                return false;
            }

            chainInfo = String.format("%d;%d;%d",
                    chainData[CHAIN_LINK_COUNT], chainData[CHAIN_ASSEMBLY_LINK_COUNT], chainLength);
        }

        int length = (int)abs(beEnd.position() - beStart.position());

        // if either variant already has foldback info set, favour
        // a) simple inversions then
        // b) shortest length

        boolean skipFoldback = false;
        if(varEnd.getFoldbackBreakend(beEndUsesStart) != null && !singleSV && varEnd.getFoldbackLength(beEndUsesStart) < length)
        {
            skipFoldback = true;
        }
        else if(varStart.getFoldbackBreakend(beStartUsesStart) != null && !singleSV && varStart.getFoldbackLength(beStartUsesStart) < length)
        {
            skipFoldback = true;
        }

        if(skipFoldback)
            return false;

        if(varEnd.getFoldbackBreakend(beEndUsesStart) != null)
        {
            final SvBreakend otherBreakend = varEnd.getFoldbackBreakend(beEndUsesStart);
            otherBreakend.getSV().setFoldbackLink(otherBreakend.usesStart(), null, -1, "");
        }

        if(varStart.getFoldbackBreakend(beStartUsesStart) != null)
        {
            final SvBreakend otherBreakend = varStart.getFoldbackBreakend(beStartUsesStart);
            otherBreakend.getSV().setFoldbackLink(otherBreakend.usesStart(), null, -1, "");
        }

        varEnd.setFoldbackLink(beEndUsesStart, beStart, length, chainInfo);
        varStart.setFoldbackLink(beStartUsesStart, beEnd, length, chainInfo);

        if(varEnd.equals(varStart))
        {
            LOGGER.debug("cluster({}) foldback inversion SV({}) length({})",
                    cluster1.id(), varEnd.posId(), length);
        }
        else
        {
            LOGGER.debug("cluster({}) foldback be1({}) be2({}) length({})",
                    cluster1.id(), beEnd.toString(), beStart.toString(), length);
        }

        return true;
    }

    private void checkReplicatedBreakendFoldback(SvBreakend be)
    {
        // a special case where one ends of an SV connects to both ends of a single other variant
        // during a replication event and in doing so forms a foldback
        final SvVarData var = be.getSV();

        // if(var.getReplicatedCount() < 2)
        //     return;

        if(!var.getAssemblyMatchType(true).equals(ASSEMBLY_MATCH_MATCHED)
        && !var.getAssemblyMatchType(false).equals(ASSEMBLY_MATCH_MATCHED))
        {
            return;
        }

        // replication for chaining isn't done at this earlier stage any more so a simple check is made
        SvLinkedPair remoteTiLink = var.getLinkedPair(!be.usesStart());

        if(remoteTiLink != null)
        {
            SvBreakend remoteBreakend = var.getBreakend(!be.usesStart());
            SvBreakend otherBreakend = remoteTiLink.getOtherBreakend(remoteBreakend);
            SvVarData otherVar = otherBreakend.getSV();

            if (haveLinkedAssemblies(var, otherVar, remoteBreakend.usesStart(), !otherBreakend.usesStart()))
            {
                final String chainInfo = String.format("%d;%d;%d", 1, 1, remoteTiLink.length());

                var.setFoldbackLink(be.usesStart(), be, 0, chainInfo);

                LOGGER.debug("cluster({}) foldback translocation SV({} : {}) with own breakend({})",
                        var.getCluster().id(), var.posId(), var.type(), be.toString());
            }
        }

        /*
        // check if the replicated SV has the same linked pairing or if the variant forms both ends of the cluster's chain
        final SvCluster cluster = var.getCluster();

        for(final SvChain chain : cluster.getChains())
        {
            if(chain.getFirstSV().equals(var, true)
            && chain.getLastSV().equals(var, true)
            && chain.firstLinkOpenOnStart() == chain.lastLinkOpenOnStart())
            {
                final String chainInfo = String.format("%d;%d;%d",
                        chain.getLinkCount(), chain.getAssemblyLinkCount(), chain.getLength(false));

                boolean foldbackIsStart = chain.firstLinkOpenOnStart();

                var.setFoldbackLink(foldbackIsStart, var.getBreakend(foldbackIsStart), 0, chainInfo);

                LOGGER.info("cluster({}) foldback translocation SV({} : {}) with self on {}",
                        cluster.id(), var.posId(), var.type(), foldbackIsStart ? "start" : "end");
            }
        }
        */
    }

    private void reportOtherFeatures()
    {
        annotateTemplatedInsertions(mClusters, mClusteringMethods.getChrBreakendMap());
        // checkSkippedLOHEvents();

        // annotateFoldbacks(mClusters); // unused for now

        if(runAnnotation(mConfig.RequiredAnnotations, DOUBLE_MINUTES))
        {
            findPotentialDoubleMinuteClusters(mSampleId, mClusteringMethods.getChrBreakendMap(), mCopyNumberAnalyser, mGeneCollection);
        }
    }

    private void reportClusterFeatures(final SvCluster cluster)
    {
        // reportDoubleMinutes(cluster, mClusteringMethods.getChrBreakendMap());

        checkLooseFoldbacks(cluster); // required for cluster-level arm group info (eg COMPLEX OTHER)

        annotateChainedClusters(cluster, mClusteringMethods.getProximityDistance());

        if(runAnnotation(mConfig.RequiredAnnotations, FOLDBACK_MATCHES))
        {
            findIncompleteFoldbackCandidates(mSampleId, cluster, mClusteringMethods.getChrBreakendMap(), mCopyNumberAnalyser);
        }

        if(runAnnotation(mConfig.RequiredAnnotations, REPLICATION_REPAIR))
        {
            reportClusterRepRepairSegments(mSampleId, cluster);
        }

    }

    public void logStats()
    {
        mPcClustering.logStats();
        mPcChaining.logStats();
        mPcAnnotation.logStats();
    }

}
