package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.annotateNearestSvData;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.populateChromosomeBreakendMap;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.setSimpleVariantLengths;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.DOUBLE_MINUTES;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.FOLDBACK_MATCHES;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.annotateChainedClusters;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.annotateTemplatedInsertions;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.calcNetCopyNumberChangeAcrossCluster;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.findIncompleteFoldbackCandidates;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.runAnnotation;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CLUSTER_REASON_BE_PLOIDY_DROP;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CLUSTER_REASON_COMMON_ARMS;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CLUSTER_REASON_FOLDBACKS;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CLUSTER_REASON_LOH_CHAIN;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CLUSTER_REASON_LOOSE_OVERLAP;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CLUSTER_REASON_NET_ARM_END_PLOIDY;
import static com.hartwig.hmftools.linx.analysis.ClusteringState.CLUSTER_REASON_SATELLITE_SGL;
import static com.hartwig.hmftools.linx.analysis.SvClassification.isSimpleSingleSV;
import static com.hartwig.hmftools.linx.analysis.SimpleClustering.addClusterReasons;
import static com.hartwig.hmftools.linx.analysis.SimpleClustering.checkClusterDuplicates;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.findCentromereBreakendIndex;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.markLineCluster;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.CN_SEG_DATA_MAP_BEFORE;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;
import static com.hartwig.hmftools.linx.types.ResolvedType.NONE;
import static com.hartwig.hmftools.linx.types.ResolvedType.SIMPLE_GRP;
import static com.hartwig.hmftools.linx.types.SvCluster.areSpecificClusters;
import static com.hartwig.hmftools.linx.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.haveSameChrArms;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;
import static com.hartwig.hmftools.linx.types.SvaConstants.MAX_FOLDBACK_NEXT_CLUSTER_DISTANCE;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.linx.chaining.ChainFinder;
import com.hartwig.hmftools.linx.chaining.LinkFinder;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvArmGroup;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.LinxConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ClusterAnalyser {

    private final LinxConfig mConfig;
    private final ClusteringState mState;

    private SvFilters mFilters;
    private SimpleClustering mClusteringMethods;

    private CnDataLoader mCnDataLoader;
    private SvGeneTranscriptCollection mGeneTransCache;
    private DoubleMinuteFinder mDmFinder;

    private String mSampleId;
    private final List<SvCluster> mClusters;
    private List<SvVarData> mAllVariants;
    private ChainFinder mChainFinder;
    private LinkFinder mLinkFinder;

    private boolean mUseAllelePloidies;
    private boolean mRunValidationChecks;

    PerformanceCounter mPcClustering;
    PerformanceCounter mPcChaining;
    PerformanceCounter mPcAnnotation;

    public static int SMALL_CLUSTER_SIZE = 3;

    private static final Logger LOGGER = LogManager.getLogger(ClusterAnalyser.class);

    public ClusterAnalyser(final LinxConfig config)
    {
        mConfig = config;
        mState = new ClusteringState();

        mFilters = new SvFilters(mState);
        mClusteringMethods = new SimpleClustering(mState);

        mCnDataLoader = null;
        mGeneTransCache = null;
        mClusters = Lists.newArrayList();
        mSampleId = "";
        mAllVariants = Lists.newArrayList();
        mLinkFinder = new LinkFinder();
        mChainFinder = new ChainFinder();
        mDmFinder = new DoubleMinuteFinder();
        mUseAllelePloidies = false;

        if(mConfig.hasMultipleSamples())
        {
            mChainFinder.getDiagnostics().setOutputDir(mConfig.OutputDataPath, mConfig.LogChainingMaxSize);
            mDmFinder.setOutputDir(mConfig.OutputDataPath);
        }

        mChainFinder.setLogVerbose(mConfig.LogVerbose);
        mLinkFinder.setLogVerbose(mConfig.LogVerbose);

        mRunValidationChecks = false; // emabled in unit tests and after changes to merging-rule flow

        mPcClustering = new PerformanceCounter("Clustering");
        mPcChaining = new PerformanceCounter("Chaining");
        mPcAnnotation = new PerformanceCounter("Annotation");
    }

    public final ClusteringState getState() { return mState; }

    public void setCnDataLoader(CnDataLoader cnAnalyser)
    {
        mCnDataLoader = cnAnalyser;
        mState.setSampleCnEventData(mCnDataLoader.getLohData(), mCnDataLoader.getHomLossData());
        // mClusteringMethods.setChrCopyNumberMap(mCnDataLoader.getChrCopyNumberMap());

        mDmFinder.setCopyNumberAnalyser(cnAnalyser);
    }

    public void setGeneCollection(final SvGeneTranscriptCollection geneTransCache)
    {
        mGeneTransCache = geneTransCache;
        mDmFinder.setGeneTransCache(geneTransCache);
    }

    public void setUseAllelePloidies(boolean toggle)
    {
        mChainFinder.setUseAllelePloidies(toggle);
        mUseAllelePloidies = toggle;
    }

    // access for unit testing
    public final ChainFinder getChainFinder() { return mChainFinder; }
    public final LinkFinder getLinkFinder() { return mLinkFinder; }

    public void setRunValidationChecks(boolean toggle) { mRunValidationChecks = toggle; }

    public void setSampleData(final String sampleId, List<SvVarData> allVariants)
    {
        mSampleId = sampleId;
        mAllVariants = allVariants;
        mClusters.clear();
        mChainFinder.setSampleId(sampleId);
    }

    public final List<SvCluster> getClusters() { return mClusters; }

    public void preClusteringPreparation()
    {
        mState.reset();

        populateChromosomeBreakendMap(mAllVariants, mState);
        mFilters.applyFilters();
        // mFilters.reportInferredSegments(mSampleId);

        annotateNearestSvData(mState.getChrBreakendMap());

        LinkFinder.findDeletionBridges(mState.getChrBreakendMap());

        setSimpleVariantLengths(mState);
    }

    public boolean clusterAndAnalyse()
    {
        mClusters.clear();
        mDmFinder.clear();

        mPcClustering.start();
        mFilters.clusterExcludedVariants(mClusters);
        mClusteringMethods.clusterByProximity(mClusters, mConfig.ProximityDistance);
        mPcClustering.pause();

        // mark line clusters since these are excluded from most subsequent logic
        mClusters.forEach(x -> markLineCluster(x, mConfig.ProximityDistance));

        if(mRunValidationChecks)
        {
            if(!mClusteringMethods.validateClustering(mClusters))
            {
                LOGGER.info("exiting with cluster-validation errors");
                return false;
            }
        }

        mPcChaining.start();
        findLimitedChains();
        mPcChaining.pause();

        mPcClustering.resume();
        mClusteringMethods.mergeClusters(mSampleId, mClusters);
        mPcClustering.pause();

        // log basic clustering details
        mClusters.stream().filter(x -> x.getSvCount() > 1).forEach(SvCluster::logDetails);

        // INVs and other SV-pairs which make foldbacks are now used in the inconsistent clustering logic
        FoldbackFinder.markFoldbacks(mState.getChrBreakendMap());

        // subclonal clusters won't be merged any further
        mClusters.forEach(x -> x.markSubclonal());

        mPcClustering.resume();
        applyComplexClusteringRules();
        mPcClustering.stop();

        mPcChaining.resume();
        findLinksAndChains();
        dissolveSimpleGroups();
        mPcChaining.stop();

        if(mRunValidationChecks)
        {
            if(!mClusteringMethods.validateClustering(mClusters) || !checkClusterDuplicates(mClusters))
            {
                LOGGER.warn("exiting with cluster-validation errors");
                return false;
            }
        }

        mPcAnnotation.start();

        // final clean-up and analysis
        for(SvCluster cluster : mClusters)
        {
            if(!cluster.isResolved() && cluster.getResolvedType() != NONE)
            {
                // any cluster with a long DEL or DUP not merged can now be marked as resolved
                if(cluster.getSvCount() == 1 && cluster.getResolvedType().isSimple())
                    cluster.setResolved(true, cluster.getResolvedType());
            }

            cluster.cacheLinkedPairs();
            cluster.buildArmClusters();

            // isSpecificCluster(cluster);

            reportClusterFeatures(cluster);
        }

        reportOtherFeatures();

        mPcAnnotation.stop();
        return true;
    }

    public void findLimitedChains()
    {
        // chain small clusters and only assembled links in larger ones
        for(SvCluster cluster : mClusters)
        {
            if(isSimpleSingleSV(cluster))
            {
                setClusterResolvedState(cluster, false);
                continue;
            }

            // more complicated clusters for now
            boolean isSimple = cluster.getSvCount() <= SMALL_CLUSTER_SIZE && cluster.isConsistent() && !cluster.hasVariedPloidy();

            mLinkFinder.findAssembledLinks(cluster);
            cluster.setPloidyReplication(mConfig.ChainingSvLimit);

            if(isSimple)
                mDmFinder.analyseCluster(cluster);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster, !isSimple);

            if(isSimple)
            {
                setClusterResolvedState(cluster, false);

                if(cluster.isFullyChained(true))
                {
                    LOGGER.debug("cluster({}) simple and consistent with {} SVs", cluster.id(), cluster.getSvCount());
                }
            }
        }
    }

    private void findLinksAndChains()
    {
        for (SvCluster cluster : mClusters)
        {
            // isSpecificCluster(cluster);

            if (cluster.isResolved() || cluster.getResolvedType() == LINE) // only simple assembly links for LINE clusters
                continue;

            // these are either already chained or no need to chain
            if (isSimpleSingleSV(cluster) || cluster.isFullyChained(false) || cluster.getSvCount() < 2)
            {
                setClusterResolvedState(cluster, true);
                continue;
            }

            cluster.dissolveLinksAndChains();

            // look for and mark clusters has DM candidates, which can subsequently affect chaining
            mDmFinder.analyseCluster(cluster, true);

            cluster.setPloidyReplication(mConfig.ChainingSvLimit);

            // no need to re-find assembled TIs

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster, false);

            setClusterResolvedState(cluster, true);
            cluster.logDetails();
        }
    }

    private void dissolveSimpleGroups()
    {
        List<SvCluster> simpleGroups = mClusters.stream()
                .filter(x -> x.getResolvedType() == SIMPLE_GRP)
                .collect(Collectors.toList());

        for(SvCluster cluster : simpleGroups)
        {
            mClusters.remove(cluster);

            for(SvVarData var : cluster.getSVs())
            {
                SvCluster newCluster = new SvCluster(mState.getNextClusterId());
                newCluster.addVariant(var);

                mDmFinder.analyseCluster(newCluster);

                setClusterResolvedState(newCluster, true);
                mClusters.add(newCluster);
            }
        }
    }

    private void setClusterResolvedState(SvCluster cluster, boolean isFinal)
    {
        SvClassification.setClusterResolvedState(cluster, isFinal,
                mState.getDelCutoffLength(), mState.getDupCutoffLength(), mConfig.ProximityDistance);
    }


    public void applyComplexClusteringRules()
    {
        // second round of cluster merging on more complex criteria and inconsistencies:
        // merge on foldbacks on the same arm
        // merge on links between common arms
        // merge if one cluster has footprints which overlap unresolved complex SVs
        // merge clusters which resolve another's LOH DUP

        long longDelDupCutoffLength = max(mState.getDelCutoffLength(), mState.getDupCutoffLength());

        // first collect the clusters for which these complex rules apply
        List<SvCluster> complexClusters = mClusters.stream()
                .filter(x -> !x.isResolved())
                .filter(x -> !x.isSubclonal())
                .collect(Collectors.toList());

        int iterations = 1;
        boolean foundMerges = true;

        while(foundMerges)
        {
            foundMerges = false;

            int index1 = 0;
            while(index1 < complexClusters.size())
            {
                SvCluster cluster1 = complexClusters.get(index1);

                if(cluster1.isResolved())
                {
                    ++index1;
                    continue;
                }

                int index2 = index1 + 1;
                while(index2 < complexClusters.size())
                {
                    SvCluster cluster2 = complexClusters.get(index2);

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

                    foundMerges = true;
                    cluster1.mergeOtherCluster(cluster2, false);

                    complexClusters.remove(index2);
                    mClusters.remove(cluster2);
                }

                ++index1;
            }

            if(mergeTraversingClusters(complexClusters))
                foundMerges = true;

            if(mergeLOHResolvingClusters(complexClusters))
                foundMerges = true;

            if(mergePloidyResolvingClusters(complexClusters))
                foundMerges = true;

            if(mergeNetCopyNumberResolvingClusters(complexClusters))
                foundMerges = true;

            if(mergeSingleSatelliteRepeats(complexClusters))
                foundMerges = true;

            ++iterations;

            if(iterations > 20)
            {
                LOGGER.warn("sample({}) reached {} iterations of clustering merging", mSampleId, iterations);
                break;
            }
       }
    }

    private void findChains(SvCluster cluster, boolean assembledLinksOnly)
    {
        if(!assembledLinksOnly)
            isSpecificCluster(cluster);

        int svCount = cluster.getSvCount();

        if(mConfig.ChainingSvLimit > 0 && svCount > mConfig.ChainingSvLimit)
        {
            LOGGER.info("sample({}) skipping large cluster({}) with SV counts: unique({}) replicated({})",
                    mSampleId, cluster.id(), cluster.getSvCount(), svCount);
            return;
        }

        cluster.getChains().clear();
        mChainFinder.initialise(cluster);
        mChainFinder.formChains(assembledLinksOnly);
        mChainFinder.addChains(cluster);

        if(!assembledLinksOnly)
            mChainFinder.getDiagnostics().diagnoseChains();

        cluster.setValidAllelePloidySegmentPerc(mChainFinder.getValidAllelePloidySegmentPerc());
        mChainFinder.clear(); // release any refs to clusters and SVs
    }

    private boolean mergePloidyResolvingClusters(List<SvCluster> clusters)
    {
        if(!mUseAllelePloidies)
            return false;

        List<SvCluster> mergedClusters = Lists.newArrayList();

        int clusterIndex = 0;
        while(clusterIndex < clusters.size())
        {
            SvCluster cluster = clusters.get(clusterIndex);

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

                List<SvBreakend> fullBreakendList = mState.getChrBreakendMap().get(entry.getKey());

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

                        SvCluster nextCluster = nextBreakend.getCluster();
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

                            addClusterReasons(breakend.getSV(), nextBreakend.getSV(), CLUSTER_REASON_BE_PLOIDY_DROP);

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

        mergedClusters.forEach(x -> clusters.remove(x));
        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }

    private boolean mergeLOHResolvingClusters(List<SvCluster> clusters)
    {
        // merge clusters if one resolves another's LOH event with a DUP on one side
        List<SvCluster> clustersWithLohEvents = clusters.stream()
                .filter(x -> !x.getLohEvents().isEmpty())
                .filter(x -> !x.hasLinkingLineElements())
                .collect(Collectors.toList());

        List<SvCluster> mergedClusters = Lists.newArrayList();

        for(SvCluster lohCluster : clustersWithLohEvents)
        {
            if(mergedClusters.contains(lohCluster)) // if this has LOH events they will have been added to the parent cluster
                continue;

            List<LohEvent> lohEvents = lohCluster.getLohEvents();

            int lohIndex = 0; // used since merging another cluster can add more LOH events
            while(lohIndex < lohEvents.size())
            {
                LohEvent lohEvent = lohEvents.get(lohIndex);

                if(!lohEvent.isValid())
                {
                    ++lohIndex;
                    continue;
                }

                for(int be = SE_START; be <= SE_END; ++be)
                {
                    SvBreakend lohBreakend = lohEvent.getBreakend(isStart(be));

                    if(lohBreakend == null || lohBreakend.getSV().type() != DUP)
                        continue;

                    // it's possible that the breakends for this LOH are not clustered, eg if one is LINE
                    if(lohBreakend.getCluster() != lohCluster)
                        continue;

                    // walk towards the LOH from the other end of this DUP to see if it can find a resolving event within the cluster
                    List<SvBreakend> fullBreakendList = mState.getChrBreakendMap().get(lohBreakend.chromosome());

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

                            addClusterReasons(lohBreakend.getSV(), resolvingBreakend.getSV(), CLUSTER_REASON_LOH_CHAIN);

                            resolvingCluster.addClusterReason(CLUSTER_REASON_LOH_CHAIN);
                            lohCluster.addClusterReason(CLUSTER_REASON_LOH_CHAIN);

                            lohCluster.mergeOtherCluster(resolvingCluster);

                            mergedClusters.add(resolvingCluster);
                            break;
                        }

                        if(nextBreakend.orientation() == otherBreakend.orientation())
                            continue;

                        SvCluster otherCluster = nextBreakend.getCluster();

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

        mergedClusters.forEach(x -> clusters.remove(x));
        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }

    private boolean mergeNetCopyNumberResolvingClusters(List<SvCluster> clusters)
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
        while(clusterIndex < clusters.size())
        {
            SvCluster cluster = clusters.get(clusterIndex);

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
                List<SvCNData> cnDataList = mCnDataLoader.getChrCnDataMap().get(chromosome);

                if(cnDataList == null || cnDataList.isEmpty())
                    continue;

                List<SvBreakend> fullBreakendList = mState.getChrBreakendMap().get(entry.getKey());

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

                    double[] centromereCNData = mCnDataLoader.getCentromereCopyNumberData(chromosome, arm.equals(CHROMOSOME_ARM_P));

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

                            SvCluster otherCluster = nextBreakend.getCluster();

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

                                addClusterReasons(clusterBreakend.getSV(), nextBreakend.getSV(), CLUSTER_REASON_NET_ARM_END_PLOIDY);

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

        mergedClusters.forEach(x -> clusters.remove(x));
        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }

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
                for (int be1 = SE_START; be1 <= SE_END; ++be1)
                {
                    boolean v1Start = isStart(be1);

                    if (be1 == SE_END && var1.type() != BND)
                        continue;

                    if (var1.getFoldbackBreakend(v1Start) == null)
                        continue;

                    for (final SvVarData var2 : cluster2Foldbacks)
                    {
                        for (int be2 = SE_START; be2 <= SE_END; ++be2)
                        {
                            boolean v2Start = isStart(be2);

                            if (be2 == SE_END && var2.type() != BND)
                                continue;

                            if (var2.getFoldbackBreakend(v2Start) == null)
                                continue;

                            if (!var1.chromosome(v1Start).equals(var2.chromosome(v2Start)) || !var1.arm(v1Start).equals(var2.arm(v2Start)))
                                continue;

                            LOGGER.debug("cluster({}) SV({}) and cluster({}) SV({}) have foldbacks on same arm",
                                    cluster1.id(), var1.posId(), cluster2.id(), var2.posId());

                            addClusterReasons(var1, var2, CLUSTER_REASON_FOLDBACKS);

                            cluster1.addClusterReason(CLUSTER_REASON_FOLDBACKS);
                            cluster2.addClusterReason(CLUSTER_REASON_FOLDBACKS);
                            return true;
                        }
                    }
                }
            }
        }

        final Map<String, List<SvBreakend>> chrBreakendMap = mState.getChrBreakendMap();

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
                || nextBreakend.getCluster() != otherCluster)
                    continue;

                if (abs(nextBreakend.position() - foldbackBreakend.position()) > MAX_FOLDBACK_NEXT_CLUSTER_DISTANCE)
                    continue;

                double fbPloidy = foldbackBreakend.ploidy();
                double nbPloidy = nextBreakend.ploidy();

                if (nbPloidy < fbPloidy && !copyNumbersEqual(nbPloidy, fbPloidy))
                    continue;

                LOGGER.debug("cluster({}) foldback breakend({}) faces cluster({}) breakend({})",
                        foldbackCluster.id(), foldbackBreakend.toString(), otherCluster.id(), nextBreakend.toString());

                addClusterReasons(var, nextBreakend.getSV(), CLUSTER_REASON_FOLDBACKS);

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
        final SvCluster fbCluster = foldbackBreakend.getCluster();

        int index = startIndex;

        while(index >= 0 && index < breakendList.size())
        {
            final SvBreakend breakend = breakendList.get(index);
            final SvCluster cluster = breakend.getCluster();

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

                        addClusterReasons(var1, var2, CLUSTER_REASON_COMMON_ARMS);

                        cluster1.addClusterReason(CLUSTER_REASON_COMMON_ARMS);
                        cluster2.addClusterReason(CLUSTER_REASON_COMMON_ARMS);
                        return true;
                    }
                }
            }
        }

        return false;
    }

    private boolean mergeTraversingClusters(List<SvCluster> clusters)
    {
        // merge clusters if one is straddled by opposing breakends in the other
        List<SvCluster> mergedClusters = Lists.newArrayList();

        int clusterIndex = 0;
        while(clusterIndex < clusters.size())
        {
            SvCluster cluster = clusters.get(clusterIndex);

            if(mergedClusters.contains(cluster) || cluster.isResolved())
            {
                ++clusterIndex;
                continue;
            }

            boolean mergedOtherClusters = false;

            for (final Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
            {
                List<SvBreakend> breakendList = entry.getValue();

                List<SvBreakend> fullBreakendList = mState.getChrBreakendMap().get(entry.getKey());

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

                            if (lowerBreakend.getSV().getFoldbackBreakend(lowerBreakend.usesStart()) != null
                            && lowerBreakend.getSV().getFoldbackBreakend(lowerBreakend.usesStart())
                                    == nextBreakend.getSV().getFoldbackBreakend(nextBreakend.usesStart()))
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

                        final SvCluster otherCluster = otherBreakend.getCluster();

                        if (otherCluster == cluster || otherCluster.isResolved() || mergedClusters.contains(otherCluster))
                            continue;

                        LOGGER.debug("cluster({}) breakends({} & {}) overlap cluster({}) breakend({})",
                                cluster.id(), lowerBreakend.toString(), upperBreakend.toString(), otherCluster.id(), otherBreakend.toString());

                        addClusterReasons(otherBreakend.getSV(), lowerBreakend.getSV(), CLUSTER_REASON_LOOSE_OVERLAP);

                        otherCluster.addClusterReason(CLUSTER_REASON_LOOSE_OVERLAP);
                        cluster.addClusterReason(CLUSTER_REASON_LOOSE_OVERLAP);
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

        mergedClusters.forEach(x -> clusters.remove(x));
        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }

    private boolean mergeSingleSatelliteRepeats(List<SvCluster> clusters)
    {
        // merge any cluster with less than or 1 non SGL and non INF breakend with any other cluster which contains a SGL on the same
        // chromosome with matching repeat class or type marked as satellite

        // To protect against false positives and joining complex clusters which both touch repeats, but otherwise don’t appear to overlap,
        // we avoid clustering 2 clusters which already have multiple non SGL breakends
        List<SvCluster> clustersWithSatelliteRepeats = clusters.stream()
                .filter(x -> x.getSVs().stream().anyMatch(y -> y.sglToSatelliteRepeats()))
                .collect(Collectors.toList());

        if(clustersWithSatelliteRepeats.isEmpty())
            return false;

        List<SvCluster> clustersWithOneOrNoSgls = clustersWithSatelliteRepeats.stream()
                .filter(x -> x.getTypeCount(SGL) <= 1)
                .collect(Collectors.toList());

        List<SvCluster> mergedClusters = Lists.newArrayList();

        for(SvCluster srCluster : clustersWithSatelliteRepeats)
        {
            if(mergedClusters.contains(srCluster))
                continue;

            List<String> satelliteChromosomes1 = srCluster.getSVs().stream()
                    .filter(x -> x.sglToSatelliteRepeats())
                    .map(x -> x.chromosome(true))
                    .collect(Collectors.toList());

            int index = 0;
            while(index < clustersWithOneOrNoSgls.size())
            {
                SvCluster sglCluster = clustersWithOneOrNoSgls.get(index);

                if(sglCluster == srCluster)
                {
                    ++index;
                    continue;
                }


                boolean merged = false;

                for(SvVarData var : sglCluster.getSVs())
                {
                    if(var.sglToSatelliteRepeats() && satelliteChromosomes1.contains(var.chromosome(true)))
                    {
                        final String chromosome = var.chromosome(true);
                        LOGGER.debug("cluster({}) has same chromosome({}) link with satellite cluster({}) SV({})",
                                srCluster.id(), chromosome, sglCluster.id(), var.id());

                        // find the other linking SGL
                        final SvVarData otherSV = srCluster.getSVs().stream()
                                .filter(x -> x.sglToSatelliteRepeats() && x.chromosome(true).equals(chromosome))
                                .findFirst().get();

                        addClusterReasons(otherSV, var, CLUSTER_REASON_SATELLITE_SGL);
                        srCluster.addClusterReason(CLUSTER_REASON_SATELLITE_SGL);

                        srCluster.mergeOtherCluster(sglCluster);
                        mergedClusters.add(sglCluster);
                        merged = true;
                        break;
                    }
                }

                if(merged)
                    clustersWithOneOrNoSgls.remove(index);
                else
                    ++index;
            }
        }

        if(mergedClusters.isEmpty())
            return false;

        mergedClusters.forEach(x -> clusters.remove(x));
        mergedClusters.forEach(x -> mClusters.remove(x));
        return true;
    }


    private void reportOtherFeatures()
    {
        annotateTemplatedInsertions(mClusters, mState.getChrBreakendMap());
        // checkSkippedLOHEvents();

        // annotateFoldbacks(mClusters); // unused for now
    }

    private void reportClusterFeatures(final SvCluster cluster)
    {
        annotateChainedClusters(cluster, mConfig.ProximityDistance);

        if(runAnnotation(mConfig.RequiredAnnotations, FOLDBACK_MATCHES))
        {
            findIncompleteFoldbackCandidates(mSampleId, cluster, mState.getChrBreakendMap(), mCnDataLoader);
        }

        if(runAnnotation(mConfig.RequiredAnnotations, DOUBLE_MINUTES))
            mDmFinder.reportCluster(mSampleId, cluster);
    }

    public void close()
    {
        mDmFinder.close();
        mChainFinder.getDiagnostics().close();
    }

    public void logStats()
    {
        mPcClustering.logStats();
        mPcChaining.logStats();
        mPcAnnotation.logStats();
    }

}
