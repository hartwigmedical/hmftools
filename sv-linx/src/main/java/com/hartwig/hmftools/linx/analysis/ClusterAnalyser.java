package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.DOUBLE_MINUTES;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.FOLDBACK_MATCHES;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.REPLICATION_REPAIR;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.annotateChainedClusters;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.annotateTemplatedInsertions;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.calcNetCopyNumberChangeAcrossCluster;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.findIncompleteFoldbackCandidates;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.reportClusterRepRepairSegments;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.runAnnotation;
import static com.hartwig.hmftools.linx.analysis.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.analysis.LinkFinder.haveLinkedAssemblies;
import static com.hartwig.hmftools.linx.analysis.SvClassification.isSimpleSingleSV;
import static com.hartwig.hmftools.linx.analysis.SvClusteringMethods.CLUSTER_REASON_BE_PLOIDY_DROP;
import static com.hartwig.hmftools.linx.analysis.SvClusteringMethods.CLUSTER_REASON_COMMON_ARMS;
import static com.hartwig.hmftools.linx.analysis.SvClusteringMethods.CLUSTER_REASON_FOLDBACKS;
import static com.hartwig.hmftools.linx.analysis.SvClusteringMethods.CLUSTER_REASON_LOH_CHAIN;
import static com.hartwig.hmftools.linx.analysis.SvClusteringMethods.CLUSTER_REASON_LOOSE_OVERLAP;
import static com.hartwig.hmftools.linx.analysis.SvClusteringMethods.CLUSTER_REASON_NET_ARM_END_PLOIDY;
import static com.hartwig.hmftools.linx.analysis.SvClusteringMethods.addClusterReasons;
import static com.hartwig.hmftools.linx.analysis.SvClusteringMethods.checkClusterDuplicates;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.findCentromereBreakendIndex;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.markLineCluster;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.CN_SEG_DATA_MAP_BEFORE;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;
import static com.hartwig.hmftools.linx.types.ResolvedType.NONE;
import static com.hartwig.hmftools.linx.types.ResolvedType.SIMPLE_GRP;
import static com.hartwig.hmftools.linx.types.SvChain.CHAIN_ASSEMBLY_LINK_COUNT;
import static com.hartwig.hmftools.linx.types.SvChain.CHAIN_LENGTH;
import static com.hartwig.hmftools.linx.types.SvChain.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNONTATION_DM;
import static com.hartwig.hmftools.linx.types.SvCluster.areSpecificClusters;
import static com.hartwig.hmftools.linx.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.haveSameChrArms;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;
import static com.hartwig.hmftools.linx.types.SvaConstants.DEFAULT_CHAINING_SV_LIMIT;
import static com.hartwig.hmftools.linx.types.SvaConstants.MAX_FOLDBACK_CHAIN_LENGTH;
import static com.hartwig.hmftools.linx.types.SvaConstants.MAX_FOLDBACK_NEXT_CLUSTER_DISTANCE;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvArmGroup;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvChain;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.types.LinxConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ClusterAnalyser {

    private final LinxConfig mConfig;
    private SvClusteringMethods mClusteringMethods;
    private CnDataLoader mCnDataLoader;
    private SvGeneTranscriptCollection mGeneTransCache;
    private DoubleMinuteFinder mDmFinder;

    String mSampleId;
    List<SvCluster> mClusters;
    private ChainFinder mChainFinder;
    private LinkFinder mLinkFinder;

    private boolean mUseAllelePloidies;
    private boolean mRunValidationChecks;

    PerformanceCounter mPcClustering;
    PerformanceCounter mPcChaining;
    PerformanceCounter mPcAnnotation;

    public static int SMALL_CLUSTER_SIZE = 3;

    private static final Logger LOGGER = LogManager.getLogger(ClusterAnalyser.class);

    public ClusterAnalyser(final LinxConfig config, SvClusteringMethods clusteringMethods)
    {
        mConfig = config;
        mClusteringMethods = clusteringMethods;
        mCnDataLoader = null;
        mGeneTransCache = null;
        mClusters = Lists.newArrayList();
        mSampleId = "";
        mLinkFinder = new LinkFinder();
        mChainFinder = new ChainFinder();
        mDmFinder = new DoubleMinuteFinder();
        mChainFinder.setLogVerbose(mConfig.LogVerbose);
        mLinkFinder.setLogVerbose(mConfig.LogVerbose);
        mUseAllelePloidies = false;

        if(mConfig.hasMultipleSamples())
        {
            mChainFinder.getDiagnostics().setOutputDir(mConfig.OutputDataPath, mConfig.LogChainingMaxSize);
            mDmFinder.setOutputDir(mConfig.OutputDataPath);
        }

        mRunValidationChecks = false; // emabled in unit tests and after changes to merging-rule flow

        mPcClustering = new PerformanceCounter("Clustering");
        mPcChaining = new PerformanceCounter("Chaining");
        mPcAnnotation = new PerformanceCounter("Annotation");
    }

    public void setCnDataLoader(CnDataLoader cnAnalyser)
    {
        mCnDataLoader = cnAnalyser;
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
        mClusters.clear();
    }

    public final List<SvCluster> getClusters() { return mClusters; }

    public boolean clusterAndAnalyse()
    {
        mClusters.clear();
        mDmFinder.clear();

        mPcClustering.start();
        mClusteringMethods.clusterExcludedVariants(mClusters);
        mClusteringMethods.clusterByProximity(mClusters);
        mPcClustering.pause();

        // mark line clusters since these are exluded from most subsequent logic
        mClusters.forEach(x -> markLineCluster(x, mClusteringMethods.getProximityDistance()));

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
        markFoldbacks();

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
                mDmFinder.analyseCluster(cluster);

                setClusterResolvedState(cluster, false);
                continue;
            }

            // more complicated clusters for now
            boolean isSimple = cluster.getSvCount() <= SMALL_CLUSTER_SIZE && cluster.isConsistent() && !cluster.hasVariedPloidy();

            mLinkFinder.findAssembledLinks(cluster);
            applySvPloidyReplication(cluster);

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

            if (cluster.isResolved() && cluster.getResolvedType() != LINE)
                continue;

            // these are either already chained or no need to chain
            if (isSimpleSingleSV(cluster) || cluster.isFullyChained(false) || cluster.getSvCount() < 2)
            {
                setClusterResolvedState(cluster, true);
                continue;
            }

            cluster.dissolveLinksAndChains();

            cluster.removeReplicatedSvs();

            // look for and mark clusters has DM candidates, which can subsequently affect chaining
            mDmFinder.analyseCluster(cluster);

            applySvPloidyReplication(cluster);

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
                SvCluster newCluster = new SvCluster(mClusteringMethods.getNextClusterId());
                newCluster.addVariant(var);
                setClusterResolvedState(newCluster, true);
                mClusters.add(newCluster);
            }
        }
    }

    private void setClusterResolvedState(SvCluster cluster, boolean isFinal)
    {
        SvClassification.setClusterResolvedState(cluster, isFinal,
                mClusteringMethods.getDelCutoffLength(), mClusteringMethods.getDupCutoffLength(), mConfig.ProximityDistance);
    }

    private void applySvPloidyReplication(SvCluster cluster)
    {
        if(!cluster.hasVariedPloidy() && !cluster.requiresReplication())
            return;

        // use the relative copy number change to replicate some SVs within a cluster
        // isSpecificCluster(cluster);

        // first establish the lowest copy number change
        double clusterMinPloidy = cluster.getMinPloidy();
        double clusterMaxPloidy = cluster.getMaxPloidy();

        if(clusterMinPloidy <= 0)
        {
            LOGGER.debug("cluster({}) warning: invalid ploidy variation(min={} max={})",
                    cluster.id(), clusterMinPloidy, clusterMaxPloidy);
            return;
        }

        // check for samples with a broad range of ploidies, not just concentrated in a few SVs
        int totalReplicationCount = 0;
        double replicationFactor = 1;

        for(SvVarData var : cluster.getSVs())
        {
            int svPloidy = var.getImpliedPloidy();
            int svMultiple = (int)max(round(svPloidy / clusterMinPloidy),1);
            totalReplicationCount += svMultiple;
        }

        int replicationCap = mConfig.ChainingSvLimit > 0 ? min(mConfig.ChainingSvLimit, DEFAULT_CHAINING_SV_LIMIT) : DEFAULT_CHAINING_SV_LIMIT;
        if(totalReplicationCount > replicationCap)
        {
            LOGGER.debug("cluster({}) totalRepCount({}) vs svCount({}) with cluster ploidy(min={} min={}) will be scaled vs limit({})",
                    cluster.id(), totalReplicationCount, cluster.getSvCount(), clusterMinPloidy, clusterMaxPloidy, replicationCap);

            replicationFactor = replicationCap / (double)totalReplicationCount;
        }

        // replicate the SVs which have a higher copy number than their peers
        for(SvVarData var : cluster.getSVs())
        {
            int svPloidy = var.getImpliedPloidy();
            int maxAssemblyBreakends = var.getMaxAssembledBreakend();

            int svMultiple = (int)round(svPloidy / clusterMinPloidy);

            if(maxAssemblyBreakends > 1)
                svMultiple = max(svMultiple, maxAssemblyBreakends);

            svMultiple = max((int)round(svMultiple * replicationFactor), 1);

            if(svMultiple <= 1)
                continue;

            LOGGER.debug("cluster({}) replicating SV({}) {} times, copyNumChg({} vs min={})",
                    cluster.id(), var.posId(), svMultiple, svPloidy, clusterMinPloidy);

            var.setReplicatedCount(svMultiple);

            for(int j = 1; j < svMultiple; ++j)
            {
                SvVarData newVar = new SvVarData(var, true);
                cluster.addVariant(newVar);
            }
        }
    }

    public void applyComplexClusteringRules()
    {
        // second round of cluster merging on more complex criteria and inconsistencies:
        // merge on foldbacks on the same arm
        // merge on links between common arms
        // merge if one cluster has footprints which overlap unresolved complex SVs
        // merge clusters which resolve another's LOH DUP

        long longDelDupCutoffLength = max(mClusteringMethods.getDelCutoffLength(), mClusteringMethods.getDupCutoffLength());

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
        isSpecificCluster(cluster);
        if(mConfig.ChainingSvLimit > 0 && cluster.getSvCount(true) > mConfig.ChainingSvLimit)
        {
            LOGGER.info("sample({}) skipping large cluster({}) with SV counts: unique({}) replicated({})",
                    mSampleId, cluster.id(), cluster.getSvCount(), cluster.getSvCount(true));
            return;
        }

        cluster.getChains().clear();
        mChainFinder.initialise(cluster);
        mChainFinder.formChains(assembledLinksOnly);
        mChainFinder.addChains(cluster);

        if(!assembledLinksOnly)
            mChainFinder.getDiagnostics().diagnoseChains(mSampleId);

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

                    if (var1.getFoldbackLink(v1Start).isEmpty())
                        continue;

                    for (final SvVarData var2 : cluster2Foldbacks)
                    {
                        for (int be2 = SE_START; be2 <= SE_END; ++be2)
                        {
                            boolean v2Start = isStart(be2);

                            if (be2 == SE_END && var2.type() != BND)
                                continue;

                            if (var2.getFoldbackLink(v2Start).isEmpty())
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

                // isSpecificSV(breakend.getSV());

                SvBreakend nextBreakend = breakendList.get(i + 1);

                SvBreakend beFront = null; // the lower position for orientation +1 and vice versa
                SvBreakend beBack = null;

                int skippedAssemblies = 0;

                int j = i + 1;
                while(j < breakendList.size())
                {
                    nextBreakend = breakendList.get(j);

                    // first skip over any breakends in a DB with the initial breakend
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

        final SvCluster cluster = varEnd.getCluster();

        boolean singleSV = varEnd.equals(varStart);

        if(singleSV)
        {
            if(varStart.type() != INV)
                return false;
        }
        else
        {
            // must be same cluster
            if(varStart.getCluster() != cluster)
                return false;
        }

        // skip unclustered DELs & DUPs, reciprocal INV or reciprocal BNDs
        if(cluster.isResolved())
            return false;

        boolean beEndUsesStart = beEnd.usesStart();
        boolean beStartUsesStart = beStart.usesStart();

        String chainInfo = "0;0;0";

        if(singleSV)
        {
            // constraint is that the ends of this INV don't link to BND taking the path off this chromosome
            final SvChain chain = cluster.findChain(varEnd);

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
            if(varEnd.getReplicatedCount() != varStart.getReplicatedCount())
                return false;

            final SvChain chain1 = cluster.findChain(varEnd);
            final SvChain chain2 = cluster.findChain(varStart);

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
            LOGGER.debug("cluster({}) foldback inversion SV({}) length({})", cluster.id(), varEnd.posId(), length);
        }
        else
        {
            LOGGER.debug("cluster({}) foldback be1({}) be2({}) length({})", cluster.id(), beEnd.toString(), beStart.toString(), length);
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

                LOGGER.debug("cluster({}) foldback SV({} : {}) with own breakend({})",
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
    }

    private void reportClusterFeatures(final SvCluster cluster)
    {
        annotateChainedClusters(cluster, mClusteringMethods.getProximityDistance());

        if(runAnnotation(mConfig.RequiredAnnotations, FOLDBACK_MATCHES))
        {
            findIncompleteFoldbackCandidates(mSampleId, cluster, mClusteringMethods.getChrBreakendMap(), mCnDataLoader);
        }

        if(runAnnotation(mConfig.RequiredAnnotations, REPLICATION_REPAIR))
        {
            reportClusterRepRepairSegments(mSampleId, cluster);
        }

        if(runAnnotation(mConfig.RequiredAnnotations, DOUBLE_MINUTES) && cluster.hasAnnotation(CLUSTER_ANNONTATION_DM))
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
