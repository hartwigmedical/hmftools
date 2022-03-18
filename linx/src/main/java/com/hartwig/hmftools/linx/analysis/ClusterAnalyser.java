package com.hartwig.hmftools.linx.analysis;

import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.UNDER_CLUSTERING;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.annotateClusterChains;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.annotateClusterDeletions;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.annotateReplicationBeforeRepair;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.annotateTemplatedInsertions;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.reportUnderclustering;
import static com.hartwig.hmftools.linx.analysis.ClusterAnnotations.runAnnotation;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.isFilteredResolvedType;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.isSimpleSingleSV;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.annotateNearestSvData;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.associateBreakendCnEvents;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.populateChromosomeBreakendMap;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.setSimpleVariantLengths;
import static com.hartwig.hmftools.linx.analysis.SimpleClustering.checkClusterDuplicates;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.DELETED_TOTAL;
import static com.hartwig.hmftools.linx.chaining.ChainJcnLimits.RANGE_TOTAL;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.createAssemblyLinkedPairs;
import static com.hartwig.hmftools.linx.types.ArmCluster.buildArmClusters;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;
import static com.hartwig.hmftools.linx.types.ResolvedType.NONE;
import static com.hartwig.hmftools.linx.types.ResolvedType.SIMPLE_GRP;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_DM;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.linx.CohortDataWriter;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.annotators.LineElementAnnotator;
import com.hartwig.hmftools.linx.chaining.ChainFinder;
import com.hartwig.hmftools.linx.chaining.LinkFinder;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class ClusterAnalyser {

    private final LinxConfig mConfig;
    private final ClusteringState mState;

    private final SvFiltering mFilters;
    private final SimpleClustering mSimpleClustering;
    private final ComplexClustering mComplexClustering;

    private CnDataLoader mCnDataLoader;
    private final DoubleMinuteFinder mDmFinder;
    private final BfbFinder mBfbFinder;
    private LineElementAnnotator mLineElementAnnotator;

    private String mSampleId;
    private final List<SvCluster> mClusters;
    private final List<SvVarData> mAllVariants;
    private final ChainFinder mChainFinder;

    private boolean mRunValidationChecks;

    PerformanceCounter mPcClustering;
    PerformanceCounter mPcChaining;

    private static final int SMALL_CLUSTER_SIZE = 3;

    public ClusterAnalyser(final LinxConfig config, final CohortDataWriter cohortDataWriter)
    {
        mConfig = config;
        mState = new ClusteringState();
        mClusters = Lists.newArrayList();

        mFilters = new SvFiltering(mState);
        mSimpleClustering = new SimpleClustering(mConfig, mState, cohortDataWriter);
        mComplexClustering = new ComplexClustering(mState, mClusters, mSimpleClustering);

        mCnDataLoader = null;
        mLineElementAnnotator = null;
        mSampleId = "";
        mAllVariants = Lists.newArrayList();
        mChainFinder = new ChainFinder(cohortDataWriter);
        mDmFinder = new DoubleMinuteFinder(config, cohortDataWriter, mState.getChrBreakendMap());
        mBfbFinder = new BfbFinder();

        if(mConfig.hasMultipleSamples())
            mChainFinder.initialiseOutput(mConfig);

        mChainFinder.setUseAllelePloidies(true); // can probably remove and assume always in place
        mChainFinder.setLogVerbose(mConfig.LogVerbose);

        mRunValidationChecks = false; // enabled in unit tests and after changes to merging-rule flow

        mPcClustering = new PerformanceCounter("Clustering");
        mPcChaining = new PerformanceCounter("Chaining");
    }

    public final ClusteringState getState() { return mState; }

    public void setLineAnnotator(final LineElementAnnotator lineAnnotator)
    {
        mLineElementAnnotator = lineAnnotator;
    }

    public void setCnDataLoader(CnDataLoader cnDataLoader)
    {
        mCnDataLoader = cnDataLoader;
        mState.setSampleCnEventData(mCnDataLoader.getLohData(), mCnDataLoader.getHomLossData());

        mDmFinder.setCopyNumberAnalyser(cnDataLoader);
        mBfbFinder.setCopyNumberAnalyser(cnDataLoader);
    }

    public void setGeneCollection(final EnsemblDataCache geneDataCache)
    {
        mDmFinder.setGeneTransCache(geneDataCache);
    }

    // access for unit testing
    public final ChainFinder getChainFinder() { return mChainFinder; }
    public final DoubleMinuteFinder getDoubleMinuteFinder() { return mDmFinder; }

    public void setRunValidationChecks(boolean toggle) { mRunValidationChecks = toggle; }

    public void setSampleData(final String sampleId, List<SvVarData> allVariants)
    {
        mSampleId = sampleId;
        mAllVariants.clear();
        mAllVariants.addAll(allVariants);
        mClusters.clear();
        mSimpleClustering.initialise(sampleId);
        mChainFinder.setSampleId(sampleId);
    }

    public final List<SvCluster> getClusters() { return mClusters; }

    public void preClusteringPreparation()
    {
        mState.reset();

        populateChromosomeBreakendMap(mAllVariants, mState);
        mFilters.applyFilters();

        annotateNearestSvData(mState.getChrBreakendMap());

        LinkFinder.findDeletionBridges(mState.getChrBreakendMap());

        setSimpleVariantLengths(mState);
    }

    public boolean clusterAndAnalyse()
    {
        if(mConfig.IsGermline)
            return clusterAndAnalyseGermline();

        mClusters.clear();
        mDmFinder.clear();

        mPcClustering.start();
        mFilters.clusterExcludedVariants(mClusters);
        mSimpleClustering.clusterByProximity(mClusters);
        mClusters.stream().filter(x -> x.getSvCount() > 1).forEach(x -> x.updateClusterDetails());
        mPcClustering.pause();

        // mark line clusters since these are excluded from most subsequent logic
        mClusters.forEach(x -> mLineElementAnnotator.markLineCluster(x));

        associateBreakendCnEvents(mSampleId, mState);

        if(mRunValidationChecks)
        {
            if(!mSimpleClustering.validateClustering(mClusters))
            {
                LNX_LOGGER.info("exiting with cluster-validation errors");
                return false;
            }
        }

        mPcChaining.start();
        findLimitedChains();
        mPcChaining.pause();

        mPcClustering.resume();
        mSimpleClustering.mergeClusters(mClusters);
        mPcClustering.pause();

        // log basic clustering details
        mClusters.stream().filter(x -> x.getSvCount() > 1).forEach(SvCluster::logDetails);

        // INVs and other SV-pairs which make foldbacks are now used in the inconsistent clustering logic
        FoldbackFinder.markFoldbacks(mState.getChrBreakendMap());

        mPcClustering.resume();
        mComplexClustering.applyRules(mSampleId, false);
        mSimpleClustering.mergeLongDelDupClusters(mClusters);
        mPcClustering.stop();

        mPcChaining.resume();
        findLinksAndChains();
        dissolveSimpleGroups();
        mPcChaining.stop();

        if(mRunValidationChecks)
        {
            if(!mSimpleClustering.validateClustering(mClusters) || !checkClusterDuplicates(mClusters))
            {
                LNX_LOGGER.warn("exiting with cluster-validation errors");
                return false;
            }
        }

        // final clean-up and analysis

        // take note of DM clusters so they can be rechecked after any final cluster merging
        final List<Integer> dmClusterIds = mClusters.stream().filter(x -> x.hasAnnotation(CLUSTER_ANNOT_DM))
                .map(x -> Integer.valueOf(x.id())).collect(Collectors.toList());

        // re-check foldbacks amongst newly formed chains and then DM status
        if(FoldbackFinder.markFoldbacks(mState.getChrBreakendMap(), true))
        {
            mComplexClustering.applyRules(mSampleId, true);
        }

        for(SvCluster cluster : mClusters)
        {
            if(isFilteredResolvedType(cluster.getResolvedType()))
                continue;

            if(!cluster.getDoubleMinuteSVs().isEmpty() || dmClusterIds.contains(cluster.id()))
                mDmFinder.analyseCluster(cluster, true);

            if(!cluster.hasAnnotation(CLUSTER_ANNOT_DM) && !cluster.isResolved())
                mBfbFinder.analyseCluster(cluster);

            if(!cluster.isResolved())
            {
                if(cluster.getResolvedType() != NONE)
                {
                    // any cluster with a long DEL or DUP not merged can now be marked as resolved
                    if(cluster.getSvCount() == 1 && cluster.getResolvedType().isSimple())
                        cluster.setResolved(true, cluster.getResolvedType());
                }
                else
                {
                    setClusterResolvedState(cluster, true);
                    cluster.logDetails();
                }
            }

            cluster.cacheLinkedPairs();
            buildArmClusters(cluster);
        }

        return true;
    }

    public boolean clusterAndAnalyseGermline()
    {
        mClusters.clear();

        mPcClustering.start();
        mFilters.clusterExcludedVariants(mClusters);
        mSimpleClustering.clusterByProximity(mClusters);
        mClusters.stream().filter(x -> x.getSvCount() > 1).forEach(x -> x.updateClusterDetails());
        mPcClustering.pause();

        // mark line clusters since these are excluded from most subsequent logic
        mClusters.forEach(x -> mLineElementAnnotator.markLineCluster(x));

        mPcChaining.start();
        findLimitedChains();

        // log basic clustering details
        mClusters.stream().filter(x -> x.getSvCount() > 1).forEach(SvCluster::logDetails);

        FoldbackFinder.markFoldbacks(mState.getChrBreakendMap());

        dissolveSimpleGroups();
        mPcChaining.stop();

        for(SvCluster cluster : mClusters)
        {
            if(!cluster.isResolved())
            {
                if(cluster.getResolvedType() != NONE)
                {
                    // any cluster with a long DEL or DUP not merged can now be marked as resolved
                    if(cluster.getSvCount() == 1 && cluster.getResolvedType().isSimple())
                        cluster.setResolved(true, cluster.getResolvedType());
                }
                else
                {
                    setClusterResolvedState(cluster, true);
                    cluster.logDetails();
                }
            }

            cluster.cacheLinkedPairs();
            buildArmClusters(cluster);
        }

        return true;
    }

    private void findLimitedChains()
    {
        // chain small clusters and only assembled links in larger ones
        boolean checkDMs = !mConfig.IsGermline;

        for(SvCluster cluster : mClusters)
        {
            if(isFilteredResolvedType(cluster.getResolvedType()))
                continue;

            if(checkDMs && isSimpleSingleSV(cluster))
            {
                mDmFinder.analyseCluster(cluster);
                setClusterResolvedState(cluster, false);
                continue;
            }

            // more complicated clusters for now
            boolean isSimple = cluster.getSvCount() <= SMALL_CLUSTER_SIZE && cluster.isConsistent() && !cluster.hasVariedJcn();

            cluster.setAssemblyLinkedPairs(createAssemblyLinkedPairs(cluster));
            cluster.determineRequiresReplication();

            if(checkDMs)
                mDmFinder.analyseCluster(cluster);

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster, !isSimple);

            if(isSimple)
            {
                setClusterResolvedState(cluster, false);

                if(cluster.isFullyChained(true))
                {
                    LNX_LOGGER.debug("cluster({}) simple and consistent with {} SVs", cluster.id(), cluster.getSvCount());
                }
            }
        }
    }

    private void findLinksAndChains()
    {
        for (SvCluster cluster : mClusters)
        {
            if (cluster.getResolvedType() == LINE) // only simple assembly links for LINE clusters
                continue;

            if (cluster.getResolvedType() != NONE) // any cluster previously resolved and not modified does not need to be chained again
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

            cluster.determineRequiresReplication();

            // no need to re-find assembled TIs

            // then look for fully-linked clusters, ie chains involving all SVs
            findChains(cluster, false);

            setClusterResolvedState(cluster, true);
            cluster.logDetails();
        }
    }

    private void dissolveSimpleGroups()
    {
        // break apart any clusters of simple SVs which aren't likely or required to be chained
        List<SvCluster> simpleGroups = mClusters.stream().filter(x -> x.getResolvedType() == SIMPLE_GRP).collect(Collectors.toList());

        final List<SvVarData> discardSVs = Lists.newArrayList();

        for(SvCluster cluster : simpleGroups)
        {
            LNX_LOGGER.debug("cluster({}: {}) de-merging {} simple SVs", cluster.id(), cluster.getDesc(), discardSVs.size());
            mClusters.remove(cluster);
            discardSVs.addAll(cluster.getSVs());
        }

        for(SvVarData var : discardSVs)
        {
            SvCluster newCluster = new SvCluster(mState.getNextClusterId());
            var.clearClusteringData();
            newCluster.addVariant(var);

            mDmFinder.analyseCluster(newCluster);

            setClusterResolvedState(newCluster, true);
            mClusters.add(newCluster);
        }
    }

    private void setClusterResolvedState(SvCluster cluster, boolean isFinal)
    {
        ClusterClassification.setClusterResolvedState(cluster, isFinal, mConfig.IsGermline,
                mState.getDelCutoffLength(), mState.getDupCutoffLength(), mState.getChrBreakendMap());
    }

    private void findChains(SvCluster cluster, boolean assembledLinksOnly)
    {
        if(mConfig.ChainingSvLimit > 0 && cluster.getSvCount() > mConfig.ChainingSvLimit)
        {
            LNX_LOGGER.debug("sample({}) skipping chaining large cluster({}) with SV count({})",
                    mSampleId, cluster.id(), cluster.getSvCount());
            return;
        }

        cluster.getChains().clear();
        mChainFinder.initialise(cluster);
        mChainFinder.formChains(assembledLinksOnly);
        mChainFinder.addChains(cluster);

        if(!assembledLinksOnly)
            mChainFinder.getDiagnostics().diagnoseChains();

        final long[] rangeData = mChainFinder.calcRangeData();

        if(rangeData != null)
        {
            cluster.getMetrics().ValidAlleleJcnSegmentPerc = mChainFinder.getValidAllelePloidySegmentPerc();
            cluster.getMetrics().TraversedRange = rangeData[RANGE_TOTAL];
            cluster.getMetrics().TotalDeleted = rangeData[DELETED_TOTAL];
        }

        mChainFinder.clear(); // release any refs to clusters and SVs
    }

    public void annotateClusters()
    {
        // final clean-up and analysis
        mClusters.forEach(x -> annotateTemplatedInsertions(x, mState.getChrBreakendMap()));

        mClusters.forEach(this::reportClusterFeatures);

        if(runAnnotation(mConfig.RequiredAnnotations, UNDER_CLUSTERING))
        {
            reportUnderclustering(mSampleId, mClusters, mState.getChrBreakendMap());
        }
    }

    private void reportClusterFeatures(final SvCluster cluster)
    {
        annotateClusterChains(cluster);
        annotateClusterDeletions(cluster, mState.getChrBreakendMap());
        annotateReplicationBeforeRepair(cluster);

        ClusterMetrics metrics = cluster.getMetrics();

        if(metrics.TotalDeleted == 0)
            metrics.TotalDeleted = metrics.TotalDBLength;

        mDmFinder.reportCluster(mSampleId, cluster);
    }

    public void close()
    {
        mDmFinder.close();
    }

    public void logStats()
    {
        mPcClustering.logStats();
        mPcChaining.logStats();
    }

}
