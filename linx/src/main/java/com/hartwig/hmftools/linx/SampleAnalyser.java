package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.purple.Gender.MALE;
import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.SvFileLoader.createSvData;
import static com.hartwig.hmftools.linx.SvFileLoader.loadVariantsFromVcf;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.getClusterCategory;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.linkSglMappedInferreds;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.asStr;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.linx.LinxGermlineSv;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.linx.ImmutableLinxCluster;
import com.hartwig.hmftools.common.linx.ImmutableLinxLink;
import com.hartwig.hmftools.common.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxLink;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.linx.analysis.ClusterAnalyser;
import com.hartwig.hmftools.linx.annotators.LineElementType;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.CnSegmentBuilder;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.drivers.DriverGeneAnnotator;
import com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser;
import com.hartwig.hmftools.linx.fusion.FusionResources;
import com.hartwig.hmftools.linx.gene.BreakendGenePrep;
import com.hartwig.hmftools.linx.types.ArmCluster;
import com.hartwig.hmftools.common.purple.ChromosomeArm;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumber;
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomain;
import com.hartwig.hmftools.linx.visualiser.file.VisSampleData;
import com.hartwig.hmftools.linx.visualiser.file.VisSegment;
import com.hartwig.hmftools.linx.visualiser.file.VisSvData;

public class SampleAnalyser implements Callable
{
    private final int mId;
    private final LinxConfig mConfig;
    private List<String> mSampleIds;

    private final ClusterAnalyser mAnalyser;
    private final EnsemblDataCache mEnsemblDataCache;

    private final VisSampleData mVisSampleData;
    private final CohortDataWriter mCohortDataWriter;

    private final SvAnnotators mSvAnnotators;
    private final CnDataLoader mCnDataLoader;

    // data per run (ie sample)
    private String mCurrentSampleId;
    private final List<SvVarData> mAllVariants; // full list of SVs

    private final DriverGeneAnnotator mDriverGeneAnnotator;
    private final FusionDisruptionAnalyser mFusionAnalyser;

    private boolean mIsValid;

    private final Map<String,PerformanceCounter> mPerfCounters;

    public static final String PERF_COUNTER_TOTAL = "Total";
    public static final String PERF_COUNTER_PREP = "Preparation";
    public static final String PERF_COUNTER_CLUSTER = "ClusterAndAnalyse";
    public static final String PERF_COUNTER_ANNOTATE = "Annotation";
    public static final String PERF_COUNTER_WRITE = "WriteAndUpload";

    public SampleAnalyser(
            int instanceId, final LinxConfig config, final SvAnnotators svAnnotators,
            final EnsemblDataCache ensemblDataCache, final FusionResources fusionResources, final CohortDataWriter cohortDataWriter)
    {
        mId = instanceId;
        mConfig = config;
        mSampleIds = Lists.newArrayList();

        mVisSampleData = new VisSampleData();
        mCurrentSampleId = "";

        mEnsemblDataCache = ensemblDataCache;

        mCohortDataWriter = cohortDataWriter;
        mSvAnnotators = svAnnotators;

        mAnalyser = new ClusterAnalyser(config, mCohortDataWriter);

        mCnDataLoader = new CnDataLoader(config.PurpleDataPath);

        mDriverGeneAnnotator = mConfig.RunDrivers ?
                new DriverGeneAnnotator(ensemblDataCache, config, mCnDataLoader, cohortDataWriter, mVisSampleData) : null;

        mFusionAnalyser = new FusionDisruptionAnalyser(config, ensemblDataCache, fusionResources, cohortDataWriter, mVisSampleData);

        mAllVariants = Lists.newArrayList();

        mIsValid = true;

        mAnalyser.setLineAnnotator(mSvAnnotators.LineElementAnnotator);
        mAnalyser.setCnDataLoader(mCnDataLoader);
        mAnalyser.setGeneCollection(ensemblDataCache);

        mPerfCounters = Maps.newHashMap();
        mPerfCounters.put(PERF_COUNTER_TOTAL, new PerformanceCounter(PERF_COUNTER_TOTAL));
        mPerfCounters.put(PERF_COUNTER_PREP, new PerformanceCounter(PERF_COUNTER_PREP));
        mPerfCounters.put(PERF_COUNTER_CLUSTER, new PerformanceCounter(PERF_COUNTER_CLUSTER));
        mPerfCounters.put(PERF_COUNTER_WRITE, new PerformanceCounter(PERF_COUNTER_WRITE));
        mPerfCounters.put(PERF_COUNTER_ANNOTATE, new PerformanceCounter(PERF_COUNTER_ANNOTATE));

        if(mConfig.RunFusions)
        {
            PerformanceCounter fusionPc = mFusionAnalyser.getPerfCounter();
            mPerfCounters.put(fusionPc.getName(), fusionPc);
        }

        if(mDriverGeneAnnotator != null)
        {
            PerformanceCounter driverPc = mDriverGeneAnnotator.getPerfCounter();
            mPerfCounters.put(driverPc.getName(), driverPc);
        }
    }

    public Map<String,PerformanceCounter> getPerfCounters() { return mPerfCounters; }

    public void setSampleIds(final List<String> sampleIds)
    {
        mSampleIds.clear();
        mSampleIds.addAll(sampleIds);
    }

    @Override
    public Long call()
    {
        processSamples();
        return (long)1;
    }

    public void processSamples()
    {
        if(mSampleIds.size() == 1)
        {
            mPerfCounters.get(PERF_COUNTER_TOTAL).start();
            processSample(mSampleIds.get(0));
            mPerfCounters.get(PERF_COUNTER_TOTAL).stop();
            return;
        }

        LNX_LOGGER.info("{}: processing {} samples", mId, mSampleIds.size());

        for(int i = 0; i < mSampleIds.size(); ++i)
        {
            mPerfCounters.get(PERF_COUNTER_TOTAL).start();

            try
            {
                processSample(mSampleIds.get(i));
            }
            catch(Exception e)
            {
                LNX_LOGGER.error("sample({}) processing failed: {}", mSampleIds.get(i), e.toString());

                if(mConfig.FailOnMissing || mConfig.isSingleSample())
                {
                    e.printStackTrace();
                    System.exit(1);
                }
            }

            if(i > 10 && (i % 10) == 0)
            {
                LNX_LOGGER.info("{}: processed {} samples", mId, i);
            }

            mPerfCounters.get(PERF_COUNTER_TOTAL).stop();

            if(!inValidState())
                break;
        }

        LNX_LOGGER.info("{}: completed processing of {} samples", mId, mSampleIds.size());
    }

    private void processSample(final String sampleId)
    {
        mCurrentSampleId = sampleId;
        mVisSampleData.setSampleId(sampleId);

        List<StructuralVariantData> svRecords = loadVariantsFromVcf(mConfig, mCurrentSampleId);

        final List<SvVarData> svDataList = createSvData(svRecords, mConfig);

        if(svDataList.isEmpty())
        {
            LNX_LOGGER.info("sample({}) has no passing SVs", mCurrentSampleId);

            if(mConfig.isSingleSample())
                writeSampleWithNoSVs();

            return;
        }

        if(!mConfig.IsGermline)
            mCnDataLoader.loadSampleData(mCurrentSampleId, svRecords);

        setSampleSVs(svDataList);

        if(mEnsemblDataCache != null)
        {
            Map<String,Integer> specificPreGeneDistances = mFusionAnalyser.getSpecialFusions().getSpecificPreGeneDistances();
            int preGeneDefaultDistance = mConfig.RunFusions ? PRE_GENE_PROMOTOR_DISTANCE : 0;

            BreakendGenePrep.setSvGeneData(mAllVariants, mEnsemblDataCache, preGeneDefaultDistance, specificPreGeneDistances);
        }

        analyse();

        if(!inValidState())
        {
            LNX_LOGGER.info("exiting after sample({}), in invalid state", mCurrentSampleId);
            return;
        }

        // when matching RNA, allow all transcripts regardless of their viability for fusions
        boolean purgeInvalidTranscripts = !mFusionAnalyser.hasRnaSampleData();
        mFusionAnalyser.annotateTranscripts(mAllVariants, purgeInvalidTranscripts);

        annotate();

        if(mDriverGeneAnnotator != null)
            mDriverGeneAnnotator.annotateSVs(mCurrentSampleId, getChrBreakendMap());

        if(mConfig.RunFusions || mConfig.IsGermline)
            mFusionAnalyser.run(mCurrentSampleId, mAllVariants, mAnalyser.getClusters(), getChrBreakendMap());

        writeOutput();
        close();

        LNX_LOGGER.info("sample({}) procesed {} SVs", sampleId, svDataList.size());
    }

    public final List<SvVarData> getVariants() { return mAllVariants; }

    public boolean inValidState() { return mIsValid; }
    public final Map<String, List<SvBreakend>> getChrBreakendMap() { return mAnalyser.getState().getChrBreakendMap(); }

    public void setSampleSVs(final List<SvVarData> variants)
    {
        mAllVariants.clear();

        if(variants.isEmpty())
            return;

        mAllVariants.addAll(variants);

        linkSglMappedInferreds(mAllVariants);

        if(!mConfig.IsGermline)
        {
            // look-up and cache relevant CN data into each SV
            SvCNData.setSvCopyNumberData(
                    mAllVariants,
                    mCnDataLoader.getSvJcnCalcMap(),
                    mCnDataLoader.getSvIdCnDataMap(),
                    mCnDataLoader.getChrCnDataMap());
        }

        LNX_LOGGER.debug("sample({}) loaded {} SVs", mCurrentSampleId, mAllVariants.size());
    }

    public void analyse()
    {
        if(mAllVariants.isEmpty())
            return;

        LNX_LOGGER.debug("sample({}) analysing {} variants", mCurrentSampleId, mAllVariants.size());

        mPerfCounters.get(PERF_COUNTER_PREP).start();

        annotateAndFilterVariants();

        mAnalyser.setSampleData(mCurrentSampleId, mAllVariants);

        mAnalyser.preClusteringPreparation();

        if(mConfig.IsGermline)
        {
            buildGermlineCopyNumberData();
        }

        mPerfCounters.get(PERF_COUNTER_PREP).stop();

        mPerfCounters.get(PERF_COUNTER_CLUSTER).start();

        mIsValid = mAnalyser.clusterAndAnalyse();

        mPerfCounters.get(PERF_COUNTER_CLUSTER).stop();
    }

    public void annotate()
    {
        mPerfCounters.get(PERF_COUNTER_ANNOTATE).start();

        mAnalyser.annotateClusters();

        mSvAnnotators.PseudoGeneFinder.checkPseudoGeneAnnotations(mAnalyser.getClusters(), mVisSampleData);

        mPerfCounters.get(PERF_COUNTER_ANNOTATE).stop();
    }

    private void writeOutput()
    {
        // if processing a single sample, write flat-files and optionally load the same data to the DB
        // if running in batch mode, skip flat-file generation and DB load, and instead write verbose batch output files
        if(mConfig.OutputDataPath.isEmpty())
            return;

        mPerfCounters.get(PERF_COUNTER_WRITE).start();

        boolean prepareSampleData = mConfig.isSingleSample();

        List<SvCluster> allClusters = mAnalyser.getAllClusters();

        final List<LinxSvAnnotation> linxSvData = prepareSampleData ? generateSvDataOutput() : null;
        final List<LinxCluster> clusterData = prepareSampleData ? generateClusterOutput(allClusters) : null;
        final List<LinxLink> linksData = prepareSampleData ? generateLinksOutput() : null;

        if(mCohortDataWriter.writeCohortFiles())
        {
            mCohortDataWriter.writeSvData(mCurrentSampleId, mAllVariants);
            mCohortDataWriter.writeLinksData(mCurrentSampleId, mAnalyser.getClusters());
            mCohortDataWriter.writeClusterData(mCurrentSampleId, allClusters);
        }

        mCohortDataWriter.getVisWriter().writeOutput(
                mVisSampleData, mAnalyser.getClusters(), mAllVariants, mCnDataLoader.getChrCnDataMap());

        List<DriverCatalog> driverCatalogs = Lists.newArrayList();
        List<LinxDriver> linxDrivers = Lists.newArrayList();

        if(!mConfig.IsGermline)
        {
            if(mDriverGeneAnnotator != null)
            {
                driverCatalogs.addAll(mDriverGeneAnnotator.getDriverCatalog());
                linxDrivers.addAll(mDriverGeneAnnotator.getDriverOutputList());
            }

            // include any reportable disruptions for genes which do not already have homozgous distruptins
            mFusionAnalyser.getDisruptionFinder().addReportableDisruptions(driverCatalogs);
        }

        if(mConfig.isSingleSample())
        {
            try
            {
                // write per-sample DB-style output
                LinxSvAnnotation.write(LinxSvAnnotation.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), linxSvData);
                LinxCluster.write(LinxCluster.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), clusterData);
                LinxLink.write(LinxLink.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), linksData);

                if(!mConfig.IsGermline)
                {
                    final String driverCatalogFile = LinxDriver.generateCatalogFilename(mConfig.OutputDataPath, mCurrentSampleId, true);
                    DriverCatalogFile.write(driverCatalogFile, driverCatalogs);

                    final String driversFile = LinxDriver.generateFilename(mConfig.OutputDataPath, mCurrentSampleId);
                    LinxDriver.write(driversFile, linxDrivers);
                }
            }
            catch (IOException e)
            {
                LNX_LOGGER.error("failed to write sample SV data: {}", e.toString());
                e.printStackTrace();
                System.exit(1);
            }
        }

        mPerfCounters.get(PERF_COUNTER_WRITE).stop();
    }

    public void writeSampleWithNoSVs()
    {
        try
        {
            LinxSvAnnotation.write(LinxSvAnnotation.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), Collections.EMPTY_LIST);
            LinxCluster.write(LinxCluster.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), Collections.EMPTY_LIST);
            LinxLink.write(LinxLink.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), Collections.EMPTY_LIST);
            LinxBreakend.write(LinxBreakend.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), Collections.EMPTY_LIST);

            // write out the Purple drivers again as would usually be done with SVs
            List<DriverCatalog> purpleDrivers = Lists.newArrayList();

            if(!mConfig.IsGermline) // germline doesn't re-write Purple drivers, and somatic may soon follow suit
            {
                String purpleDriverFile = DriverCatalogFile.generateSomaticFilename(mConfig.PurpleDataPath, mCurrentSampleId);
                purpleDrivers.addAll(DriverCatalogFile.read(purpleDriverFile));
            }

            DriverCatalogFile.write(LinxDriver.generateCatalogFilename(mConfig.OutputDataPath, mCurrentSampleId, !mConfig.IsGermline), purpleDrivers);

            if(mConfig.IsGermline)
            {
                LinxGermlineSv.write(LinxGermlineSv.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Collections.EMPTY_LIST);
            }
            else
            {
                LinxFusion.write(LinxFusion.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Collections.EMPTY_LIST);
                LinxDriver.write(LinxDriver.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Collections.EMPTY_LIST);

                VisSvData.write(VisSvData.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), Collections.EMPTY_LIST);
                VisCopyNumber.write(VisCopyNumber.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), Collections.EMPTY_LIST);
                VisGeneExon.write(VisGeneExon.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), Collections.EMPTY_LIST);
                VisSegment.write(VisSegment.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), Collections.EMPTY_LIST);
                VisFusion.write(VisFusion.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), Collections.EMPTY_LIST);
                VisProteinDomain.write(VisProteinDomain.generateFilename(mConfig.OutputDataPath, mCurrentSampleId, mConfig.IsGermline), Collections.EMPTY_LIST);
            }
        }
        catch (IOException e)
        {
            LNX_LOGGER.error("failed to write sample SV data: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void annotateAndFilterVariants()
    {
        int currentIndex = 0;

        while(currentIndex < mAllVariants.size())
        {
            SvVarData var = mAllVariants.get(currentIndex);

            var.setFragileSites(
                    mSvAnnotators.FragileSiteAnnotator.isFragileSite(var, true),
                    mSvAnnotators.FragileSiteAnnotator.isFragileSite(var, false));

            mSvAnnotators.LineElementAnnotator.setKnownLineElements(var);

            ChromosomeArm startArm = getChromosomalArm(var.chromosome(true), var.position(true));

            ChromosomeArm endArm;
            if(!var.isSglBreakend())
                endArm = getChromosomalArm(var.chromosome(false), var.position(false));
            else
                endArm = ChromosomeArm.P_ARM;

            var.setChromosomalArms(startArm, endArm);

            ++currentIndex;
        }
    }

    private void buildGermlineCopyNumberData()
    {
        CnSegmentBuilder cnSegmentBuilder = new CnSegmentBuilder();
        cnSegmentBuilder.setAllelePloidies(1, 1);

        cnSegmentBuilder.createGermlineCopyNumberData(mCnDataLoader, mAnalyser.getState().getChrBreakendMap());

        cnSegmentBuilder.setSamplePurity(mCnDataLoader, 1, 2, MALE);

        mCnDataLoader.createChrCopyNumberMap();

        SvCNData.setSvCopyNumberData(
                mAllVariants,
                mCnDataLoader.getSvJcnCalcMap(),
                mCnDataLoader.getSvIdCnDataMap(),
                mCnDataLoader.getChrCnDataMap());
    }

    private List<LinxSvAnnotation> generateSvDataOutput()
    {
        final List<LinxSvAnnotation> linxSvData = Lists.newArrayList();

        for(final SvVarData var : mAllVariants)
        {
            final SvCluster cluster = var.getCluster();

            if(cluster == null)
            {
                LNX_LOGGER.error("SV({}) not assigned to a cluster", var.posId());
                continue;
            }

            final ArmCluster armClusterStart = cluster.findArmCluster(var.getBreakend(true));

            if(armClusterStart == null)
            {
                LNX_LOGGER.error("start breakend({}) not assigned to arm cluster", var.getBreakend(true));
                continue;
            }

            final ArmCluster armClusterEnd = !var.isSglBreakend() ? cluster.findArmCluster(var.getBreakend(false)) : null;

            // suppress SGL-INF artificial SVs, but use their details in place of the SGL
            if(var.getLinkedSVs() != null && var.isSglBreakend())
                continue;

            linxSvData.add(ImmutableLinxSvAnnotation.builder()
                    .vcfId(var.getSvData().vcfId())
                    .svId(var.id())
                    .clusterId(cluster.id())
                    .clusterReason(var.getClusterReason())
                    .fragileSiteStart(var.isFragileSite(true))
                    .fragileSiteEnd(var.isFragileSite(false))
                    .isFoldback(var.isFoldback())
                    .lineTypeStart(LineElementType.toString(var.getLineElement(true)))
                    .lineTypeEnd(LineElementType.toString(var.getLineElement(false)))
                    .junctionCopyNumberMin(var.jcnMin())
                    .junctionCopyNumberMax(var.jcnMax())
                    .geneStart(var.getGeneInBreakend(true, false))
                    .geneEnd(var.getGeneInBreakend(false, false))
                    .localTopologyIdStart(armClusterStart.id())
                    .localTopologyIdEnd(armClusterEnd != null ? armClusterEnd.id() : -1)
                    .localTopologyStart(armClusterStart.getTypeStr())
                    .localTopologyEnd(armClusterEnd != null ? armClusterEnd.getTypeStr() : "")
                    .localTICountStart(armClusterStart.getTICount())
                    .localTICountEnd(armClusterEnd != null ? armClusterEnd.getTICount() : 0)
                    .build());
        }

        return linxSvData;
    }

    private List<LinxCluster> generateClusterOutput(final List<SvCluster> allClusters)
    {
        final List<LinxCluster> clusterData = Lists.newArrayList();

        for(final SvCluster cluster : allClusters)
        {
            final String superType = getClusterCategory(cluster);

                clusterData.add(ImmutableLinxCluster.builder()
                        .clusterId(cluster.id())
                        .category(superType)
                        .synthetic(cluster.isSyntheticType())
                        .resolvedType(cluster.getResolvedType().toString())
                        .clusterCount(cluster.getSvCount())
                        .clusterDesc(cluster.getDesc())
                        .build());
        }

        return clusterData;
    }

    private List<LinxLink> generateLinksOutput()
    {
        final List<LinxLink> linksData = Lists.newArrayList();

        for(final SvCluster cluster : mAnalyser.getClusters())
        {
            List<SvChain> chains = cluster.getChains();

            for(final SvChain chain : chains)
            {
                int chainSvCount = chain.getSvCount();

                List<LinkedPair> uniquePairs = Lists.newArrayList();
                final List<LinkedPair> chainLinks = chain.getLinkedPairs();
                boolean isDoubleMinute = chain.isDoubleMinute();

                for(int chainIndex = 0; chainIndex < chainLinks.size(); ++chainIndex)
                {
                    final LinkedPair pair = chainLinks.get(chainIndex);

                    if(uniquePairs.stream().anyMatch(x -> x.matches(pair)))
                        continue;

                    uniquePairs.add(pair);

                    String chainIndexStr = String.valueOf(chainIndex);

                    for(int j = chainIndex + 1; j < chainLinks.size(); ++j)
                    {
                        if(chainLinks.get(j).matches(pair))
                        {
                            chainIndexStr = appendStr(chainIndexStr, String.valueOf(j), ';');
                        }
                    }

                    final SvBreakend beStart = pair.getBreakend(true);
                    final SvBreakend beEnd = pair.getBreakend(false);

                    final String linkArm = beStart.arm() == beEnd.arm() ? asStr(beStart.arm()) : asStr(beStart.arm()) + "_" + asStr(beEnd.arm());

                    linksData.add(ImmutableLinxLink.builder()
                            .clusterId(cluster.id())
                            .chainId(chain.id())
                            .chainCount(chainSvCount)
                            .chainIndex(chainIndexStr)
                            .lowerSvId(beStart.getSV().id())
                            .upperSvId(beEnd.getSV().id())
                            .lowerBreakendIsStart(beStart.usesStart())
                            .upperBreakendIsStart(beEnd.usesStart())
                            .chromosome(beStart.chromosome())
                            .arm(linkArm)
                            .assembled(pair.isAssembled())
                            .traversedSVCount(pair.getTraversedSVCount())
                            .length(pair.baseLength())
                            .junctionCopyNumber(Doubles.round(chain.jcn(), 4))
                            .junctionCopyNumberUncertainty(Doubles.round(chain.jcnUncertainty(), 4))
                            .pseudogeneInfo(pair.getExonMatchData())
                            .ecDna(isDoubleMinute)
                            .build());
                }
            }
        }

        return linksData;
    }

    public void close()
    {
        mAnalyser.close();
        mFusionAnalyser.close();

        /* TODO
        if(mConfig.hasMultipleSamples() || LNX_LOGGER.isDebugEnabled())
        {
            // log perf stats
            mPcPrep.logStats();
            mPcClusterAnalyse.logStats();
            mAnalyser.logStats();
            mPcAnnotation.logStats();
            mPcWrite.logStats();
        }
        */
    }
}
