package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.purple.Gender.MALE;
import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.SvFileLoader.createSvData;
import static com.hartwig.hmftools.linx.SvFileLoader.loadSampleSvDataFromFile;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.getClusterCategory;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.linkSglMappedInferreds;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.linx.drivers.DriverGeneAnnotator.LINX_DRIVER_CATALOG;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.asStr;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxCluster;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxLink;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxCluster;
import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxLink;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.sv.linx.LinxViralInsertion;
import com.hartwig.hmftools.linx.analysis.ClusterAnalyser;
import com.hartwig.hmftools.linx.analysis.VariantPrep;
import com.hartwig.hmftools.linx.annotators.LineElementType;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.CnSegmentBuilder;
import com.hartwig.hmftools.linx.drivers.DriverGeneAnnotator;
import com.hartwig.hmftools.linx.fusion.FusionDisruptionAnalyser;
import com.hartwig.hmftools.linx.fusion.FusionResources;
import com.hartwig.hmftools.linx.types.ArmCluster;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumberFile;
import com.hartwig.hmftools.linx.visualiser.file.VisFusionFile;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExonFile;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomainFile;
import com.hartwig.hmftools.linx.visualiser.file.VisSampleData;
import com.hartwig.hmftools.linx.visualiser.file.VisSegmentFile;
import com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.DatabaseUtil;

public class SampleAnalyser implements Callable
{
    private final int mId;
    private final LinxConfig mConfig;
    private List<String> mSampleIds;

    private final ClusterAnalyser mAnalyser;
    private final DatabaseAccess mDbAccess;
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

    private final PerformanceCounter mPcTotal;
    private final PerformanceCounter mPcPrep;
    private final PerformanceCounter mPcClusterAnalyse;
    private final PerformanceCounter mPcAnnotation;
    private final PerformanceCounter mPcWrite;

    public SampleAnalyser(
            int instanceId, final LinxConfig config, final DatabaseAccess dbAccess, final SvAnnotators svAnnotators,
            final EnsemblDataCache ensemblDataCache, final FusionResources fusionResources, final CohortDataWriter cohortDataWriter)
    {
        mId = instanceId;
        mConfig = config;
        mSampleIds = Lists.newArrayList();

        mVisSampleData = new VisSampleData();
        mCurrentSampleId = "";

        mDbAccess = dbAccess;
        mEnsemblDataCache = ensemblDataCache;

        mCohortDataWriter = cohortDataWriter;
        mSvAnnotators = svAnnotators;

        mAnalyser = new ClusterAnalyser(config, mCohortDataWriter);

        mCnDataLoader = new CnDataLoader(config.PurpleDataPath, dbAccess);

        mDriverGeneAnnotator = mConfig.RunDrivers ?
                new DriverGeneAnnotator(dbAccess, ensemblDataCache, config, mCnDataLoader, cohortDataWriter, mVisSampleData) : null;

        mFusionAnalyser = new FusionDisruptionAnalyser(
                config.CmdLineArgs, config, ensemblDataCache, dbAccess, fusionResources, cohortDataWriter, mVisSampleData);

        mAllVariants = Lists.newArrayList();

        mIsValid = true;

        mAnalyser.setLineAnnotator(mSvAnnotators.LineElementAnnotator);
        mAnalyser.setCnDataLoader(mCnDataLoader);
        mAnalyser.setGeneCollection(ensemblDataCache);

        mPcTotal = new PerformanceCounter("Total");
        mPcPrep = new PerformanceCounter("Preparation");
        mPcClusterAnalyse = new PerformanceCounter("ClusterAndAnalyse");
        mPcAnnotation = new PerformanceCounter("Annotation");
        mPcWrite = new PerformanceCounter("WriteCSV");
    }

    public void setSampleIds(final List<String> sampleIds)
    {
        mSampleIds.clear();
        mSampleIds.addAll(sampleIds);
    }

    @Override
    public Long call()
    {
        // processSample();
        return (long)1;
    }

    public void processSamples()
    {
        if(mSampleIds.size() == 1)
        {
            mPcTotal.start();
            processSample(mSampleIds.get(0));
            mPcTotal.stop();
            return;
        }

        LNX_LOGGER.info("%d: processing {} samples", mId, mSampleIds.size());

        for(int i = 0; i < mSampleIds.size(); ++i)
        {
            mPcTotal.start();

            processSample(mSampleIds.get(i));

            if(i > 10 && (i % 10) == 0)
            {
                LNX_LOGGER.info("%d: procesed {} samples", mId, i);
            }

            mPcTotal.stop();

            if(!inValidState())
                break;
        }
    }

    private void processSample(final String sampleId)
    {
        mCurrentSampleId = sampleId;
        mVisSampleData.setSampleId(sampleId);

        final List<StructuralVariantData> svRecords = mConfig.loadSampleDataFromFile() ?
                loadSampleSvDataFromFile(mConfig, mCurrentSampleId, mConfig.CmdLineArgs) : mDbAccess.readStructuralVariantData(mCurrentSampleId);

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
            VariantPrep.setSvGeneData(svDataList, mEnsemblDataCache, mConfig.RunFusions, mConfig.breakendGeneLoading());
        }

        analyse();

        if(!inValidState())
        {
            LNX_LOGGER.info("exiting after sample({}), in invalid state", mCurrentSampleId);
            return;
        }

        // when matching RNA, allow all transcripts regardless of their viability for fusions
        boolean purgeInvalidTranscripts = !mFusionAnalyser.hasRnaSampleData();
        mFusionAnalyser.annotateTranscripts(svDataList, purgeInvalidTranscripts);

        annotate();

        if(mDriverGeneAnnotator != null)
            mDriverGeneAnnotator.annotateSVs(mCurrentSampleId, getChrBreakendMap());

        if(mConfig.RunFusions || mConfig.IsGermline)
            mFusionAnalyser.run(mCurrentSampleId, svDataList, getClusters(), getChrBreakendMap());

        writeOutput();
        close();

        LNX_LOGGER.info("sample({}) procesed {} SVs", sampleId, svDataList.size());
    }

    public final List<SvVarData> getVariants() { return mAllVariants; }
    public final List<SvCluster> getClusters() { return mAnalyser.getClusters(); }
    public boolean inValidState() { return mIsValid; }
    public final Map<String, List<SvBreakend>> getChrBreakendMap() { return mAnalyser.getState().getChrBreakendMap(); }

    public void setSampleSVs(final List<SvVarData> variants)
    {
        if (variants.isEmpty())
            return;

        mAllVariants.addAll(variants);

        if(!mConfig.IsGermline)
        {
            // look-up and cache relevant CN data into each SV
            VariantPrep.setSvCopyNumberData(
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

        mPcPrep.start();

        linkSglMappedInferreds(mAllVariants);
        annotateAndFilterVariants();

        mAnalyser.setSampleData(mCurrentSampleId, mAllVariants);

        mAnalyser.preClusteringPreparation();

        if(mConfig.IsGermline)
        {
            buildGermlineCopyNumberData();
        }

        mSvAnnotators.KataegisAnnotator.annotateVariants(mCurrentSampleId, mAnalyser.getState().getChrBreakendMap());
        mSvAnnotators.ReplicationOriginAnnotator.setReplicationOrigins(mAnalyser.getState().getChrBreakendMap());

        mPcPrep.stop();

        mPcClusterAnalyse.start();

        mIsValid = mAnalyser.clusterAndAnalyse();

        mPcClusterAnalyse.stop();
    }

    public void annotate()
    {
        mPcAnnotation.start();

        mAnalyser.annotateClusters();

        mSvAnnotators.PseudoGeneFinder.checkPseudoGeneAnnotations(getClusters(), mVisSampleData);

        if(mSvAnnotators.IndelAnnotator != null)
        {
            mSvAnnotators.IndelAnnotator.loadIndels(mCurrentSampleId);

            if(!mSvAnnotators.IndelAnnotator.exceedsThresholds())
                getClusters().forEach(x -> mSvAnnotators.IndelAnnotator.annotateCluster(x));
        }

        mPcAnnotation.stop();
    }

    private void writeOutput()
    {
        // if processing a single sample, write flat-files and optionally load the same data to the DB
        // if running in batch mode, skip flat-file generation and DB load, and instead write verbose batch output files
        if(mConfig.OutputDataPath.isEmpty())
            return;

        mPcWrite.start();

        boolean prepareSampleData = mConfig.isSingleSample() || mConfig.UploadToDB;

        final List<LinxSvAnnotation> linxSvData = prepareSampleData ? generateSvDataOutput() : null;
        final List<LinxCluster> clusterData = prepareSampleData ? generateClusterOutput() : null;
        final List<LinxLink> linksData = prepareSampleData ? generateLinksOutput() : null;

        if(mCohortDataWriter.writeCohortFiles())
        {
            mCohortDataWriter.writeSvData(mCurrentSampleId, mAllVariants);
            mCohortDataWriter.writeLinksData(mCurrentSampleId, mAnalyser.getClusters());
            mCohortDataWriter.writeClusterData(mCurrentSampleId, mAnalyser.getClusters());
        }

        mCohortDataWriter.getVisWriter().writeOutput(
                mVisSampleData, mAnalyser.getClusters(), mAllVariants, mCnDataLoader.getChrCnDataMap());

        if(mConfig.isSingleSample())
        {
            try
            {
                // write per-sample DB-style output
                LinxSvAnnotation.write(LinxSvAnnotation.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), linxSvData);
                LinxCluster.write(LinxCluster.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), clusterData);
                LinxLink.write(LinxLink.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), linksData);
            }
            catch (IOException e)
            {
                LNX_LOGGER.error("failed to write sample SV data: {}", e.toString());
            }
        }

        if(mConfig.UploadToDB && mDbAccess != null)
        {
            writeToDatabase(mCurrentSampleId, mDbAccess, linxSvData, clusterData, linksData);
        }

        mPcWrite.stop();
    }

    private synchronized static void writeToDatabase(
            final String sampleId, final DatabaseAccess dbAccess,
            final List<LinxSvAnnotation> linxSvData, final List<LinxCluster> clusterData, final List<LinxLink> linksData)
    {
        dbAccess.writeSvLinxData(sampleId, linxSvData);
        dbAccess.writeSvClusters(sampleId, clusterData);
        dbAccess.writeSvLinks(sampleId, linksData);
    }

    public void writeSampleWithNoSVs()
    {
        try
        {
            LinxSvAnnotation.write(LinxSvAnnotation.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());
            LinxCluster.write(LinxCluster.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());
            LinxLink.write(LinxLink.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());

            LinxFusion.write(LinxFusion.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());
            LinxBreakend.write(LinxBreakend.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());
            LinxDriver.write(LinxDriver.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());

            final String driverCatalogFile = mConfig.OutputDataPath + mCurrentSampleId + LINX_DRIVER_CATALOG;
            DriverCatalogFile.write(driverCatalogFile, Lists.newArrayList());

            VisSvDataFile.write(VisSvDataFile.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());
            VisCopyNumberFile.write(VisCopyNumberFile.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());
            VisGeneExonFile.write(VisGeneExonFile.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());
            VisSegmentFile.write(VisSegmentFile.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());
            VisFusionFile.write(VisFusionFile.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());
            VisProteinDomainFile.write(VisProteinDomainFile.generateFilename(mConfig.OutputDataPath, mCurrentSampleId), Lists.newArrayList());
        }
        catch (IOException e)
        {
            LNX_LOGGER.error("failed to write sample SV data: {}", e.toString());
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

        cnSegmentBuilder.createIndependentCopyNumberData(mCnDataLoader, mAnalyser.getState().getChrBreakendMap());

        cnSegmentBuilder.setSamplePurity(mCnDataLoader, 1, 2, MALE);

        mCnDataLoader.createChrCopyNumberMap();

        VariantPrep.setSvCopyNumberData(
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
                    .replicationTimingStart(var.getReplicationOrigin(true))
                    .replicationTimingEnd(var.getReplicationOrigin(false))
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

    private List<LinxCluster> generateClusterOutput()
    {
        final List<LinxCluster> clusterData = Lists.newArrayList();

        for(final SvCluster cluster : getClusters())
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

        for(final SvCluster cluster : getClusters())
        {
            List<SvChain> chains = cluster.getChains();

            for (final SvChain chain : chains)
            {
                int chainSvCount = chain.getSvCount();

                List<LinkedPair> uniquePairs = Lists.newArrayList();
                final List<LinkedPair> chainLinks = chain.getLinkedPairs();
                boolean isDoubleMinute = chain.isDoubleMinute();

                for (int chainIndex = 0; chainIndex < chainLinks.size(); ++chainIndex)
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
                            .length(pair.length())
                            .junctionCopyNumber(DatabaseUtil.decimal(chain.jcn()))
                            .junctionCopyNumberUncertainty(DatabaseUtil.decimal(chain.jcnUncertainty()))
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

        mAnalyser.close();
    }
}
