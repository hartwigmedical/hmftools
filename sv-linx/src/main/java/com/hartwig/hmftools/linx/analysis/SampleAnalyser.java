package com.hartwig.hmftools.linx.analysis;

import static com.hartwig.hmftools.common.purple.gender.Gender.MALE;
import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.getClusterCategory;
import static com.hartwig.hmftools.linx.analysis.ClusteringPrep.linkSglMappedInferreds;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.linx.annotators.ViralInsertAnnotator.VH_ID;
import static com.hartwig.hmftools.linx.annotators.ViralInsertAnnotator.VH_NAME;
import static com.hartwig.hmftools.linx.drivers.DriverGeneAnnotator.LINX_DRIVER_CATALOG;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.types.ChromosomeArm.asStr;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxLink;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxLink;
import com.hartwig.hmftools.common.variant.structural.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.linx.annotators.IndelAnnotator;
import com.hartwig.hmftools.linx.annotators.KataegisAnnotator;
import com.hartwig.hmftools.linx.annotators.LineElementAnnotator;
import com.hartwig.hmftools.linx.annotators.LineElementType;
import com.hartwig.hmftools.linx.annotators.PseudoGeneFinder;
import com.hartwig.hmftools.linx.annotators.ReplicationOriginAnnotator;
import com.hartwig.hmftools.linx.annotators.ViralInsertAnnotator;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.CnSegmentBuilder;
import com.hartwig.hmftools.linx.cn.JcnCalcData;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.types.ArmCluster;
import com.hartwig.hmftools.linx.types.ChromosomeArm;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumberFile;
import com.hartwig.hmftools.linx.visualiser.file.VisFusionFile;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExonFile;
import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomainFile;
import com.hartwig.hmftools.linx.visualiser.file.VisSegmentFile;
import com.hartwig.hmftools.linx.visualiser.file.VisSvDataFile;
import com.hartwig.hmftools.linx.visualiser.file.VisualiserWriter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.DatabaseUtil;

public class SampleAnalyser
{
    private final LinxConfig mConfig;
    private final ClusterAnalyser mAnalyser;
    private final CohortDataWriter mCohortDataWriter;

    // data per run (ie sample)
    private String mSampleId;
    private final List<SvVarData> mAllVariants; // full list of SVs

    private final FragileSiteAnnotator mFragileSiteAnnotator;
    private final LineElementAnnotator mLineElementAnnotator;
    private final ReplicationOriginAnnotator mReplicationOriginAnnotator;
    private final ViralInsertAnnotator mViralInsertAnnotator;
    private final KataegisAnnotator mKataegisAnnotator;
    private CnDataLoader mCnDataLoader;
    private final PseudoGeneFinder mPseudoGeneFinder;
    private final IndelAnnotator mIndelAnnotator;

    private boolean mIsValid;

    private final PerformanceCounter mPcPrep;
    private final PerformanceCounter mPcClusterAnalyse;
    private final PerformanceCounter mPcAnnotation;
    private final PerformanceCounter mPcWrite;

    public SampleAnalyser(final LinxConfig config, DatabaseAccess dbAccess)
    {
        mConfig = config;
        mSampleId = "";

        mAnalyser = new ClusterAnalyser(config);

        mAllVariants = Lists.newArrayList();
        mCnDataLoader = null;

        mIsValid = true;

        mCohortDataWriter = new CohortDataWriter(config, mAnalyser);

        mFragileSiteAnnotator = new FragileSiteAnnotator();
        mFragileSiteAnnotator.loadFragileSitesFile(mConfig.FragileSiteFile);

        mLineElementAnnotator = new LineElementAnnotator(mConfig.ProximityDistance);
        mLineElementAnnotator.loadLineElementsFile(mConfig.LineElementFile);
        mAnalyser.setLineAnnotator(mLineElementAnnotator);

        mPseudoGeneFinder = new PseudoGeneFinder(mCohortDataWriter.getVisWriter());
        mLineElementAnnotator.setPseudoGeneFinder(mPseudoGeneFinder);

        mReplicationOriginAnnotator = new ReplicationOriginAnnotator();
        mReplicationOriginAnnotator.loadReplicationOrigins(mConfig.ReplicationOriginsFile);

        mViralInsertAnnotator = new ViralInsertAnnotator();
        mViralInsertAnnotator.loadViralHostData(mConfig.ViralHostsFile);

        mKataegisAnnotator = new KataegisAnnotator(mConfig.OutputDataPath);
        mKataegisAnnotator.loadKataegisData(mConfig.KataegisFile);

        mIndelAnnotator = config.IndelAnnotation ? new IndelAnnotator(dbAccess, config) : null;

        mPcPrep = new PerformanceCounter("Preparation");
        mPcClusterAnalyse = new PerformanceCounter("ClusterAndAnalyse");
        mPcAnnotation = new PerformanceCounter("Annotation");
        mPcWrite = new PerformanceCounter("WriteCSV");
    }

    public final List<SvVarData> getVariants() { return mAllVariants; }
    public final List<SvCluster> getClusters() { return mAnalyser.getClusters(); }
    public boolean inValidState() { return mIsValid; }
    public final Map<String, List<SvBreakend>> getChrBreakendMap() { return mAnalyser.getState().getChrBreakendMap(); }
    public final VisualiserWriter getVisWriter() { return mCohortDataWriter.getVisWriter(); }

    public void setCnDataLoader(CnDataLoader cnAnalyser)
    {
        mCnDataLoader = cnAnalyser;
        mAnalyser.setCnDataLoader(cnAnalyser);
    }

    public void setGeneCollection(final EnsemblDataCache geneDataCache)
    {
        mAnalyser.setGeneCollection(geneDataCache);
        mPseudoGeneFinder.setGeneTransCache(geneDataCache);
    }

    private void clearState()
    {
        if(mSampleId.isEmpty())
            return;

        mSampleId = "";
        mAllVariants.clear();
    }

    public void setSampleId(final String sampleId)
    {
        clearState();

        mSampleId = sampleId;
        mCohortDataWriter.getVisWriter().setSampleId(sampleId);
    }

    public void setSampleSVs(final List<SvVarData> variants)
    {
        if (variants.isEmpty())
            return;

        mAllVariants.addAll(variants);

        if(!mConfig.IsGermline)
        {
            // look-up and cache relevant CN data into each SV
            setSvCopyNumberData(
                    mAllVariants,
                    mCnDataLoader.getSvJcnCalcMap(),
                    mCnDataLoader.getSvIdCnDataMap(),
                    mCnDataLoader.getChrCnDataMap());
        }

        LNX_LOGGER.debug("loaded {} SVs", mAllVariants.size());
    }

    public static void setSvCopyNumberData(List<SvVarData> svList, final Map<Integer, JcnCalcData> svJcnCalcDataMap,
            final Map<Integer,SvCNData[]> svIdCnDataMap, final Map<String,List<SvCNData>> chrCnDataMap)
    {
        if((svJcnCalcDataMap == null || svJcnCalcDataMap.isEmpty()) && svIdCnDataMap.isEmpty())
            return;

        List<SvCNData> cnDataList = null;
        String currentChromosome = "";
        for(final SvVarData var : svList)
        {
            if(svJcnCalcDataMap != null)
            {
                final JcnCalcData jcnData = svJcnCalcDataMap.get(var.id());
                if (jcnData != null)
                {
                    double estJcn = jcnData.JcnEstimate;
                    double estUncertainty = jcnData.JcnUncertainty;
                    var.setJcnRecalcData(estJcn - estUncertainty, estJcn + estUncertainty);
                }
            }

            final SvCNData[] cnDataPair = svIdCnDataMap.get(var.id());

            if(cnDataPair == null)
                continue;

            for(int be = SE_START; be <= SE_END; ++be)
            {
                if(var.isSglBreakend() && be == SE_END)
                    continue;

                boolean isStart = isStart(be);

                if(!currentChromosome.equals(var.chromosome(isStart)))
                {
                    currentChromosome = var.chromosome(isStart);
                    cnDataList = chrCnDataMap.get(currentChromosome);
                }

                SvCNData cnDataPost = cnDataPair[be];

                if(cnDataList == null || cnDataPost == null)
                    continue;

                SvCNData cnDataPrev = cnDataList.get(cnDataPost.getIndex() - 1);

                var.setCopyNumberData(isStart, cnDataPrev, cnDataPost);
            }
        }
    }

    public static void setSvGeneData(
            final List<SvVarData> svList, final EnsemblDataCache ensemblDataCache, boolean applyPromotorDistance, boolean loadBreakendGenes)
    {
        int upstreamDistance = applyPromotorDistance ? PRE_GENE_PROMOTOR_DISTANCE : 0;

        if (loadBreakendGenes)
        {
            // only load transcript info for the genes covered
            final List<String> restrictedGeneIds = Lists.newArrayList();

            for (final SvVarData var : svList)
            {
                for (int be = SE_START; be <= SE_END; ++be)
                {
                    if (be == SE_END && var.isSglBreakend())
                    {
                        // special case of looking for mappings to locations containing genes so hotspot fusions can be found
                        for(final SglMapping mapping : var.getSglMappings())
                        {
                            ensemblDataCache.populateGeneIdList(restrictedGeneIds, mapping.Chromosome, mapping.Position, upstreamDistance);
                        }
                    }
                    else
                    {
                        boolean isStart = isStart(be);
                        ensemblDataCache.populateGeneIdList(restrictedGeneIds, var.chromosome(isStart), var.position(isStart), upstreamDistance);
                    }
                }
            }

            ensemblDataCache.getAlternativeGeneData().stream().filter(x -> !restrictedGeneIds.contains(x.GeneId))
                    .forEach(x -> restrictedGeneIds.add(x.GeneId));

            ensemblDataCache.loadTranscriptData(restrictedGeneIds);
        }

        // associate breakends with transcripts
        for (final SvVarData var : svList)
        {
            for (int be = SE_START; be <= SE_END; ++be)
            {
                boolean isStart = isStart(be);
                final List<BreakendGeneData> genesList = var.getGenesList(isStart);

                if (be == SE_END && var.isSglBreakend())
                {
                    // special case of looking for mappings to locations containing genes so hotspot fusions can be found
                    for(final SglMapping mapping : var.getSglMappings())
                    {
                        final List<BreakendGeneData> mappingGenes = ensemblDataCache.findGeneAnnotationsBySv(
                                var.id(), isStart, mapping.Chromosome, mapping.Position, mapping.Orientation, upstreamDistance);

                        mappingGenes.forEach(x -> x.setType(var.type()));

                        genesList.addAll(mappingGenes);
                    }
                }
                else
                {
                    genesList.addAll(ensemblDataCache.findGeneAnnotationsBySv(
                            var.id(), isStart, var.chromosome(isStart), var.position(isStart), var.orientation(isStart), upstreamDistance));

                    for (BreakendGeneData gene : genesList)
                    {
                        gene.setSvData(var.getSvData(), var.jcn());
                    }
                }
            }
        }
    }

    public void analyse()
    {
        if(mAllVariants.isEmpty())
            return;

        LNX_LOGGER.debug("sample({}) analysing {} variants", mSampleId, mAllVariants.size());

        mPcPrep.start();

        linkSglMappedInferreds(mAllVariants);
        annotateAndFilterVariants();

        mAnalyser.setSampleData(mSampleId, mAllVariants);

        mAnalyser.preClusteringPreparation();

        if(mConfig.IsGermline)
        {
            buildGermlineCopyNumberData();
        }

        mKataegisAnnotator.annotateVariants(mSampleId, mAnalyser.getState().getChrBreakendMap());
        mReplicationOriginAnnotator.setReplicationOrigins(mAnalyser.getState().getChrBreakendMap());

        mPcPrep.stop();

        mPcClusterAnalyse.start();

        mIsValid = mAnalyser.clusterAndAnalyse();

        mPcClusterAnalyse.stop();
    }

    public void annotate()
    {
        mPcAnnotation.start();

        mAnalyser.annotateClusters();

        mPseudoGeneFinder.checkPseudoGeneAnnotations(getClusters());

        if(mIndelAnnotator != null)
        {
            mIndelAnnotator.loadIndels(mSampleId);

            if(!mIndelAnnotator.exceedsThresholds())
                getClusters().forEach(x -> mIndelAnnotator.annotateCluster(x));
        }

        mPcAnnotation.stop();
    }

    public void writeOutput(final DatabaseAccess dbAccess)
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
        final List<LinxViralInsertion> viralInserts = prepareSampleData ? generateViralInserts() : null;

        if(mCohortDataWriter.writeCohortFiles())
        {
            mCohortDataWriter.writeSvData(mSampleId, mAllVariants, mViralInsertAnnotator);
            mCohortDataWriter.writeLinksData(mSampleId, mAnalyser.getClusters());
            mCohortDataWriter.writeClusterData(mSampleId, mAnalyser.getClusters());
        }

        mCohortDataWriter.getVisWriter().writeOutput(mAnalyser.getClusters(), mAllVariants, mCnDataLoader.getChrCnDataMap());

        if(mConfig.isSingleSample())
        {
            try
            {
                // write per-sample DB-style output
                LinxSvAnnotation.write(LinxSvAnnotation.generateFilename(mConfig.OutputDataPath, mSampleId), linxSvData);
                LinxCluster.write(LinxCluster.generateFilename(mConfig.OutputDataPath, mSampleId), clusterData);
                LinxLink.write(LinxLink.generateFilename(mConfig.OutputDataPath, mSampleId), linksData);
                LinxViralInsertion.write(LinxViralInsertion.generateFilename(mConfig.OutputDataPath, mSampleId), viralInserts);

            }
            catch (IOException e)
            {
                LNX_LOGGER.error("failed to write sample SV data: {}", e.toString());
            }
        }

        if(mConfig.UploadToDB && dbAccess != null)
        {
            dbAccess.writeSvLinxData(mSampleId, linxSvData);
            dbAccess.writeSvClusters(mSampleId, clusterData);
            dbAccess.writeSvLinks(mSampleId, linksData);
            dbAccess.writeSvViralInserts(mSampleId, viralInserts);
        }

        mPcWrite.stop();
    }

    public void writeSampleWithNoSVs()
    {
        try
        {
            LinxSvAnnotation.write(LinxSvAnnotation.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());
            LinxCluster.write(LinxCluster.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());
            LinxLink.write(LinxLink.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());
            LinxViralInsertion.write(LinxViralInsertion.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());

            LinxFusion.write(LinxFusion.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());
            LinxBreakend.write(LinxBreakend.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());
            LinxDriver.write(LinxDriver.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());

            final String driverCatalogFile = mConfig.OutputDataPath + mSampleId + LINX_DRIVER_CATALOG;
            DriverCatalogFile.write(driverCatalogFile, Lists.newArrayList());

            VisSvDataFile.write(VisSvDataFile.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());
            VisCopyNumberFile.write(VisCopyNumberFile.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());
            VisGeneExonFile.write(VisGeneExonFile.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());
            VisSegmentFile.write(VisSegmentFile.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());
            VisFusionFile.write(VisFusionFile.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());
            VisProteinDomainFile.write(VisProteinDomainFile.generateFilename(mConfig.OutputDataPath, mSampleId), Lists.newArrayList());
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
                    mFragileSiteAnnotator.isFragileSite(var, true),
                    mFragileSiteAnnotator.isFragileSite(var, false));

            mLineElementAnnotator.setKnownLineElements(var);

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

        setSvCopyNumberData(
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

    private List<LinxViralInsertion> generateViralInserts()
    {
        List<LinxViralInsertion> viralInserts = Lists.newArrayList();

        for(final SvVarData var : mAllVariants)
        {
            final String[] viralInsertData = mViralInsertAnnotator.matchesViralInsert(var);

            if(viralInsertData != null)
            {
                viralInserts.add(new LinxViralInsertion(mSampleId, var.id(), viralInsertData[VH_ID], viralInsertData[VH_NAME]));
            }
        }

        return viralInserts;
    }

    public void close()
    {
        if(mConfig.hasMultipleSamples() || LNX_LOGGER.isDebugEnabled())
        {
            // log perf stats
            mPcPrep.logStats();
            mPcClusterAnalyse.logStats();
            mAnalyser.logStats();
            mPcAnnotation.logStats();
            mPcWrite.logStats();
        }

        mCohortDataWriter.close();

        mKataegisAnnotator.close();
        mAnalyser.close();

        if(mIndelAnnotator != null)
            mIndelAnnotator.close();
    }
}
