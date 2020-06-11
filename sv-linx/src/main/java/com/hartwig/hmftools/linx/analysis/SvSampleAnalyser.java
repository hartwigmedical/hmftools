package com.hartwig.hmftools.linx.analysis;

import static com.hartwig.hmftools.common.purple.gender.Gender.MALE;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.SvClassification.getSuperType;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.linx.annotators.ViralInsertAnnotator.VH_ID;
import static com.hartwig.hmftools.linx.annotators.ViralInsertAnnotator.VH_NAME;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;
import static com.hartwig.hmftools.linx.types.ChromosomeArm.asStr;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_COMPLEX_LINE;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_COMPLEX_OTHER;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_FOLDBACK;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_FOLDBACK_DSB;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_ISOLATED_BE;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_SAME_ORIENT;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_SIMPLE_DUP;
import static com.hartwig.hmftools.linx.types.SvArmCluster.ARM_CL_TI_ONLY;
import static com.hartwig.hmftools.linx.types.SvArmCluster.getArmClusterData;
import static com.hartwig.hmftools.linx.types.SvConstants.NO_DB_MARKER;
import static com.hartwig.hmftools.linx.types.SvConstants.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.ReportableDisruptionFile;
import com.hartwig.hmftools.common.fusion.ReportableGeneFusionFile;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxLink;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxSvData;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakendFile;
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.LinxClusterFile;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriverFile;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusionFile;
import com.hartwig.hmftools.common.variant.structural.linx.LinxLink;
import com.hartwig.hmftools.common.variant.structural.linx.LinxLinkFile;
import com.hartwig.hmftools.common.variant.structural.linx.LinxSvData;
import com.hartwig.hmftools.common.variant.structural.linx.LinxSvDataFile;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertFile;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.linx.annotators.IndelAnnotator;
import com.hartwig.hmftools.linx.annotators.KataegisAnnotator;
import com.hartwig.hmftools.linx.annotators.LineElementAnnotator;
import com.hartwig.hmftools.linx.annotators.PseudoGeneFinder;
import com.hartwig.hmftools.linx.annotators.ReplicationOriginAnnotator;
import com.hartwig.hmftools.linx.annotators.ViralInsertAnnotator;
import com.hartwig.hmftools.linx.chaining.ChainMetrics;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.CnSegmentBuilder;
import com.hartwig.hmftools.linx.cn.JcnCalcData;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.types.ChromosomeArm;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvArmCluster;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisualiserWriter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.DatabaseUtil;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public class SvSampleAnalyser {

    private final LinxConfig mConfig;
    private final ClusterAnalyser mAnalyser;

    // data per run (ie sample)
    private String mSampleId;
    private final List<SvVarData> mAllVariants; // full list of SVs

    private BufferedWriter mSvFileWriter;
    private BufferedWriter mClusterFileWriter;
    private BufferedWriter mLinksFileWriter;
    private final VisualiserWriter mVisWriter;

    private final FragileSiteAnnotator mFragileSiteAnnotator;
    private final LineElementAnnotator mLineElementAnnotator;
    private final ReplicationOriginAnnotator mReplicationOriginAnnotator;
    private final ViralInsertAnnotator mViralInsertAnnotator;
    private final KataegisAnnotator mKataegisAnnotator;
    private CnDataLoader mCnDataLoader;
    private final PseudoGeneFinder mPseudoGeneFinder;
    private IndelAnnotator mIndelAnnotator;

    private boolean mIsValid;

    private final PerformanceCounter mPcPrep;
    private final PerformanceCounter mPcClusterAnalyse;
    private final PerformanceCounter mPcAnnotation;
    private final PerformanceCounter mPcWrite;

    private static final Logger LOGGER = LogManager.getLogger(SvSampleAnalyser.class);

    public SvSampleAnalyser(final LinxConfig config, DatabaseAccess dbAccess)
    {
        mConfig = config;
        mSampleId = "";

        mAnalyser = new ClusterAnalyser(config);
        mVisWriter = new VisualiserWriter(config.OutputDataPath, config.Output.WriteVisualisationData, config.hasMultipleSamples());

        mAllVariants = Lists.newArrayList();
        mCnDataLoader = null;

        mIsValid = true;

        mSvFileWriter = null;
        mLinksFileWriter = null;
        mClusterFileWriter = null;

        mFragileSiteAnnotator = new FragileSiteAnnotator();
        mFragileSiteAnnotator.loadFragileSitesFile(mConfig.FragileSiteFile);

        mLineElementAnnotator = new LineElementAnnotator();
        mLineElementAnnotator.loadLineElementsFile(mConfig.LineElementFile);
        mAnalyser.setLineAnnotator(mLineElementAnnotator);

        mPseudoGeneFinder = new PseudoGeneFinder(mVisWriter);
        mLineElementAnnotator.setPseudoGeneFinder(mPseudoGeneFinder);

        mReplicationOriginAnnotator = new ReplicationOriginAnnotator();
        mReplicationOriginAnnotator.loadReplicationOrigins(mConfig.ReplicationOriginsFile);

        mViralInsertAnnotator = new ViralInsertAnnotator();
        mViralInsertAnnotator.loadViralHostData(mConfig.ViralHostsFile);

        mKataegisAnnotator = new KataegisAnnotator(mConfig.OutputDataPath);
        mKataegisAnnotator.loadKataegisData(mConfig.KataegisFile);

        if(config.IndelAnnotation)
        {
            mIndelAnnotator = new IndelAnnotator(dbAccess, config);
        }

        mPcPrep = new PerformanceCounter("Preparation");
        mPcClusterAnalyse = new PerformanceCounter("ClusterAndAnalyse");
        mPcAnnotation = new PerformanceCounter("Annotation");
        mPcWrite = new PerformanceCounter("WriteCSV");

        createOutputFiles();
    }

    public final List<SvVarData> getVariants() { return mAllVariants; }
    public final List<SvCluster> getClusters() { return mAnalyser.getClusters(); }
    public boolean inValidState() { return mIsValid; }
    public final Map<String, List<SvBreakend>> getChrBreakendMap() { return mAnalyser.getState().getChrBreakendMap(); }
    public final VisualiserWriter getVisWriter() { return mVisWriter; }

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

    public void setSampleSVs(final String sampleId, final List<SvVarData> variants)
    {
        clearState();

        if (variants.isEmpty())
            return;

        mSampleId = sampleId;
        mAllVariants.addAll(variants);
        mVisWriter.setSampleId(sampleId);

        if(!mConfig.IsGermline)
        {
            // look-up and cache relevant CN data into each SV
            setSvCopyNumberData(
                    mAllVariants,
                    mCnDataLoader.getSvJcnCalcMap(),
                    mCnDataLoader.getSvIdCnDataMap(),
                    mCnDataLoader.getChrCnDataMap());
        }

        LOGGER.debug("loaded {} SVs", mAllVariants.size());
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
            final List<SvVarData> svList, final EnsemblDataCache ensemblDataCache,
            boolean applyPromotorDistance, boolean selectiveLoading)
    {
        int upstreamDistance = applyPromotorDistance ? PRE_GENE_PROMOTOR_DISTANCE : 0;

        if (selectiveLoading)
        {
            // only load transcript info for the genes covered
            List<String> restrictedGeneIds = Lists.newArrayList();

            for (final SvVarData var : svList)
            {
                for (int be = SE_START; be <= SE_END; ++be)
                {
                    if (be == SE_END && var.isSglBreakend())
                        continue;

                    boolean isStart = isStart(be);

                    ensemblDataCache.populateGeneIdList(restrictedGeneIds, var.chromosome(isStart), var.position(isStart), upstreamDistance);
                }
            }

            ensemblDataCache.loadTranscriptData(restrictedGeneIds);
        }

        // associate breakends with transcripts
        for (final SvVarData var : svList)
        {
            for (int be = SE_START; be <= SE_END; ++be)
            {
                if (be == SE_END && var.isSglBreakend())
                    continue;

                boolean isStart = isStart(be);

                List<GeneAnnotation> genesList = ensemblDataCache.findGeneAnnotationsBySv(
                        var.id(), isStart, var.chromosome(isStart), var.position(isStart), var.orientation(isStart), upstreamDistance);

                if (genesList.isEmpty())
                    continue;

                for (GeneAnnotation gene : genesList)
                {
                    gene.setSvData(var.getSvData());
                }

                var.setGenesList(genesList, isStart);
            }
        }
    }

    public void analyse()
    {
        if(mAllVariants.isEmpty())
            return;

        LOGGER.debug("sample({}) analysing {} variants", mSampleId, mAllVariants.size());

        mPcPrep.start();

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

        List<LinxSvData> linxSvData = prepareSampleData ? Lists.newArrayList() : null;
        List<LinxCluster> clusterData = prepareSampleData ? Lists.newArrayList() : null;
        List<LinxLink> linksData = prepareSampleData ? Lists.newArrayList() : null;
        List<LinxViralInsertFile> viralInserts = prepareSampleData ? generateViralInserts() : null;

        generateSvDataOutput(linxSvData);
        generateClusterOutput(clusterData);
        generateLinksOutput(linksData);

        mVisWriter.writeOutput(mAnalyser.getClusters(), mAllVariants, mCnDataLoader.getChrCnDataMap());

        if(mConfig.isSingleSample())
        {
            try
            {
                // write per-sample DB-style output
                LinxSvDataFile.write(LinxSvDataFile.generateFilename(mConfig.OutputDataPath, mSampleId), linxSvData);
                LinxClusterFile.write(LinxClusterFile.generateFilename(mConfig.OutputDataPath, mSampleId), clusterData);
                LinxLinkFile.write(LinxLinkFile.generateFilename(mConfig.OutputDataPath, mSampleId), linksData);
                LinxViralInsertFile.write(LinxViralInsertFile.generateFilename(mConfig.OutputDataPath, mSampleId), viralInserts);

            } catch (IOException e)
            {
                LOGGER.error("failed to write sample SV data: {}", e.toString());
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

    public void writeSampleWithNoSVs(final String sampleId)
    {
        try
        {
            LinxSvDataFile.write(LinxSvDataFile.generateFilename(mConfig.OutputDataPath, sampleId), Lists.newArrayList());
            LinxClusterFile.write(LinxClusterFile.generateFilename(mConfig.OutputDataPath, sampleId), Lists.newArrayList());
            LinxLinkFile.write(LinxLinkFile.generateFilename(mConfig.OutputDataPath, sampleId), Lists.newArrayList());
            LinxViralInsertFile.write(LinxViralInsertFile.generateFilename(mConfig.OutputDataPath, sampleId), Lists.newArrayList());

            LinxFusionFile.write(LinxFusionFile.generateFilename(mConfig.OutputDataPath, sampleId), Lists.newArrayList());
            LinxBreakendFile.write(LinxBreakendFile.generateFilename(mConfig.OutputDataPath, sampleId), Lists.newArrayList());
            LinxDriverFile.write(LinxDriverFile.generateFilename(mConfig.OutputDataPath, sampleId), Lists.newArrayList());

            ReportableDisruptionFile.write(ReportableDisruptionFile.generateFilename(mConfig.OutputDataPath, sampleId), Lists.newArrayList());
            ReportableGeneFusionFile.write(ReportableGeneFusionFile.generateFilename(mConfig.OutputDataPath, sampleId), Lists.newArrayList());
        }
        catch (IOException e)
        {
            LOGGER.error("failed to write sample SV data: {}", e.toString());
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

            var.setLineElement(mLineElementAnnotator.isLineElement(var, true), true);
            var.setLineElement(mLineElementAnnotator.isLineElement(var, false), false);

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

        cnSegmentBuilder.createCopyNumberData(mCnDataLoader, mAnalyser.getState().getChrBreakendMap());

        cnSegmentBuilder.setSamplePurity(mCnDataLoader, 1, 2, MALE);

        mCnDataLoader.createChrCopyNumberMap();

        setSvCopyNumberData(
                mAllVariants,
                mCnDataLoader.getSvJcnCalcMap(),
                mCnDataLoader.getSvIdCnDataMap(),
                mCnDataLoader.getChrCnDataMap());
    }

    private void createOutputFiles()
    {
        if(!mConfig.hasMultipleSamples())
            return;

        if(mSvFileWriter != null || mClusterFileWriter != null || mLinksFileWriter != null)
            return;

        // open and write headers for multi-sample output files
        createSvDataFile();
        createClusterFile();
        createLinksFile();
    }

    private void createSvDataFile()
    {
        try
        {
            String outputFileName = mConfig.OutputDataPath + "LNX_SVS.csv";

            mSvFileWriter = createBufferedWriter(outputFileName, false);

            // definitional fields
            mSvFileWriter.write("SampleId,Id,Type,ClusterId,ClusterCount");
            mSvFileWriter.write(",ChrStart,PosStart,OrientStart,ArmStart,ChrEnd,PosEnd,OrientEnd,ArmEnd");

            // position and copy number
            mSvFileWriter.write(",CNStart,CNChgStart,CNEnd,CNChgEnd,Jcn,JcnMin,JcnMax");

            // cluster info
            mSvFileWriter.write(",ClusterReason,ClusterDesc,ResolvedType");

            mSvFileWriter.write(",FSStart,FSEnd,LEStart,LEEnd");

            // linked pair info
            mSvFileWriter.write(",LnkSvStart,LnkLenStart,LnkSvEnd,LnkLenEnd,AsmbStart,AsmbEnd");

            // chain info
            mSvFileWriter.write(",ChainId,ChainCount,ChainIndex");

            // proximity info and other link info
            mSvFileWriter.write(",NearestLen,NearestType,DBLenStart,DBLenEnd");

            // proximity info and other link info
            mSvFileWriter.write(",FoldbackLnkStart,FoldbackLenStart,FoldbackInfoStart,FoldbackLnkEnd,FoldbackLenEnd,FoldbackInfoEnd");

            // local topology from arm cluster
            mSvFileWriter.write(",LocTopIdStart,LocTopTypeStart,LocTopTIStart,LocTopIdEnd,LocTopTypeEnd,LocTopTIEnd");

            // gene & replication info
            mSvFileWriter.write(",GeneStart,GeneEnd,RepOriginStart,RepOriginEnd,VirusName");

            if(mConfig.Output.WriteSvData)
            {
                // extra copy number info
                mSvFileWriter.write(",MinorAPStartPrev,MinorAPStartPost,MinorAPEndPrev,MinorAPEndPost,AFStart,AFEnd");

                // SV table info
                mSvFileWriter.write(",HomologyStart,HomologyEnd,InsertSeq,Imprecise,QualScore");
                mSvFileWriter.write(",RefContextStart,RefContextEnd,InsSeqAlignments");
                mSvFileWriter.write(",Recovered,RepeatClass,RepeatType,AnchorStart,AnchorEnd");
                mAnalyser.writeComponentSvHeaders(mSvFileWriter);
            }

            mSvFileWriter.newLine();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to open and write output file headers");
        }
    }

    private static final int INF_DB_MARKER = -2000;

    private void generateSvDataOutput(@Nullable  List<LinxSvData> linxSvData)
    {
        try
        {
            for(final SvVarData var : mAllVariants)
            {
                final SvCluster cluster = var.getCluster();

                if(cluster == null)
                {
                    LOGGER.error("SV({}) not assigned to any cluster", var.posId());
                    continue;
                }

                final StructuralVariantData dbData = var.getSvData();

                final SvArmCluster armClusterStart = cluster.findArmCluster(var.getBreakend(true));

                final SvArmCluster armClusterEnd = !var.isSglBreakend() ? cluster.findArmCluster(var.getBreakend(false)) : null;

                if(mSvFileWriter != null)
                {
                    mSvFileWriter.write(String.format("%s,%d,%s,%d,%d",
                            mSampleId, var.id(), var.typeStr(), cluster.id(), cluster.getSvCount()));

                    mSvFileWriter.write(String.format(",%s,%d,%d,%s,%s,%d,%d,%s",
                            var.chromosome(true), var.position(true), var.orientation(true), asStr(var.arm(true)),
                            var.chromosome(false), var.position(false), var.orientation(false), asStr(var.arm(false))));

                    mSvFileWriter.write(String.format(",%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f",
                            dbData.adjustedStartCopyNumber(), dbData.adjustedStartCopyNumberChange(),
                            dbData.adjustedEndCopyNumber(), dbData.adjustedEndCopyNumberChange(),
                            dbData.junctionCopyNumber(), var.jcnMin(), var.jcnMax()));

                    mSvFileWriter.write(String.format(",%s,%s,%s",
                            var.getClusterReason(), cluster.getDesc(), cluster.getResolvedType()));

                    mSvFileWriter.write(String.format(",%s,%s,%s,%s",
                            var.isFragileSite(true), var.isFragileSite(false),
                            var.getLineElement(true), var.getLineElement(false)));

                    // linked pair info
                    for (int be = SE_START; be <= SE_END; ++be)
                    {
                        boolean isStart = isStart(be);
                        final SvLinkedPair link = var.getLinkedPair(isStart);
                        if (link != null)
                        {
                            mSvFileWriter.write(String.format(",%s,%d",
                                    link.first() == var ? link.second().id() : link.first().id(), link.length()));
                        }
                        else
                        {
                            mSvFileWriter.write(",,-1");
                        }
                    }

                    // assembly info
                    mSvFileWriter.write(String.format(",%s,%s",
                            var.getAssemblyData(true), var.getAssemblyData(false)));

                    // chain info
                    final SvChain chain = cluster.findChain(var);
                    String chainStr = "";

                    if (chain != null)
                    {
                        chainStr = String.format(",%d,%d,%s", chain.id(), chain.getSvCount(), chain.getSvIndices(var));
                    }
                    else
                    {
                        chainStr = String.format(",%d,0,", cluster.getChainId(var));
                    }

                    mSvFileWriter.write(chainStr);

                    // only log DB lengths if the partner is in the cluster
                    final SvLinkedPair dbStart = var.getDBLink(true);
                    final SvLinkedPair dbEnd = var.getDBLink(false);

                    int dbLenStart = NO_DB_MARKER;
                    int dbLenEnd = NO_DB_MARKER;

                    if(dbStart != null && cluster.getSVs().contains(dbStart.getOtherSV(var)))
                    {
                        dbLenStart = (!var.isInferredSgl() && !dbStart.getOtherSV(var).isInferredSgl()) ? dbStart.length() : INF_DB_MARKER;
                    }

                    if(dbEnd != null && cluster.getSVs().contains(dbEnd.getOtherSV(var)))
                    {
                        dbLenEnd = (!var.isInferredSgl() && !dbEnd.getOtherSV(var).isInferredSgl()) ? dbEnd.length() : INF_DB_MARKER;
                    }

                    mSvFileWriter.write(String.format(",%d,%s,%d,%d",
                            var.getNearestSvDistance(), var.getNearestSvRelation(), dbLenStart, dbLenEnd));

                    mSvFileWriter.write(String.format(",%d,%d,%s,%d,%d,%s",
                            var.getFoldbackId(true), var.getFoldbackLength(true), var.getFoldbackInfo(true),
                            var.getFoldbackId(false), var.getFoldbackLength(false), var.getFoldbackInfo(false)));

                    for (int be = SE_START; be <= SE_END; ++be)
                    {
                        SvArmCluster armCluster = be == SE_START ? armClusterStart : armClusterEnd;

                        if (armCluster != null)
                            mSvFileWriter.write(String.format(",%d,%s,%d", armCluster.id(), armCluster.getTypeStr(), armCluster.getTICount()));
                        else
                            mSvFileWriter.write(",-1,,0");
                    }

                    String virusName = "";

                    if(mViralInsertAnnotator != null)
                    {
                        final String[] viralInsertData = mViralInsertAnnotator.matchesViralInsert(var);

                        if (viralInsertData != null)
                            virusName = viralInsertData[VH_NAME];
                    }

                    mSvFileWriter.write(String.format(",%s,%s,%.4f,%.4f,%s",
                            var.getGeneInBreakend(true, true), var.getGeneInBreakend(false, true),
                            var.getReplicationOrigin(true), var.getReplicationOrigin(false), virusName));

                    if(mConfig.Output.WriteSvData)
                    {
                        mSvFileWriter.write(String.format(",%.2f,%.2f,%.2f,%.2f,%.2f,%.2f",
                                var.getBreakend(true).minorAlleleJcn(true),
                                var.getBreakend(true).minorAlleleJcn(false),
                                !var.isSglBreakend() ? var.getBreakend(false).minorAlleleJcn(true) : 0,
                                !var.isSglBreakend() ? var.getBreakend(false).minorAlleleJcn(false) : 0,
                                dbData.adjustedStartAF(), dbData.adjustedEndAF()));

                        final String insSeqAlignments = dbData.insertSequenceAlignments().replaceAll(",", ";");

                        mSvFileWriter.write(String.format(",%s,%s,%s,%s,%.0f,%s,%s,%s",
                                dbData.startHomologySequence(), dbData.endHomologySequence(),
                                dbData.insertSequence(), dbData.imprecise(), dbData.qualityScore(),
                                dbData.startRefContext(), dbData.endRefContext(), insSeqAlignments));

                        mSvFileWriter.write(String.format(",%s,%s,%s,%d,%d",
                                dbData.recovered(), dbData.insertSequenceRepeatClass(), dbData.insertSequenceRepeatType(),
                                dbData.startAnchoringSupportDistance(), dbData.endAnchoringSupportDistance()));

                        mAnalyser.writeComponentSvData(mSvFileWriter, var);
                    }

                    mSvFileWriter.newLine();
                }

                if(linxSvData != null)
                {
                    linxSvData.add(ImmutableLinxSvData.builder()
                            .svId(var.id())
                            .clusterId(cluster.id())
                            .clusterReason(var.getClusterReason())
                            .fragileSiteStart(var.isFragileSite(true))
                            .fragileSiteEnd(var.isFragileSite(false))
                            .isFoldback(var.isFoldback())
                            .lineTypeStart(var.getLineElement(true))
                            .lineTypeEnd(var.getLineElement(false))
                            .junctionCopyNumberMin(var.jcnMin())
                            .junctionCopyNumberMax(var.jcnMax())
                            .geneStart(var.getGeneInBreakend(true, false))
                            .geneEnd(var.getGeneInBreakend(true, false))
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
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    private void createClusterFile()
    {
        try
        {
            String outputFileName = mConfig.OutputDataPath + "LNX_CLUSTERS.csv";

            mClusterFileWriter = createBufferedWriter(outputFileName, false);

            mClusterFileWriter.write("SampleId,ClusterId,ClusterDesc,ClusterCount,SuperType,ResolvedType,Synthetic,Subclonal,FullyChained,ChainCount");
            mClusterFileWriter.write(",DelCount,DupCount,InsCount,InvCount,BndCount,SglCount,InfCount");
            mClusterFileWriter.write(",ClusterReasons,Consistency,IsLINE,Replication,MinJcn,MaxJcn,Foldbacks");
            mClusterFileWriter.write(",ArmCount,OriginArms,FragmentArms,ConsistentArms,ComplexArms,Annotations,AlleleValidPerc");

            mClusterFileWriter.write(",TotalTIs,AssemblyTIs,ShortTIs,IntTIs,ExtTIs,IntShortTIs,ExtShortTIs,IntTIsCnGain");
            mClusterFileWriter.write(",ExtTIsCnGain,OverlapTIs,ChainEndsFace,ChainEndsAway,UnchainedSVs");

            mClusterFileWriter.write(",DBs,ShortDBs,TotalDBLength,TotalDeleted,TravDelCount,TravDelLength");
            mClusterFileWriter.write(",TotalRange,ChainedLength,ImpliedTIs");

            mClusterFileWriter.write(",ArmClusterCount,AcTotalTIs,AcIsolatedBE,AcTIOnly,AcDsb,AcSimpleDup");
            mClusterFileWriter.write(",AcSingleFb,AcFbDsb,AcComplexFb,AcComplexLine,AcSameOrient,AcComplexOther");

            if(mIndelAnnotator != null)
                mClusterFileWriter.write(",IndelCount,IndelProb");

            mClusterFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing cluster-data to outputFile: {}", e.toString());
        }
    }

    private void generateClusterOutput(@Nullable List<LinxCluster> clusterData)
    {
        try
        {
            for(final SvCluster cluster : getClusters())
            {
                int clusterSvCount = cluster.getSvCount();

                if(clusterSvCount == 1 && !mConfig.Output.WriteSingleSVClusters)
                    continue;

                ResolvedType resolvedType = cluster.getResolvedType();

                final String superType = getSuperType(cluster);

                if(mClusterFileWriter != null)
                {
                    mClusterFileWriter.write(String.format("%s,%d,%s,%d,%s,%s,%s,%s,%s,%d",
                            mSampleId, cluster.id(), cluster.getDesc(), clusterSvCount,
                            superType, resolvedType, cluster.isSyntheticType(),
                            cluster.isSubclonal(), cluster.isFullyChained(false), cluster.getChains().size()));

                    mClusterFileWriter.write(String.format(",%d,%d,%d,%d,%d,%d,%d",
                            cluster.getTypeCount(DEL), cluster.getTypeCount(DUP), cluster.getTypeCount(INS),
                            cluster.getTypeCount(INV), cluster.getTypeCount(BND), cluster.getTypeCount(SGL), cluster.getTypeCount(INF)));

                    double foldbackCount = 0;

                    for (final SvVarData var : cluster.getFoldbacks())
                    {
                        // avoid double-count chained foldbacks
                        if (var.getFoldbackBreakend(true) != null)
                            foldbackCount += 0.5;
                        if (var.getFoldbackBreakend(false) != null)
                            foldbackCount += 0.5;
                    }

                    mClusterFileWriter.write(String.format(",%s,%d,%s,%s,%.1f,%.1f,%.0f",
                            cluster.getClusteringReasons(), cluster.getConsistencyCount(), cluster.hasLinkingLineElements(),
                            cluster.requiresReplication(), cluster.getMinJcn(), cluster.getMaxJcn(), foldbackCount));

                    final ClusterMetrics metrics = cluster.getMetrics();

                    mClusterFileWriter.write(String.format(",%d,%d,%d,%d,%d,%s,%.2f",
                            cluster.getArmCount(), cluster.getOriginArms(), cluster.getFragmentArms(), cluster.getConsistentArms(),
                            cluster.getComplexArms(), cluster.getAnnotations(), metrics.ValidAlleleJcnSegmentPerc));

                    long shortTIs = cluster.getLinkedPairs().stream().filter(x -> x.length() <= SHORT_TI_LENGTH).count();

                    mClusterFileWriter.write(String.format(",%d,%d,%d",
                            cluster.getLinkedPairs().size(), cluster.getAssemblyLinkedPairs().size(), shortTIs));

                    final ChainMetrics chainMetrics = cluster.getLinkMetrics();

                    mClusterFileWriter.write(String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
                            chainMetrics.InternalTIs, chainMetrics.ExternalTIs, chainMetrics.InternalShortTIs, chainMetrics.ExternalShortTIs,
                            chainMetrics.InternalTICnGain, chainMetrics.ExternalTICnGain, chainMetrics.OverlappingTIs,
                            chainMetrics.ChainEndsFace, chainMetrics.ChainEndsAway, cluster.getUnlinkedSVs().size()));

                    mClusterFileWriter.write(String.format(",%d,%d,%d,%d,%d,%d",
                            metrics.DBCount, metrics.ShortDBCount, metrics.TotalDBLength,
                            metrics.TotalDeleted, metrics.TraversedDelCount, metrics.TraversedDelLength));

                    mClusterFileWriter.write(String.format(",%d,%d,%d",
                            metrics.TotalRange, metrics.ChainedLength, metrics.ImpliedTICount));

                    final int[] armClusterData = getArmClusterData(cluster);
                    long armClusterTIs = cluster.getArmClusters().stream().mapToInt(x -> x.getTICount()).sum();

                    mClusterFileWriter.write(String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
                            cluster.getArmClusters().size(), armClusterTIs, armClusterData[ARM_CL_ISOLATED_BE],
                            armClusterData[ARM_CL_TI_ONLY], armClusterData[ARM_CL_DSB], armClusterData[ARM_CL_SIMPLE_DUP],
                            armClusterData[ARM_CL_FOLDBACK], armClusterData[ARM_CL_FOLDBACK_DSB], armClusterData[ARM_CL_COMPLEX_FOLDBACK],
                            armClusterData[ARM_CL_COMPLEX_LINE], armClusterData[ARM_CL_SAME_ORIENT], armClusterData[ARM_CL_COMPLEX_OTHER]));

                    if(mIndelAnnotator != null)
                    {
                        mClusterFileWriter.write(String.format(",%d,%f", metrics.IndelCount, metrics.IndelProbability));
                    }

                    mClusterFileWriter.newLine();
                }

                if(clusterData != null)
                {
                    clusterData.add(ImmutableLinxCluster.builder()
                            .clusterId(cluster.id())
                            .resolvedType(superType)
                            .synthetic(cluster.isSyntheticType())
                            .subClonal(cluster.isSubclonal())
                            .subType(cluster.getResolvedType().toString())
                            .clusterCount(clusterSvCount)
                            .clusterDesc(cluster.getDesc())
                            .build());
                }
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing cluster-data to outputFile: {}", e.toString());
        }
    }

    private void createLinksFile()
    {
        if(!mConfig.Output.WriteLinks)
            return;

        try
        {
            String outputFileName = mConfig.OutputDataPath + "LNX_LINKS.csv";

            mLinksFileWriter = createBufferedWriter(outputFileName, false);

            mLinksFileWriter.write("SampleId,ClusterId,ClusterCount,ResolvedType");
            mLinksFileWriter.write(",ChainId,ChainCount,ChainConsistent,LinkReason,LinkIndex,ChainIndex,Jcn,JcnUncertainty");
            mLinksFileWriter.write(",IsAssembled,TILength,NextSvDist,NextClusteredSvDist,TraversedSVCount");
            mLinksFileWriter.write(",LocationType,OverlapCount,CopyNumberGain");
            mLinksFileWriter.write(",Id1,Id2,ChrArm,PosStart,PosEnd,LocTopTypeStart,LocTopTypeEnd,GeneStart,GeneEnd,ExonMatch");
            mLinksFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing links to outputFile: {}", e.toString());
        }
    }

    private void generateLinksOutput(@Nullable List<LinxLink> linksData)
    {
        if(!mConfig.Output.WriteLinks)
            return;

        try
        {
            for(final SvCluster cluster : getClusters())
            {
                int clusterSvCount = cluster.getSvCount();

                List<SvChain> chains = cluster.getChains();

                for (final SvChain chain : chains)
                {
                    int chainSvCount = chain.getSvCount();
                    boolean chainConsistent = chain.isConsistent();

                    List<SvLinkedPair> uniquePairs = Lists.newArrayList();
                    final List<SvLinkedPair> chainLinks = chain.getLinkedPairs();

                    for (int chainIndex = 0; chainIndex < chainLinks.size(); ++chainIndex)
                    {
                        final SvLinkedPair pair = chainLinks.get(chainIndex);

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

                        if(mLinksFileWriter != null)
                        {
                            mLinksFileWriter.write(String.format("%s,%d,%d,%s",
                                    mSampleId, cluster.id(), clusterSvCount, cluster.getResolvedType()));

                            mLinksFileWriter.write(String.format(",%d,%d,%s,%s,%d,%s,%s,%.3f",
                                    chain.id(), chainSvCount, chainConsistent, pair.getLinkReason(), pair.getLinkIndex(),
                                    chainIndexStr, formatJcn(chain.jcn()), chain.jcnUncertainty()));

                            mLinksFileWriter.write(String.format(",%s,%d,%d,%d,%d,%s,%d,%s",
                                    pair.isAssembled(), pair.length(),
                                    pair.getNextSvDistance(), pair.getNextClusteredSvDistance(), pair.getTraversedSVCount(),
                                    pair.locationType(), pair.overlapCount(), pair.hasCopyNumberGain()));

                            SvArmCluster acStart = cluster.findArmCluster(beStart);
                            SvArmCluster acEnd = cluster.findArmCluster(beEnd);

                            mLinksFileWriter.write(String.format(",%d,%d,%s,%d,%d,%s,%s,%s,%s,%s",
                                    beStart.getSV().id(), beEnd.getSV().id(),
                                    beStart.getChrArm(), beStart.position(), beEnd.position(),
                                    acStart != null ? acStart.getTypeStr() : "", acEnd != null ? acEnd.getTypeStr() : "",
                                    beStart.getSV().getGeneInBreakend(beStart.usesStart(), false),
                                    beEnd.getSV().getGeneInBreakend(beEnd.usesStart(), false), pair.getExonMatchData()));

                            mLinksFileWriter.newLine();
                        }

                        if(linksData != null)
                        {
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
                                    .build());
                        }
                    }
                }
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing links to outputFile: {}", e.toString());
        }
    }

    private List<LinxViralInsertFile> generateViralInserts()
    {
        List<LinxViralInsertFile> viralInserts = Lists.newArrayList();

        for(final SvVarData var : mAllVariants)
        {
            final String[] viralInsertData = mViralInsertAnnotator.matchesViralInsert(var);

            if(viralInsertData != null)
            {
                viralInserts.add(new LinxViralInsertFile(mSampleId, var.id(), viralInsertData[VH_ID], viralInsertData[VH_NAME]));
            }
        }

        return viralInserts;
    }

    public void close()
    {
        if(mConfig.hasMultipleSamples() || LOGGER.isDebugEnabled())
        {
            // log perf stats
            mPcPrep.logStats();
            mPcClusterAnalyse.logStats();
            mAnalyser.logStats();
            mPcAnnotation.logStats();
            mPcWrite.logStats();
        }

        closeBufferedWriter(mSvFileWriter);
        closeBufferedWriter(mClusterFileWriter);
        closeBufferedWriter(mLinksFileWriter);
        mKataegisAnnotator.close();
        mVisWriter.close();
        mAnalyser.close();

        if(mIndelAnnotator != null)
            mIndelAnnotator.close();
    }
}
