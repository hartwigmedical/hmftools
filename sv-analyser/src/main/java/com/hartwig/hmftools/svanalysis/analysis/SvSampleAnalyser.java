package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.NO_DB_MARKER;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_LINE;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_OTHER;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_FOLDBACK;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_FOLDBACK_DSB;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_FOLDBACK_PAIR_FACING;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_FOLDBACK_PAIR_OPPOSING;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_FOLDBACK_PAIR_SAME;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_FOLDBACK_TI;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_MULTIPLE_DSBS;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_REMOTE_TI;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_SINGLE;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.getArmClusterData;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_EXT_TI_CN_GAIN;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_DB;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_SHORT_DB;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_INT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_INT_TI_CN_GAIN;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isSpecificSV;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.ReplicationOriginAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.VisualiserWriter;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvaConfig;
import com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class SvSampleAnalyser {

    private final SvaConfig mConfig;
    private final ClusterAnalyser mAnalyser;

    // data per run (ie sample)
    private String mSampleId;
    private List<SvVarData> mAllVariants; // the original list to analyse

    private BufferedWriter mSvFileWriter;
    private BufferedWriter mClusterFileWriter;
    private BufferedWriter mLinksFileWriter;
    private VisualiserWriter mVisWriter;

    private FragileSiteAnnotator mFragileSiteAnnotator;
    private LineElementAnnotator mLineElementAnnotator;
    private ReplicationOriginAnnotator mReplicationOriginAnnotator;
    private SvClusteringMethods mClusteringMethods;
    private CNAnalyser mCopyNumberAnalyser;

    private PerformanceCounter mPcPrep;
    private PerformanceCounter mPcClusterAnalyse;
    private PerformanceCounter mPcWrite;

    private static final Logger LOGGER = LogManager.getLogger(SvSampleAnalyser.class);

    public SvSampleAnalyser(final SvaConfig config)
    {
        mConfig = config;
        mSampleId = "";

        mClusteringMethods = new SvClusteringMethods(mConfig.ProximityDistance);
        mAnalyser = new ClusterAnalyser(config, mClusteringMethods);
        mVisWriter = new VisualiserWriter(config.OutputCsvPath, config.WriteVisualisationData);

        mAnalyser.setUseAllelePloidies(true);

        mAllVariants = Lists.newArrayList();
        mCopyNumberAnalyser = null;

        mSvFileWriter = null;
        mLinksFileWriter = null;
        mClusterFileWriter = null;

        mFragileSiteAnnotator = new FragileSiteAnnotator();
        mFragileSiteAnnotator.loadFragileSitesFile(mConfig.FragileSiteFile);

        mLineElementAnnotator = new LineElementAnnotator();
        mLineElementAnnotator.loadLineElementsFile(mConfig.LineElementFile);

        mReplicationOriginAnnotator = new ReplicationOriginAnnotator();
        mReplicationOriginAnnotator.loadReplicationOrigins(mConfig.ReplicationOriginsFile);

        mPcPrep = new PerformanceCounter("Preparation");
        mPcClusterAnalyse = new PerformanceCounter("ClusterAndAnalyse");
        mPcWrite = new PerformanceCounter("WriteCSV");
    }

    public final List<SvVarData> getVariants() { return mAllVariants; }
    public final List<SvCluster> getClusters() { return mAnalyser.getClusters(); }
    public final Map<String, List<SvBreakend>> getChrBreakendMap() { return mClusteringMethods.getChrBreakendMap(); }
    public final VisualiserWriter getVisWriter() { return mVisWriter; }
    public void setCopyNumberAnalyser(CNAnalyser cnAnalyser)
    {
        mCopyNumberAnalyser = cnAnalyser;
        mAnalyser.setCopyNumberAnalyser(cnAnalyser);
        mClusteringMethods.setSampleLohData(mCopyNumberAnalyser.getSampleLohData());
        mClusteringMethods.setChrCopyNumberMap(mCopyNumberAnalyser.getChrCopyNumberMap());
    }

    public void setGeneCollection(SvGeneTranscriptCollection geneCollection) { mAnalyser.setGeneCollection(geneCollection); }

    private void clearState()
    {
        mClusteringMethods.clearLOHBreakendData(mSampleId);

        // reduce maps with already processed sample data for the larger data sources
        if(!mSampleId.isEmpty())
        {
            if (mCopyNumberAnalyser.getSampleSvPloidyCalcMap() != null)
                mCopyNumberAnalyser.getSampleSvPloidyCalcMap().remove(mSampleId); // shrink the data source to make future look-ups faster

            if (mClusteringMethods.getSampleLohData() != null)
                mClusteringMethods.getSampleLohData().remove(mSampleId);
        }

        mSampleId = "";
        mAllVariants.clear();
    }

    public void loadFromDatabase(final String sampleId, final List<SvVarData> variants)
    {
        clearState();

        if (variants.isEmpty())
            return;

        mSampleId = sampleId;
        mAllVariants = Lists.newArrayList(variants);
        mVisWriter.setSampleId(sampleId);

        // look-up and cache relevant CN data into each SV
        final Map<String,double[]> svPloidyCalcDataMap = mCopyNumberAnalyser.getSampleSvPloidyCalcMap().get(mSampleId);
        final Map<String,SvCNData[]> svIdCnDataMap = mCopyNumberAnalyser.getSvIdCnDataMap();
        final Map<String,List<SvCNData>> chrCnDataMap = mCopyNumberAnalyser.getChrCnDataMap();

        setSvCopyNumberData(
                mAllVariants,
                mCopyNumberAnalyser.getSampleSvPloidyCalcMap().get(mSampleId),
                mCopyNumberAnalyser.getSvIdCnDataMap(),
                mCopyNumberAnalyser.getChrCnDataMap());

        LOGGER.debug("loaded {} SVs", mAllVariants.size());
    }

    public static void setSvCopyNumberData(List<SvVarData> svList, final Map<String,double[]> svPloidyCalcDataMap,
            final Map<String,SvCNData[]> svIdCnDataMap, final Map<String,List<SvCNData>> chrCnDataMap)
    {
        if((svPloidyCalcDataMap == null || svPloidyCalcDataMap.isEmpty()) && svIdCnDataMap.isEmpty())
            return;

        List<SvCNData> cnDataList = null;
        String currentChromosome = "";
        for(final SvVarData var : svList)
        {
            if(svPloidyCalcDataMap != null)
            {
                final double[] ploidyData = svPloidyCalcDataMap.get(var.id());
                if (ploidyData != null)
                {
                    double estPloidy = ploidyData[0];
                    double estUncertainty = ploidyData[1];
                    var.setPloidyRecalcData(estPloidy - estUncertainty, estPloidy + estUncertainty);
                }
            }

            final SvCNData[] cnDataPair = svIdCnDataMap.get(var.id());

            if(cnDataPair == null)
                continue;

            for(int be = SVI_START; be <= SVI_END; ++be)
            {
                if(var.isNullBreakend() && be == SVI_END)
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

    public void analyse()
    {
        if(mAllVariants.isEmpty())
            return;

        LOGGER.debug("sample({}) clustering {} variants", mSampleId, mAllVariants.size());

        mPcPrep.start();
        annotateAndFilterVariants();
        mClusteringMethods.populateChromosomeBreakendMap(mAllVariants);
        mClusteringMethods.annotateNearestSvData();
        LinkFinder.findDeletionBridges(mClusteringMethods.getChrBreakendMap());
        mClusteringMethods.setSimpleVariantLengths(mSampleId);
        mReplicationOriginAnnotator.setReplicationOrigins(mClusteringMethods.getChrBreakendMap());
        mPcPrep.stop();

        mPcClusterAnalyse.start();

        mAnalyser.setSampleData(mSampleId, mAllVariants);
        mAnalyser.clusterAndAnalyse();

        mPcClusterAnalyse.stop();
    }

    public void writeOutput()
    {
        mPcWrite.start();

        if(!mConfig.OutputCsvPath.isEmpty())
        {
            writeClusterSVOutput();
            writeClusterData();

            if(mConfig.hasMultipleSamples())
            {
                writeClusterLinkData();
            }

            mVisWriter.writeOutput(mAnalyser.getClusters(), mAllVariants);
        }

        mPcWrite.stop();

    }

    private void annotateAndFilterVariants()
    {
        int currentIndex = 0;

        while(currentIndex < mAllVariants.size())
        {
            SvVarData var = mAllVariants.get(currentIndex);

            var.setFragileSites(mFragileSiteAnnotator.isFragileSite(var, true), mFragileSiteAnnotator.isFragileSite(var, false));
            var.setLineElement(mLineElementAnnotator.isLineElement(var, true), true);
            var.setLineElement(mLineElementAnnotator.isLineElement(var, false), false);

            String startArm = getChromosomalArm(var.chromosome(true), var.position(true));

            String endArm;
            if(!var.isNullBreakend())
                endArm = getChromosomalArm(var.chromosome(false), var.position(false));
            else
                endArm = CHROMOSOME_ARM_P;

            var.setChromosomalArms(startArm, endArm);

            ++currentIndex;
        }
    }

    private void writeClusterSVOutput()
    {
        try
        {
            if(mSvFileWriter == null)
            {
                String outputFileName = mConfig.OutputCsvPath;

                if(mConfig.hasMultipleSamples())
                    outputFileName += "SVA_SVS.csv";
                else
                    outputFileName += mSampleId + "_SVA.csv";

                mSvFileWriter = createBufferedWriter(outputFileName, false);

                // definitional fields
                mSvFileWriter.write("SampleId,Id,Type,ClusterId,ClusterCount");
                mSvFileWriter.write(",ChrStart,PosStart,OrientStart,ArmStart,ChrEnd,PosEnd,OrientEnd,ArmEnd");

                // position and copy number
                mSvFileWriter.write(",AFStart,CNStart,CNChgStart,AFEnd,CNEnd,CNChgEnd,Ploidy,PloidyMin,PloidyMax");

                // cluster info
                mSvFileWriter.write(",ClusterReason,ClusterDesc,IsResolved,ResolvedType,Consistency,ArmCount");

                // SV info
                mSvFileWriter.write(",Homology,InexactHOStart,InexactHOEnd,InsertSeq,Imprecise,QualScore");
                mSvFileWriter.write(",RefContextStart,RefContextEnd,InsSeqAlignments,Recovered");

                mSvFileWriter.write(",FSStart,FSEnd,LEStart,LEEnd"); // ,DupBEStart,DupBEEnd

                // linked pair info
                mSvFileWriter.write(",LnkSvStart,LnkLenStart,LnkSvEnd,LnkLenEnd");
                mSvFileWriter.write(",AsmbStart,AsmbEnd,AsmbMatchStart,AsmbMatchEnd");

                // chain info
                mSvFileWriter.write(",ChainId,ChainCount,ChainIndex");

                // proximity info and other link info
                mSvFileWriter.write(",NearestLen,NearestType,DBLenStart,DBLenEnd,SynDelDupLen,SynDelDupTILen");

                // proximity info and other link info
                mSvFileWriter.write(",FoldbackLnkStart,FoldbackLenStart,FoldbackInfoStart,FoldbackLnkEnd,FoldbackLenEnd,FoldbackInfoEnd");

                // gene & replication info
                mSvFileWriter.write(",GeneStart,GeneEnd,RepOriginStart,RepOriginEnd");

                // extra copy number info
                mSvFileWriter.write(",ActBafStartPrev,ActBafStartPost,ActBafEndPrev,ActBafEndPost");

                mSvFileWriter.newLine();
            }

            BufferedWriter writer = mSvFileWriter;

            int lineCount = 0;
            int svCount = 0;

            for(final SvVarData var : mAllVariants)
            {
                final SvCluster cluster = var.getCluster();

                if(cluster == null)
                {
                    LOGGER.error("SV({}) not assigned to any cluster", var.posId());
                    continue;
                }

                final StructuralVariantData dbData = var.getSvData();

                ++svCount;

                writer.write(String.format("%s,%s,%s,%d,%d",
                        mSampleId, var.id(), var.typeStr(), cluster.id(), cluster.getSvCount()));

                writer.write(String.format(",%s,%d,%d,%s,%s,%d,%d,%s",
                        var.chromosome(true), var.position(true), var.orientation(true), var.arm(true),
                        var.chromosome(false), var.position(false), var.orientation(false), var.arm(false)));

                writer.write(String.format(",%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f",
                        dbData.adjustedStartAF(), dbData.adjustedStartCopyNumber(), dbData.adjustedStartCopyNumberChange(),
                        dbData.adjustedEndAF(), dbData.adjustedEndCopyNumber(), dbData.adjustedEndCopyNumberChange(),
                        dbData.ploidy(), var.ploidyMin(), var.ploidyMax()));

                writer.write(String.format(",%s,%s,%s,%s,%d,%d",
                        var.getClusterReason(), cluster.getDesc(), cluster.isResolved(), cluster.getResolvedType(),
                        cluster.getConsistencyCount(), cluster.getArmCount()));

                final String insSeqAlignments = dbData.insertSequenceAlignments().replaceAll(",", ";");

                writer.write(String.format(",%s,%d,%d,%s,%s,%.0f,%s,%s,%s,%s",
                        dbData.insertSequence().isEmpty() && var.type() != INS ? dbData.homology() : "",
                        dbData.inexactHomologyOffsetStart(), dbData.inexactHomologyOffsetEnd(),
                        dbData.insertSequence(), dbData.imprecise(), dbData.qualityScore(),
                        dbData.startRefContext(), dbData.endRefContext(), insSeqAlignments,
                        dbData.recovered()));

                writer.write(String.format(",%s,%s,%s,%s",
                        var.isFragileSite(true), var.isFragileSite(false),
                        var.getLineElement(true), var.getLineElement(false)));

                // linked pair info
                final SvLinkedPair startLP = var.getLinkedPair(true);
                String startLinkStr = "0,-1";
                String assemblyMatchStart = var.getAssemblyMatchType(true);
                if(startLP != null)
                {
                    startLinkStr = String.format("%s,%d",
                            startLP.first().equals(var, true) ? startLP.second().origId() : startLP.first().origId(),
                            startLP.length());
                }

                final SvLinkedPair endLP = var.getLinkedPair(false);
                String endLinkStr = "0,-1";
                String assemblyMatchEnd = var.getAssemblyMatchType(false);
                if(endLP != null)
                {
                    endLinkStr = String.format("%s,%d",
                            endLP.first().equals(var, true) ? endLP.second().origId() : endLP.first().origId(),
                            endLP.length());
                }

                // assembly info
                writer.write(String.format(",%s,%s,%s,%s,%s,%s",
                        startLinkStr, endLinkStr, var.getAssemblyData(true), var.getAssemblyData(false),
                        assemblyMatchStart, assemblyMatchEnd));

                // chain info
                final SvChain chain = cluster.findChain(var);
                String chainStr = "";

                if(chain != null)
                {
                    chainStr = String.format(",%d,%d,%s", chain.id(), chain.getSvCount(), chain.getSvIndices(var));
                }
                else
                {
                    chainStr = String.format(",%d,0,", cluster.getChainId(var));
                }

                writer.write(chainStr);

                int dbLenStart = var.getDBLink(true) != null ? var.getDBLink(true).length() : NO_DB_MARKER;
                int dbLenEnd = var.getDBLink(false) != null ? var.getDBLink(false).length() : NO_DB_MARKER;

                writer.write(String.format(",%d,%s,%d,%d,%d,%d",
                        var.getNearestSvDistance(), var.getNearestSvRelation(), dbLenStart, dbLenEnd,
                        cluster.getSynDelDupLength(), cluster.getSynDelDupTILength()));

                writer.write(String.format(",%s,%d,%s,%s,%d,%s",
                        var.getFoldbackLink(true), var.getFoldbackLength(true), var.getFoldbackInfo(true),
                        var.getFoldbackLink(false), var.getFoldbackLength(false), var.getFoldbackInfo(false)));

                writer.write(String.format(",%s,%s,%.4f,%.4f",
                        var.getGeneInBreakend(true), var.getGeneInBreakend(false),
                        var.getReplicationOrigin(true), var.getReplicationOrigin(false)));

                writer.write(String.format(",%.2f,%,2f,%,2f,%,2f",
                        var.getBreakend(true).actualBaf(true), var.getBreakend(true).actualBaf(false),
                        !var.isNullBreakend() ? var.getBreakend(false).actualBaf(true) : 0,
                        !var.isNullBreakend() ? var.getBreakend(false).actualBaf(false) : 0));

                ++lineCount;
                writer.newLine();

                if(svCount != lineCount)
                {
                    LOGGER.error("inconsistent output");
                }
            }

        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    private void writeClusterData()
    {
        try
        {
            if(mClusterFileWriter == null)
            {
                String outputFileName = mConfig.OutputCsvPath;

                outputFileName += "SVA_CLUSTERS.csv";

                mClusterFileWriter = createBufferedWriter(outputFileName, false);

                mClusterFileWriter.write("SampleId,ClusterId,ClusterDesc,ClusterCount,ResolvedType,FullyChained,ChainCount");
                mClusterFileWriter.write(",DelCount,DupCount,InsCount,InvCount,BndCount,SglCount");
                mClusterFileWriter.write(",ClusterReasons,Consistency,ArmCount,OriginArms,FragmentArms,IsLINE,HasReplicated,Foldbacks");
                mClusterFileWriter.write(",IntTIs,ExtTIs,IntTIsWithGain,ExtTIsWithGain,DSBs,ShortDSBs");
                mClusterFileWriter.write(",TotalLinks,AssemblyLinks,LongDelDups,UnlinkedRemotes,ShortTIRemotes,MinCopyNumber,MaxCopyNumber");
                mClusterFileWriter.write(",SynDelDupLen,SynDelDupAvgTILen,Annotations,UnchainedSVs,AlleleValidPerc");
                mClusterFileWriter.write(",ArmClusterCount,AcSoloSv,AcRemoteTI,AcDsb,AcMultipleDsb,AcSingleFb,AcFbTI,AcFbDSB");
                mClusterFileWriter.write(",ArmFbPairSame,ArmFbPairOpp,ArmFbPairFacing,AcComplexFb,AcComplexLine,AcComplexOther");
                mClusterFileWriter.newLine();
            }

            BufferedWriter writer = mClusterFileWriter;

            for(final SvCluster cluster : getClusters())
            {
                int clusterSvCount = cluster.getSvCount();

                writer.write(String.format("%s,%d,%s,%d,%s,%s,%d",
                        mSampleId, cluster.id(), cluster.getDesc(), clusterSvCount, cluster.getResolvedType(),
                        cluster.isFullyChained(), cluster.getChains().size()));

                writer.write(String.format(",%d,%d,%d,%d,%d,%d",
                        cluster.getTypeCount(DEL), cluster.getTypeCount(DUP), cluster.getTypeCount(INS),
                        cluster.getTypeCount(INV), cluster.getTypeCount(BND), cluster.getTypeCount(SGL)));

                writer.write(String.format(",%s,%d,%d,%d,%d,%s,%s,%d",
                        cluster.getClusteringReasons(), cluster.getConsistencyCount(),
                        cluster.getArmCount(), cluster.getOriginArms(), cluster.getFragmentArms(),
                        cluster.hasLinkingLineElements(), cluster.hasReplicatedSVs(), cluster.getFoldbacks().size()));

                isSpecificCluster(cluster);

                int[] chainData = cluster.getLinkMetrics();

                writer.write(String.format(",%d,%d,%d,%d,%d,%d",
                        chainData[CM_INT_TI], chainData[CM_EXT_TI], chainData[CM_INT_TI_CN_GAIN], chainData[CM_EXT_TI_CN_GAIN],
                        chainData[CM_DB], chainData[CM_SHORT_DB]));

                /*
                final String chainInfo = cluster.getChains().stream()
                        .filter(x -> !x.getDetails().isEmpty())
                        .map(SvChain::getDetails)
                        .collect (Collectors.joining (";"));
                */

                writer.write(String.format(",%d,%d,%d,%d,%d,%.2f,%.2f",
                        cluster.getLinkedPairs().size(), cluster.getAssemblyLinkedPairs().size(), cluster.getLongDelDups().size(),
                        cluster.getUnlinkedRemoteSVs().size(), cluster.getShortTIRemoteSVs().size(),
                        cluster.getMinCNChange(), cluster.getMaxCNChange()));

                writer.write(String.format(",%d,%d,%s,%d,%.2f",
                        cluster.getSynDelDupLength(), cluster.getSynDelDupTILength(), cluster.getAnnotations(),
                        cluster.getUnlinkedSVs().size(), cluster.getValidAllelePloidySegmentPerc()));

                // ArmClusterCount,AcSoloSv,AcRemoteTI,AcDsb,AcMultipleDsb,AcSingleFb,AcFbTI,AcFbDSB
                // ArmFbPairSame,ArmFbPairOpp,ArmFbPairFacing,AcComplexFb,AcComplexLine,AcComplexOther

                final int[] armClusterData = getArmClusterData(cluster);

                writer.write(String.format(",%d,%d,%d,%d,%d,%d,%d,%d",
                        cluster.getArmClusters().size(), armClusterData[ARM_CL_SINGLE], armClusterData[ARM_CL_REMOTE_TI],
                        armClusterData[ARM_CL_DSB], armClusterData[ARM_CL_MULTIPLE_DSBS],
                        armClusterData[ARM_CL_FOLDBACK], armClusterData[ARM_CL_FOLDBACK_TI], armClusterData[ARM_CL_FOLDBACK_DSB]));

                writer.write(String.format(",%d,%d,%d,%d,%d,%d",
                        armClusterData[ARM_CL_FOLDBACK_PAIR_SAME], armClusterData[ARM_CL_FOLDBACK_PAIR_OPPOSING],
                        armClusterData[ARM_CL_FOLDBACK_PAIR_FACING], armClusterData[ARM_CL_COMPLEX_FOLDBACK],
                        armClusterData[ARM_CL_COMPLEX_LINE], armClusterData[ARM_CL_COMPLEX_OTHER]));

                writer.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing cluster-data to outputFile: {}", e.toString());
        }
    }

    private void writeClusterLinkData()
    {
        try
        {
            if(mLinksFileWriter == null)
            {
                String outputFileName = mConfig.OutputCsvPath;

                outputFileName += "SVA_LINKS.csv";

                mLinksFileWriter = createBufferedWriter(outputFileName, false);

                mLinksFileWriter.write("SampleId,ClusterId,ClusterDesc,ClusterCount,ResolvedType,IsLINE,FullyChained");
                mLinksFileWriter.write(",ChainId,ChainCount,ChainConsistent,Id1,Id2,ChrArm,IsAssembled,TILength,SynDelDupLen");
                mLinksFileWriter.write(",NextSVDistance,NextSVTraversedCount,DBLenStart,DBLenEnd,OnArmOfOrigin,TraversedSVCount");
                mLinksFileWriter.write(",PosStart,PosEnd,GeneStart,GeneEnd,ExonMatch");
                mLinksFileWriter.newLine();
            }

            BufferedWriter writer = mLinksFileWriter;

            for(final SvCluster cluster : getClusters())
            {
                int clusterSvCount = cluster.getSvCount();

                // isSpecificCluster(cluster);

                List<SvChain> chains = cluster.getChains();

                for (final SvChain chain : chains)
                {
                    int chainSvCount = chain.getSvCount();
                    boolean chainConsistent = chain.isConsistent();

                    List<SvLinkedPair> uniquePairs = Lists.newArrayList();

                    for (final SvLinkedPair pair : chain.getLinkedPairs())
                    {
                        if(pair.linkType() != LINK_TYPE_TI)
                            continue;

                        boolean isRepeat = false;

                        for (final SvLinkedPair existingPair : uniquePairs)
                        {
                            if(pair.matches(existingPair))
                            {
                                isRepeat = true;
                                break;
                            }
                        }

                        if(isRepeat)
                            continue;

                        uniquePairs.add(pair);

                        writer.write(String.format("%s,%d,%s,%d,%s,%s,%s",
                            mSampleId, cluster.id(), cluster.getDesc(), clusterSvCount, cluster.getResolvedType(),
                            cluster.hasLinkingLineElements(), cluster.isFullyChained()));

                        final SvBreakend beStart = pair.getBreakend(true);
                        final SvBreakend beEnd= pair.getBreakend(false);

                        writer.write(String.format(",%d,%d,%s,%s,%s,%s",
                                chain.id(), chainSvCount, chainConsistent,
                                beStart.getSV().origId(), beEnd.getSV().origId(), beStart.getChrArm()));

                        writer.write(String.format(",%s,%d,%d,%d,%d,%d,%d,%s,%d",
                                pair.isAssembled(), pair.length(), cluster.getSynDelDupLength(),
                                pair.getNextSVDistance(), pair.getNextSVTraversedCount(),
                                pair.getDBLenFirst(), pair.getDBLenSecond(), pair.onArmOfOrigin(),
                                pair.getTraversedSVCount()));


                        writer.write(String.format(",%d,%d,%s,%s,%s",
                                beStart.position(), beEnd.position(),
                                beStart.getSV().getGeneInBreakend(beStart.usesStart()),
                                beEnd.getSV().getGeneInBreakend(beEnd.usesStart()), pair.getExonMatchData()));

                        writer.newLine();
                    }
                }
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing links to outputFile: {}", e.toString());
        }
    }


    public void close()
    {
        closeBufferedWriter(mSvFileWriter);
        closeBufferedWriter(mClusterFileWriter);
        closeBufferedWriter(mLinksFileWriter);
        mVisWriter.close();

        // log perf stats
        mPcPrep.logStats();
        mPcClusterAnalyse.logStats();
        mPcWrite.logStats();

        mAnalyser.logStats();
    }

}
