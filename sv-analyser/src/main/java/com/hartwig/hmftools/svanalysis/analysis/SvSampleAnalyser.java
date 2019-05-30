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
import static com.hartwig.hmftools.svanalysis.analysis.SvClassification.getSuperType;
import static com.hartwig.hmftools.svanalysis.analysis.SvClassification.isSimpleType;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_LINE;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_OTHER;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_FOLDBACK;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_FOLDBACK_DSB;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_SIMPLE_DUP;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_TI_ONLY;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_ISOLATED_BE;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.getArmClusterData;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_CHAIN_ENDS_AWAY;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_CHAIN_ENDS_FACE;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_EXT_SHORT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_EXT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_EXT_TI_CN_GAIN;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_DB;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_INT_SHORT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_OVERLAPPING_TI;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_SHORT_DB;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_INT_TI;
import static com.hartwig.hmftools.svanalysis.types.SvChain.CM_INT_TI_CN_GAIN;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SE_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvaConstants.NO_DB_MARKER;
import static com.hartwig.hmftools.svanalysis.types.SvaConstants.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PSEUDO_GENE_DATA_EXON_LENGTH;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PSEUDO_GENE_DATA_EXON_MAX;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PSEUDO_GENE_DATA_EXON_RANK;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PSEUDO_GENE_DATA_TRANS_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.ReplicationOriginAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.VisualiserWriter;
import com.hartwig.hmftools.svanalysis.types.SvArmCluster;
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

    private boolean mIsValid;

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

        mIsValid = true;

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
    public boolean inValidState() { return mIsValid; }
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

            for(int be = SE_START; be <= SE_END; ++be)
            {
                if(var.isNullBreakend() && be == SE_END)
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

        mIsValid = mAnalyser.clusterAndAnalyse();

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

            var.setFragileSites(
                    mFragileSiteAnnotator.isFragileSite(var, true),
                    mFragileSiteAnnotator.isFragileSite(var, false));

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

    public void annotateWithGeneData(SvGeneTranscriptCollection geneTransCache)
    {
        checkPseudoGeneAnnotations(geneTransCache);
    }

    private void checkPseudoGeneAnnotations(SvGeneTranscriptCollection geneTransCache)
    {
        for(final SvCluster cluster : getClusters())
        {
            // isSpecificCluster(cluster);

            GeneAnnotation pseudoGene = null;
            String transcriptId = "";

            for(final SvLinkedPair pair : cluster.getLinkedPairs())
            {
                if(pair.length() > SHORT_TI_LENGTH * 8)
                    continue;

                final SvBreakend lower = pair.getBreakend(true);
                final SvBreakend upper = pair.getBreakend(false);

                // for any TI falling within the same gene, check for an exon boundary match
                if(lower.getSV().getGenesList(lower.usesStart()).isEmpty() || upper.getSV().getGenesList(upper.usesStart()).isEmpty())
                    continue;

                for(final GeneAnnotation gene1 : lower.getSV().getGenesList(lower.usesStart()))
                {
                    for(final GeneAnnotation gene2 : upper.getSV().getGenesList(upper.usesStart()))
                    {
                        if(!gene1.GeneName.equals(gene2.GeneName))
                            continue;

                        final String exonData[] = geneTransCache.getExonDetailsForPosition(gene1, lower.position(), upper.position());

                        if(exonData[PSEUDO_GENE_DATA_TRANS_ID] != null)
                        {
                            pseudoGene = gene1;
                            transcriptId = exonData[PSEUDO_GENE_DATA_TRANS_ID];

                            String exonMatchData = String.format("%s;%s;%s;%s",
                                    transcriptId, exonData[PSEUDO_GENE_DATA_EXON_RANK],
                                    exonData[PSEUDO_GENE_DATA_EXON_MAX], exonData[PSEUDO_GENE_DATA_EXON_LENGTH]);


                            pair.setExonMatchData(exonMatchData);
                        }
                    }
                }
            }

            if(pseudoGene != null)
            {
                mVisWriter.addGeneExonData(cluster.id(), pseudoGene.StableId, pseudoGene.GeneName,
                        transcriptId, pseudoGene.chromosome(), "PSEUDO");
            }
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

                mSvFileWriter.write(",FSStart,FSEnd,LEStart,LEEnd");

                // linked pair info
                mSvFileWriter.write(",LnkSvStart,LnkLenStart,LnkSvEnd,LnkLenEnd");
                mSvFileWriter.write(",AsmbStart,AsmbEnd,AsmbMatchStart,AsmbMatchEnd");

                // chain info
                mSvFileWriter.write(",ChainId,ChainCount,ChainIndex");

                // proximity info and other link info
                mSvFileWriter.write(",NearestLen,NearestType,DBLenStart,DBLenEnd");

                // proximity info and other link info
                mSvFileWriter.write(",FoldbackLnkStart,FoldbackLenStart,FoldbackInfoStart,FoldbackLnkEnd,FoldbackLenEnd,FoldbackInfoEnd");

                // local topology from arm cluster
                mSvFileWriter.write(",LocTopIdStart,LocTopTypeStart,LocTopTIStart,LocTopIdEnd,LocTopTypeEnd,LocTopTIEnd");

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
                for(int be = SE_START; be <= SE_END; ++be)
                {
                    boolean isStart = isStart(be);
                    final SvLinkedPair link = var.getLinkedPair(isStart);
                    if(link != null)
                    {
                        writer.write(String.format(",%s,%d",
                            link.first().equals(var, true) ? link.second().origId() : link.first().origId(), link.length()));
                    }
                    else
                    {
                        writer.write(",,-1");
                    }
                }

                // assembly info
                writer.write(String.format(",%s,%s,%s,%s",
                        var.getAssemblyData(true), var.getAssemblyData(false),
                        var.getAssemblyMatchType(true), var.getAssemblyMatchType(false)));

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

                writer.write(String.format(",%d,%s,%d,%d",
                        var.getNearestSvDistance(), var.getNearestSvRelation(), dbLenStart, dbLenEnd));

                writer.write(String.format(",%s,%d,%s,%s,%d,%s",
                        var.getFoldbackLink(true), var.getFoldbackLength(true), var.getFoldbackInfo(true),
                        var.getFoldbackLink(false), var.getFoldbackLength(false), var.getFoldbackInfo(false)));

                for(int be = SE_START; be <= SE_END; ++be)
                {
                    final SvArmCluster armCluster = (be == SE_START || !var.isNullBreakend()) ?
                            cluster.findArmCluster(var.getBreakend(isStart(be))) : null;

                    if(armCluster != null)
                        writer.write(String.format(",%d,%s,%d", armCluster.id(), armCluster.getTypeStr(), armCluster.getTICount()));
                    else
                        writer.write(",-1,,0");
                }

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

                mClusterFileWriter.write("SampleId,ClusterId,ClusterDesc,ClusterCount,SuperType,ResolvedType,Synthetic,Subclonal,FullyChained,ChainCount");
                mClusterFileWriter.write(",DelCount,DupCount,InsCount,InvCount,BndCount,SglCount,InfCount");
                mClusterFileWriter.write(",ClusterReasons,Consistency,ArmCount,IsLINE,Replication,MinPloidy,MaxPloidy,Foldbacks");
                mClusterFileWriter.write(",TotalTIs,AssemblyTIs,ShortTIs,IntTIs,ExtTIs,IntShortTIs,ExtShortTIs,IntTIsCnGain,ExtTIsCnGain,OverlapTIs");
                mClusterFileWriter.write(",DSBs,ShortDSBs,ChainEndsFace,ChainEndsAway,SyntheticLen,SyntheticTILen");
                mClusterFileWriter.write(",OriginArms,FragmentArms,UnchainedSVs,Annotations,AlleleValidPerc");
                mClusterFileWriter.write(",ArmClusterCount,AcTotalTIs,AcIsolatedBE,AcTIOnly,AcDsb,AcSimpleDup");
                mClusterFileWriter.write(",AcSingleFb,AcFbDsb,AcComplexFb,AcComplexLine,AcComplexOther");
                mClusterFileWriter.newLine();
            }

            BufferedWriter writer = mClusterFileWriter;

            for(final SvCluster cluster : getClusters())
            {
                int clusterSvCount = cluster.getSvCount();

                // isSpecificCluster(cluster);

                String resolvedType = cluster.getResolvedType();

                // TEMP - keep synthetic formed from SGLs separate

                int inferredCount = cluster.getInferredTypeCount();
                int sglCount = cluster.getTypeCount(SGL);

                if((inferredCount > 0 || sglCount> 0) && isSimpleType(resolvedType))
                {
                    resolvedType = "SGL_PAIR_" + resolvedType;
                }

                writer.write(String.format("%s,%d,%s,%d,%s,%s,%s,%s,%s,%d",
                        mSampleId, cluster.id(), cluster.getDesc(), clusterSvCount,
                        getSuperType(cluster), resolvedType, cluster.isSyntheticType(),
                        cluster.isSubclonal(), cluster.isFullyChained(false), cluster.getChains().size()));

                writer.write(String.format(",%d,%d,%d,%d,%d,%d,%d",
                        cluster.getTypeCount(DEL), cluster.getTypeCount(DUP), cluster.getTypeCount(INS),
                        cluster.getTypeCount(INV), cluster.getTypeCount(BND), sglCount - inferredCount, inferredCount));

                double foldbackCount = 0;

                for(final SvVarData var : cluster.getFoldbacks())
                {
                    // avoid double-count chained foldbacks
                    if(var.getFoldbackBreakend(true) != null)
                        foldbackCount += 0.5;
                    if(var.getFoldbackBreakend(false) != null)
                        foldbackCount += 0.5;
                }

                writer.write(String.format(",%s,%d,%d,%s,%s,%.0f,%.0f,%.0f",
                        cluster.getClusteringReasons(), cluster.getConsistencyCount(), cluster.getArmCount(),
                        cluster.hasLinkingLineElements(), cluster.hasReplicatedSVs(), cluster.getMinPloidy(), cluster.getMaxPloidy(),
                        foldbackCount));

                long shortTIs = cluster.getLinkedPairs().stream().filter(x -> x.length() <= SHORT_TI_LENGTH).count();

                writer.write(String.format(",%d,%d,%d",
                        cluster.getLinkedPairs().size(), cluster.getAssemblyLinkedPairs().size(), shortTIs));

                int[] chainData = cluster.getLinkMetrics();

                writer.write(String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
                        chainData[CM_INT_TI], chainData[CM_EXT_TI], chainData[CM_INT_SHORT_TI], chainData[CM_EXT_SHORT_TI],
                        chainData[CM_INT_TI_CN_GAIN], chainData[CM_EXT_TI_CN_GAIN], chainData[CM_OVERLAPPING_TI],
                        chainData[CM_DB], chainData[CM_SHORT_DB], chainData[CM_CHAIN_ENDS_FACE], chainData[CM_CHAIN_ENDS_AWAY]));

                writer.write(String.format(",%d,%d,%d,%d,%d,%s,%.2f",
                        cluster.getSyntheticLength(), cluster.getSyntheticTILength(), cluster.getOriginArms(), cluster.getFragmentArms(),
                        cluster.getUnlinkedSVs().size(), cluster.getAnnotations(), cluster.getValidAllelePloidySegmentPerc()));

                final int[] armClusterData = getArmClusterData(cluster);
                long armClusterTIs = cluster.getArmClusters().stream().mapToInt(x -> x.getTICount()).sum();

                writer.write(String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
                        cluster.getArmClusters().size(), armClusterTIs, armClusterData[ARM_CL_ISOLATED_BE], armClusterData[ARM_CL_TI_ONLY],
                        armClusterData[ARM_CL_DSB], armClusterData[ARM_CL_SIMPLE_DUP], armClusterData[ARM_CL_FOLDBACK],
                        armClusterData[ARM_CL_FOLDBACK_DSB], armClusterData[ARM_CL_COMPLEX_FOLDBACK], armClusterData[ARM_CL_COMPLEX_LINE],
                        armClusterData[ARM_CL_COMPLEX_OTHER]));

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

                mLinksFileWriter.write("SampleId,ClusterId,ClusterDesc,ClusterCount,ResolvedType,IsLINE");
                mLinksFileWriter.write(",ChainId,ChainCount,ChainConsistent,Id1,Id2,ChrArm,IsAssembled,TILength");
                mLinksFileWriter.write(",NextSvDist,NextClusteredSvDist,TraversedSVCount,DBLenStart,DBLenEnd,OnArmOfOrigin");
                mLinksFileWriter.write(",LocationType,OverlapCount,CopyNumberGain");
                mLinksFileWriter.write(",PosStart,PosEnd,LocTopTypeStart,LocTopTypeEnd,GeneStart,GeneEnd,ExonMatch");
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
                        if(uniquePairs.stream().anyMatch(x -> x.matches(pair)))
                            continue;

                        uniquePairs.add(pair);

                        writer.write(String.format("%s,%d,%s,%d,%s,%s",
                            mSampleId, cluster.id(), cluster.getDesc(), clusterSvCount, cluster.getResolvedType(),
                            cluster.hasLinkingLineElements()));

                        final SvBreakend beStart = pair.getBreakend(true);
                        final SvBreakend beEnd = pair.getBreakend(false);

                        writer.write(String.format(",%d,%d,%s,%s,%s,%s",
                                chain.id(), chainSvCount, chainConsistent,
                                beStart.getSV().origId(), beEnd.getSV().origId(), beStart.getChrArm()));

                        writer.write(String.format(",%s,%d,%d,%d,%d,%d,%d,%s,%s,%d,%s",
                                pair.isAssembled(), pair.length(),
                                pair.getNextSvDistance(), pair.getNextClusteredSvDistance(), pair.getTraversedSVCount(),
                                pair.getDBLenFirst(), pair.getDBLenSecond(), pair.onArmOfOrigin(),
                                pair.locationType(), pair.overlapCount(), pair.hasCopyNumberGain()));

                        SvArmCluster acStart = cluster.findArmCluster(beStart);
                        SvArmCluster acEnd = cluster.findArmCluster(beEnd);

                        writer.write(String.format(",%d,%d,%s,%s,%s,%s,%s",
                                beStart.position(), beEnd.position(),
                                acStart != null ? acStart.getTypeStr() : "", acEnd != null ? acEnd.getTypeStr() : "",
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
