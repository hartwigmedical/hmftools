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
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_OTHER;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_MULTIPLE_DSBS;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_REMOTE_TI;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_SIMPLE_FOLDBACK;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_SINGLE;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.getArmClusterData;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_LOW_QUALITY;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SIMPLE_SV;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLOH;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.apache.commons.math3.optim.SimpleVectorValueChecker;
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
    private BufferedWriter mVisSvsFileWriter;
    private BufferedWriter mVisSegmentsFileWriter;

    private FragileSiteAnnotator mFragileSiteAnnotator;
    private LineElementAnnotator mLineElementAnnotator;
    private SvClusteringMethods mClusteringMethods;

    private PerformanceCounter mPerfCounter;
    private PerformanceCounter mPc1;
    private PerformanceCounter mPc2;
    private PerformanceCounter mPc3;
    private PerformanceCounter mPc4;
    private PerformanceCounter mPc5;

    private static final Logger LOGGER = LogManager.getLogger(SvSampleAnalyser.class);

    public SvSampleAnalyser(final SvaConfig config)
    {
        mConfig = config;
        mClusteringMethods = new SvClusteringMethods(mConfig.ProximityDistance);
        mAnalyser = new ClusterAnalyser(config, mClusteringMethods);

        mSvFileWriter = null;
        mLinksFileWriter = null;
        mClusterFileWriter = null;
        mVisSvsFileWriter = null;
        mVisSegmentsFileWriter = null;

        mFragileSiteAnnotator = new FragileSiteAnnotator();
        mFragileSiteAnnotator.loadFragileSitesFile(mConfig.FragileSiteFile);

        mLineElementAnnotator = new LineElementAnnotator();
        mLineElementAnnotator.loadLineElementsFile(mConfig.LineElementFile);

        mPerfCounter = new PerformanceCounter("Total");

        mPc1 = new PerformanceCounter("Annotate&Filter");
        mPc2 = new PerformanceCounter("ArmsStats");
        mPc3 = new PerformanceCounter("ClusterAndAnalyse");
        // mPc4 = new PerformanceCounter("Analyse");
        mPc5 = new PerformanceCounter("WriteCSV");

        mPerfCounter.start();

        clearState();
    }

    public final List<SvCluster> getClusters() { return mAnalyser.getClusters(); }
    public final Map<String, List<SvBreakend>> getChrBreakendMap() { return mClusteringMethods.getChrBreakendMap(); }
    public final Map<String, Double> getChrCopyNumberMap() { return mClusteringMethods.getChrCopyNumberMap(); }
    public void setSampleLohData(final Map<String, List<SvLOH>> data) { mClusteringMethods.setSampleLohData(data); }

    private void clearState()
    {
        mSampleId = "";
        mAllVariants = Lists.newArrayList();
    }

    public void loadFromDatabase(final String sampleId, final List<SvVarData> variants)
    {
        clearState();

        if (variants.isEmpty())
            return;

        mSampleId = sampleId;
        mAllVariants = Lists.newArrayList(variants);

        LOGGER.debug("loaded {} SVs", mAllVariants.size());
    }

    public void analyse()
    {
        if(mAllVariants.isEmpty())
            return;

        mPc1.start();
        annotateAndFilterVariants();
        mPc1.stop();

        LOGGER.debug("sample({}) clustering {} variants", mSampleId, mAllVariants.size());

        mPc2.start();
        mClusteringMethods.setChromosomalArmStats(mAllVariants);
        mClusteringMethods.populateChromosomeBreakendMap(mAllVariants);
        mClusteringMethods.annotateNearestSvData();
        LinkFinder.findDeletionBridges(mClusteringMethods.getChrBreakendMap());
        mClusteringMethods.setSimpleVariantLengths(mSampleId);
        mPc2.stop();

        mPc3.start();

        mAnalyser.setSampleData(mSampleId, mAllVariants);
        mAnalyser.clusterAndAnalyse();

        mPc3.stop();

        // logSampleClusterInfo();
    }

    public void writeOutput()
    {
        mPc5.start();

        if(!mConfig.OutputCsvPath.isEmpty())
        {
            writeClusterSVOutput();

            if(mConfig.hasMultipleSamples())
            {
                writeClusterLinkData();
                writeClusterData();
            }

            if(mConfig.WriteVisualisationData)
            {
                writeVisualSvData();
                writeVisualSegmentData();
            }
        }

        mPc5.stop();

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

                if(!outputFileName.endsWith("/"))
                    outputFileName += File.separator;

                if(mConfig.hasMultipleSamples())
                    outputFileName += "SVA_SVS.csv";
                else
                    outputFileName += mSampleId + ".csv";

                mSvFileWriter = createBufferedWriter(outputFileName, false);

                // definitional fields
                mSvFileWriter.write("SampleId,Id,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd");

                // position and copy number
                mSvFileWriter.write(",ArmStart,AdjAFStart,AdjCNStart,AdjCNChgStart,ArmEnd,AdjAFEnd,AdjCNEnd,AdjCNChgEnd,Ploidy");

                // cluster info
                mSvFileWriter.write(",ClusterId,SubClusterId,ClusterCount,ClusterReason");
                mSvFileWriter.write(",ClusterDesc,IsResolved,ResolvedType,Consistency,ArmCount");

                // SV info
                mSvFileWriter.write(",Homology,InexactHOStart,InexactHOEnd,InsertSeq,Imprecise,QualScore,RefContextStart,RefContextEnd");

                mSvFileWriter.write(",FSStart,FSEnd,LEStart,LEEnd,DupBEStart,DupBEEnd,ArmCountStart,ArmExpStart,ArmCountEnd,ArmExpEnd");

                // linked pair info
                mSvFileWriter.write(",LnkSvStart,LnkLenStart,LnkSvEnd,LnkLenEnd");

                // GRIDDS caller info
                mSvFileWriter.write(",AsmbStart,AsmbEnd,AsmbMatchStart,AsmbMatchEnd");

                // chain info
                mSvFileWriter.write(",ChainId,ChainCount,ChainIndex");

                // proximity info and other link info
                mSvFileWriter.write(",NearestLen,NearestType,DBLenStart,DBLenEnd,SynDelDupLen,SynDelDupTILen");

                // proximity info and other link info
                mSvFileWriter.write(",FoldbackLnkStart,FoldbackLenStart,FoldbackLinkInfoStart,FoldbackLnkEnd,FoldbackLenEnd,FoldbackLinkInfoEnd");

                // gene info
                mSvFileWriter.write(",GeneStart,GeneEnd");

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

                int clusterSvCount = cluster.getUniqueSvCount();

                SvCluster subCluster = cluster;
                if(cluster.hasSubClusters())
                {
                    for(final SvCluster sc : cluster.getSubClusters())
                    {
                        if(sc.getSVs().contains(var))
                        {
                            subCluster = sc;
                            break;
                        }
                    }
                }

                final StructuralVariantData dbData = var.getSvData();

                ++svCount;

                writer.write(
                        String.format("%s,%s,%s,%s,%d,%d,%s,%d,%d",
                                mSampleId, var.id(), var.typeStr(),
                                var.chromosome(true), var.position(true), var.orientation(true),
                                var.chromosome(false), var.position(false), var.orientation(false)));

                writer.write(
                        String.format(",%s,%.2f,%.2f,%.2f,%s,%.2f,%.2f,%.2f,%.2f",
                                var.arm(true), dbData.adjustedStartAF(), dbData.adjustedStartCopyNumber(), dbData.adjustedStartCopyNumberChange(),
                                var.arm(false), dbData.adjustedEndAF(), dbData.adjustedEndCopyNumber(), dbData.adjustedEndCopyNumberChange(), dbData.ploidy()));

                writer.write(
                        String.format(",%d,%d,%d,%s",
                                cluster.id(), subCluster.id(), clusterSvCount, var.getClusterReason()));

                writer.write(
                        String.format(",%s,%s,%s,%d,%d",
                                cluster.getDesc(), cluster.isResolved(), cluster.getResolvedType(),
                                cluster.getConsistencyCount(), cluster.getArmCount()));

                int dbLenStart = var.getDBLink(true) != null ? var.getDBLink(true).length() : NO_DB_MARKER;
                int dbLenEnd = var.getDBLink(false) != null ? var.getDBLink(false).length() : NO_DB_MARKER;

                writer.write(
                        String.format(",%s,%d,%d,%s,%s,%.0f,%s,%s",
                                dbData.insertSequence().isEmpty() && var.type() != INS ? dbData.homology() : "",
                                dbData.inexactHomologyOffsetStart(), dbData.inexactHomologyOffsetEnd(),
                                dbData.insertSequence(), dbData.imprecise(), dbData.qualityScore(),
                                dbData.startRefContext(), dbData.endRefContext()));

                writer.write(
                        String.format(",%s,%s,%s,%s,%s,%s,%s",
                                var.isFragileSite(true), var.isFragileSite(false),
                                var.getLineElement(true), var.getLineElement(false),
                                var.isDupBreakend(true), var.isDupBreakend(false),
                                mClusteringMethods.getChrArmData(var)));

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
                        startLinkStr, endLinkStr, var.getAssemblyData(true), var.getAssemblyData(false), assemblyMatchStart, assemblyMatchEnd));

                // chain info
                final SvChain chain = cluster.findChain(var);
                String chainStr = "";

                if(chain != null)
                {
                    chainStr = String.format(",%d,%d,%s", chain.id(), chain.getUniqueSvCount(), chain.getSvIndices(var));
                }
                else
                {
                    chainStr = String.format(",%d,0,", cluster.getChainId(var));
                }

                writer.write(chainStr);

                writer.write(String.format(",%d,%s,%d,%d,%d,%d",
                        var.getNearestSvDistance(), var.getNearestSvRelation(), dbLenStart, dbLenEnd,
                        cluster.getSynDelDupLength(), cluster.getSynDelDupTILength()));

                writer.write(String.format(",%s,%d,%s,%s,%d,%s",
                        var.getFoldbackLink(true), var.getFoldbackLen(true), var.getFoldbackLinkInfo(true),
                        var.getFoldbackLink(false), var.getFoldbackLen(false), var.getFoldbackLinkInfo(false)));

                writer.write(String.format(",%s,%s",
                        var.getGeneInBreakend(true), var.getGeneInBreakend(false)));

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

                if(!outputFileName.endsWith("/"))
                    outputFileName += File.separator;

                outputFileName += "SVA_CLUSTERS.csv";

                mClusterFileWriter = createBufferedWriter(outputFileName, false);

                mClusterFileWriter.write("SampleId,ClusterId,ClusterDesc,ClusterCount,ResolvedType,FullyChained,ChainCount");
                mClusterFileWriter.write(",DelCount,DupCount,InsCount,InvCount,BndCount,SglCount");
                mClusterFileWriter.write(",Consistency,ArmCount,OriginArms,FragmentArms,IsLINE,HasReplicated,Foldbacks,DSBs");
                mClusterFileWriter.write(",TotalLinks,AssemblyLinks,LongDelDups,UnlinkedRemotes,ShortTIRemotes,MinCopyNumber,MaxCopyNumber");
                mClusterFileWriter.write(",SynDelDupLen,SynDelDupAvgTILen,Annotations,ChainInfo");
                mClusterFileWriter.write(",ArmClusterCount,AcSoloSv,AcRemoteTI,AcDsb,AcMultipleDsb,AcSimpleFoldback,AcComplexFoldback,AcComplexOther");
                mClusterFileWriter.newLine();
            }

            BufferedWriter writer = mClusterFileWriter;

            for(final SvCluster cluster : getClusters())
            {
                int clusterSvCount = cluster.getUniqueSvCount();

                writer.write(
                        String.format("%s,%d,%s,%d,%s,%s,%d",
                                mSampleId, cluster.id(), cluster.getDesc(), clusterSvCount, cluster.getResolvedType(),
                                cluster.isFullyChained(), cluster.getChains().size()));

                writer.write(
                        String.format(",%d,%d,%d,%d,%d,%d",
                                cluster.getTypeCount(DEL), cluster.getTypeCount(DUP), cluster.getTypeCount(INS),
                                cluster.getTypeCount(INV), cluster.getTypeCount(BND), cluster.getTypeCount(SGL)));

                writer.write(
                        String.format(",%d,%d,%d,%d,%s,%s,%d,%d",
                                cluster.getConsistencyCount(), cluster.getArmCount(), cluster.getOriginArms(), cluster.getFragmentArms(),
                                cluster.hasLinkingLineElements(), cluster.hasReplicatedSVs(),
                                cluster.getFoldbacks().size(), cluster.getClusterDBCount()));

                final String chainInfo = cluster.getChains().stream()
                        .filter(x -> !x.getDetails().isEmpty())
                        .map(SvChain::getDetails)
                        .collect (Collectors.joining (";"));

                writer.write(
                        String.format(",%d,%d,%d,%d,%d,%d,%d",
                                cluster.getLinkedPairs().size(), cluster.getAssemblyLinkedPairs().size(), cluster.getLongDelDups().size(),
                                cluster.getUnlinkedRemoteSVs().size(), cluster.getShortTIRemoteSVs().size(),
                                cluster.getMinCNChange(), cluster.getMaxCNChange()));

                writer.write(
                        String.format(",%d,%d,%s,%s",
                                cluster.getSynDelDupLength(), cluster.getSynDelDupTILength(), cluster.getAnnotations(), chainInfo));

                final int[] armClusterData = getArmClusterData(cluster);
                writer.write(
                        String.format(",%d,%d,%d,%d,%d,%d,%d,%d",
                                cluster.getArmClusters().size(), armClusterData[ARM_CL_SINGLE], armClusterData[ARM_CL_REMOTE_TI],
                                armClusterData[ARM_CL_DSB], armClusterData[ARM_CL_MULTIPLE_DSBS], armClusterData[ARM_CL_SIMPLE_FOLDBACK],
                                armClusterData[ARM_CL_COMPLEX_FOLDBACK], armClusterData[ARM_CL_COMPLEX_OTHER]));


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

                if(!outputFileName.endsWith("/"))
                    outputFileName += File.separator;

                outputFileName += "SVA_LINKS.csv";

                mLinksFileWriter = createBufferedWriter(outputFileName, false);

                mLinksFileWriter.write("SampleId,ClusterId,ClusterDesc,ClusterCount,ResolvedType,IsLINE,FullyChained");
                mLinksFileWriter.write(",ChainId,ChainCount,ChainConsistent,Id1,Id2,ChrArm,IsAssembled,TILength,SynDelDupLen");
                mLinksFileWriter.write(",NextSVDistance,NextSVTraversedCount,DBLenStart,DBLenEnd,OnArmOfOrigin,CopyNumberGain,TraversedSVCount");
                mLinksFileWriter.write(",PosStart,PosEnd,GeneStart,GeneEnd,ExonMatch");
                mLinksFileWriter.newLine();
            }

            BufferedWriter writer = mLinksFileWriter;

            for(final SvCluster cluster : getClusters())
            {
                int clusterSvCount = cluster.getUniqueSvCount();

                isSpecificCluster(cluster);

                List<SvChain> chains = cluster.getChains();

                for (final SvChain chain : chains)
                {
                    int chainSvCount = chain.getSvCount();
                    boolean chainConsistent = chain.isConsistent();

                    for (final SvLinkedPair pair : chain.getLinkedPairs())
                    {
                        if(pair.linkType() != LINK_TYPE_TI)
                            continue;

                        writer.write(
                                String.format("%s,%d,%s,%d,%s,%s,%s",
                                        mSampleId, cluster.id(), cluster.getDesc(), clusterSvCount, cluster.getResolvedType(),
                                        cluster.hasLinkingLineElements(), cluster.isFullyChained()));

                        final SvBreakend beStart = pair.getBreakend(true);
                        final SvBreakend beEnd= pair.getBreakend(false);

                        writer.write(
                                String.format(",%d,%d,%s,%s,%s,%s",
                                        chain.id(), chainSvCount, chainConsistent,
                                        beStart.getSV().origId(), beEnd.getSV().origId(), beStart.getChrArm()));

                        writer.write(
                                String.format(",%s,%d,%d,%d,%d,%d,%d,%s,%s,%d",
                                        pair.isAssembled(), pair.length(), cluster.getSynDelDupLength(),
                                        pair.getNextSVDistance(), pair.getNextSVTraversedCount(),
                                        pair.getDBLenFirst(), pair.getDBLenSecond(), pair.onArmOfOrigin(),
                                        pair.hasCopyNumberGain(), pair.getTraversedSVCount()));


                        writer.write(
                                String.format(",%d,%d,%s,%s,%s",
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

    private void writeVisualSvData()
    {
        try
        {
            if (mVisSvsFileWriter == null)
            {
                String outputFileName = mConfig.OutputCsvPath;

                if (!outputFileName.endsWith("/"))
                    outputFileName += File.separator;

                outputFileName += "SVA_VIS_SVS.csv";

                mVisSvsFileWriter = createBufferedWriter(outputFileName, false);

                // definitional fields
                mVisSvsFileWriter.write("SampleId,ClusterId,ChainId,SvId,Type,ResolvedType");
                mVisSvsFileWriter.write(",ChrStart,PosStart,OrientStart,InfoStart,ChrEnd,PosEnd,OrientEnd,InfoEnd,TraverseCount");

                mVisSvsFileWriter.newLine();
            }

            BufferedWriter writer = mVisSvsFileWriter;

            for(final SvVarData var : mAllVariants)
            {
                final SvChain chain = var.getCluster().findChain(var);
                int chainId = chain != null ? chain.id() : var.getCluster().getChainId(var);

                writer.write(
                        String.format("%s,%d,%d,%s,%s,%s",
                                mSampleId, var.getCluster().id(), chainId, var.id(),
                                var.type(), var.getCluster().getResolvedType()));

                for(int be = SVI_START; be <= SVI_END; ++be)
                {
                    boolean isStart = isStart(be);

                    if(!isStart && var.isNullBreakend())
                    {
                        writer.write(",-1,0,0,NULL");
                        continue;
                    }

                    final SvBreakend breakend = var.getBreakend(isStart);

                    writer.write(
                            String.format(",%s,%d,%d,%s",
                                    breakend.chromosome(), breakend.position(), breakend.orientation(),
                                    breakend.getSV().getFoldbackLink(isStart).isEmpty() ? "NORMAL" : "FOLDBACK"));
                }

                int repeatCount = chain != null ? max(var.getReplicatedCount(), 1) : 1;
                writer.write(String.format(",%d", repeatCount));

                writer.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to visual SVs file: {}", e.toString());
        }
    }

    private void writeVisualSegmentData()
    {
        try
        {
            if (mVisSegmentsFileWriter == null)
            {
                String outputFileName = mConfig.OutputCsvPath;

                if (!outputFileName.endsWith("/"))
                    outputFileName += File.separator;

                outputFileName += "SVA_VIS_SEGMENTS.csv";

                mVisSegmentsFileWriter = createBufferedWriter(outputFileName, false);
                mVisSegmentsFileWriter.write("SampleId,ClusterId,ChainId,Chr,PosStart,PosEnd,TraverseCount");
                mVisSegmentsFileWriter.newLine();
            }

            BufferedWriter writer = mVisSegmentsFileWriter;

            for(final SvCluster cluster : getClusters())
            {
                if(cluster.isResolved()
                && (cluster.getResolvedType() == RESOLVED_TYPE_SIMPLE_SV || cluster.getResolvedType() == RESOLVED_TYPE_LOW_QUALITY))
                {
                    continue;
                }

                for (final SvChain chain : cluster.getChains())
                {
                    List<SvLinkedPair> uniquePairs = Lists.newArrayList();

                    // log the start of the chain
                    SvBreakend breakend = chain.getFirstSV().getBreakend(chain.firstLinkOpenOnStart());
                    boolean startsOnEnd = chain.getFirstSV().equals(chain.getLastSV(), true);

                    if(breakend != null)
                    {
                        writer.write(String.format("%s,%d,%d,%s,%s,%s,%d",
                                mSampleId, cluster.id(), chain.id(), breakend.chromosome(), getPositionValue(breakend, true),
                                getPositionValue(breakend, false), startsOnEnd ? 2 : 1));

                        writer.newLine();
                    }

                    for (final SvLinkedPair pair : chain.getLinkedPairs())
                    {
                        boolean isRepeat = false;

                        for (final SvLinkedPair existingPair : uniquePairs)
                        {
                            if(pair.matches(existingPair, true))
                            {
                                isRepeat = true;
                                break;
                            }
                        }

                        if(isRepeat)
                            continue;

                        uniquePairs.add(pair);

                        writer.write(String.format("%s,%d,%d",
                                mSampleId, cluster.id(), chain.id()));

                        int pairRepeatCount = 0;

                        for (final SvLinkedPair existingPair : chain.getLinkedPairs())
                        {
                            if(pair.matches(existingPair, true))
                                ++pairRepeatCount;
                        }

                        final SvBreakend beStart = pair.getBreakend(true);
                        final SvBreakend beEnd= pair.getBreakend(false);

                        writer.write(String.format(",%s,%d,%d,%d",
                                beStart.chromosome(), beStart.position(), beEnd.position(), pairRepeatCount));

                        writer.newLine();
                    }

                    // log the end of the chain out to centromere or telomere
                    // log the start of the chain
                    breakend = chain.getLastSV().getBreakend(chain.lastLinkOpenOnStart());

                    if(breakend != null && !startsOnEnd)
                    {
                        writer.write(String.format("%s,%d,%d,%s,%s,%s,%d",
                                mSampleId, cluster.id(), chain.id(), breakend.chromosome(), getPositionValue(breakend, true),
                                getPositionValue(breakend, false), 1));

                        writer.newLine();
                    }
                }

                // finally write out all unchained SVs
                for(final SvVarData var : cluster.getUnlinkedSVs())
                {
                    if(var.isReplicatedSv())
                        continue;

                    int chainId = cluster.getChainId(var);

                    for(int be = SVI_START; be <= SVI_END; ++be)
                    {
                        final SvBreakend breakend = var.getBreakend(isStart(be));

                        if(breakend == null)
                            continue;

                        writer.write(String.format("%s,%d,%d",
                                mSampleId, cluster.id(), chainId));

                        writer.write(String.format(",%s,%s,%s,%d",
                                breakend.chromosome(), getPositionValue(breakend, true),
                                getPositionValue(breakend, false), 1));

                        writer.newLine();
                    }
                }
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to visual segments file: {}", e.toString());
        }
    }

    private static final String getPositionValue(final SvBreakend breakend, boolean isStart)
    {
        if(breakend.orientation() == 1 && breakend.arm().equals(CHROMOSOME_ARM_P))
        {
            return isStart ? "T" : breakend.position().toString();
        }
        else if(breakend.orientation() == -1 && breakend.arm().equals(CHROMOSOME_ARM_P))
        {
            return isStart ? breakend.position().toString() : "C";
        }
        else if(breakend.orientation() == -1 && breakend.arm().equals(CHROMOSOME_ARM_Q))
        {
            return isStart ? breakend.position().toString() : "T";
        }
        else
        {
            return isStart ? "C" : breakend.position().toString();
        }
    }

    private void logSampleClusterInfo()
    {
        int simpleClusterCount = 0;
        int complexClusterCount = 0;
        Map<String, Integer> armSimpleClusterCount = Maps.newHashMap();
        Map<String, Integer> armComplexClusterCount = Maps.newHashMap();

        for(final SvCluster cluster : mAnalyser.getClusters())
        {
            Map<String, Integer> targetMap;

            if(cluster.getResolvedType() == RESOLVED_TYPE_LOW_QUALITY)
                continue;

            if(cluster.getResolvedType() == RESOLVED_TYPE_SIMPLE_SV || cluster.isSyntheticSimpleType(true))
            {
                ++simpleClusterCount;
                targetMap = armSimpleClusterCount;
            }
            else
            {
                ++complexClusterCount;
                targetMap = armComplexClusterCount;
            }

            for(final SvArmGroup armGroup : cluster.getArmGroups())
            {
                if(targetMap.containsKey(armGroup))
                {
                    targetMap.put(armGroup.id(), targetMap.get(armGroup.id()) + 1);
                }
                else
                {
                    targetMap.put(armGroup.id(), 1);
                }
            }
        }

        int armsWithExcessComplexClusters = 0;

        for(final Map.Entry<String, Integer> entry : armComplexClusterCount.entrySet())
        {
            final String chrArm = entry.getKey();
            int complexCount = entry.getValue();
            Integer simpleCount = armSimpleClusterCount.containsKey(chrArm) ? armSimpleClusterCount.get(chrArm) : 0;

            if(simpleCount > 0 && complexCount > simpleCount)
            {
                LOGGER.debug("chrArm({}) clusters simple({}) vs complex({})", chrArm, simpleCount, complexCount);
                ++armsWithExcessComplexClusters;
            }
        }

        if(complexClusterCount > simpleClusterCount || armsWithExcessComplexClusters >= 2)
        {
            LOGGER.info("sample({}) clusters total({}) simple({}) complex({}) excessComplexArms({})",
                    mSampleId, mAnalyser.getClusters().size(), simpleClusterCount, complexClusterCount, armsWithExcessComplexClusters);
        }
    }

    public void close()
    {
        closeBufferedWriter(mSvFileWriter);
        closeBufferedWriter(mClusterFileWriter);
        closeBufferedWriter(mLinksFileWriter);
        closeBufferedWriter(mVisSegmentsFileWriter);
        closeBufferedWriter(mVisSvsFileWriter);

        // log perf stats
        mPerfCounter.stop();
        mPerfCounter.logStats(false);
        mPc1.logStats(false);
        mPc2.logStats(false);
        mPc3.logStats(false);
        // mPc4.logStats(false);
        mPc5.logStats(false);

        mAnalyser.logStats();
    }

}
