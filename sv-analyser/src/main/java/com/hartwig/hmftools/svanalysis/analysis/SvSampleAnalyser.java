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
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_COMPLEX_OTHER;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_MULTIPLE_DSBS;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_REMOTE_TI;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_SIMPLE_FOLDBACK;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.ARM_CL_SINGLE;
import static com.hartwig.hmftools.svanalysis.types.SvArmCluster.getArmClusterData;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.ReplicationOriginAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.VisualiserWriter;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLOH;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvaConfig;

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
        mVisWriter = new VisualiserWriter(config.OutputCsvPath, config.WriteVisualisationData);

        mAllVariants = Lists.newArrayList();

        mSvFileWriter = null;
        mLinksFileWriter = null;
        mClusterFileWriter = null;

        mFragileSiteAnnotator = new FragileSiteAnnotator();
        mFragileSiteAnnotator.loadFragileSitesFile(mConfig.FragileSiteFile);

        mLineElementAnnotator = new LineElementAnnotator();
        mLineElementAnnotator.loadLineElementsFile(mConfig.LineElementFile);

        mReplicationOriginAnnotator = new ReplicationOriginAnnotator();
        mReplicationOriginAnnotator.loadReplicationOrigins(mConfig.ReplicationOriginsFile);

        mPerfCounter = new PerformanceCounter("Total");

        // mPc1 = new PerformanceCounter("Annotate&Filter");
        mPc2 = new PerformanceCounter("Preparation");
        mPc3 = new PerformanceCounter("ClusterAndAnalyse");
        // mPc4 = new PerformanceCounter("Analyse");
        mPc5 = new PerformanceCounter("WriteCSV");

        mPerfCounter.start();

        clearState();
    }

    public final List<SvVarData> getVariants() { return mAllVariants; }
    public final List<SvCluster> getClusters() { return mAnalyser.getClusters(); }
    public final Map<String, List<SvBreakend>> getChrBreakendMap() { return mClusteringMethods.getChrBreakendMap(); }
    public void setSampleLohData(final Map<String, List<SvLOH>> data) { mClusteringMethods.setSampleLohData(data); }
    public void setChrCopyNumberMap(final Map<String, double[]> data) { mClusteringMethods.setChrCopyNumberMap(data); }
    public final VisualiserWriter getVisWriter() { return mVisWriter; }

    private void clearState()
    {
        mClusteringMethods.clearLOHBreakendData(mSampleId);
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

        LOGGER.debug("loaded {} SVs", mAllVariants.size());
    }

    public void analyse()
    {
        if(mAllVariants.isEmpty())
            return;

        LOGGER.debug("sample({}) clustering {} variants", mSampleId, mAllVariants.size());

        mPc2.start();
        annotateAndFilterVariants();
        mClusteringMethods.populateChromosomeBreakendMap(mAllVariants);
        mClusteringMethods.annotateNearestSvData();
        LinkFinder.findDeletionBridges(mClusteringMethods.getChrBreakendMap());
        mClusteringMethods.setSimpleVariantLengths(mSampleId);
        mReplicationOriginAnnotator.setReplicationOrigins(mClusteringMethods.getChrBreakendMap());
        mPc2.stop();

        mPc3.start();

        mAnalyser.setSampleData(mSampleId, mAllVariants);
        mAnalyser.clusterAndAnalyse();

        mPc3.stop();
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

            mVisWriter.writeOutput(mAnalyser.getClusters(), mAllVariants);
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

                if(mConfig.hasMultipleSamples())
                    outputFileName += "SVA_SVS.csv";
                else
                    outputFileName += mSampleId + "_SVA.csv";

                mSvFileWriter = createBufferedWriter(outputFileName, false);

                // definitional fields
                mSvFileWriter.write("SampleId,Id,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd");

                // position and copy number
                mSvFileWriter.write(",ArmStart,AdjAFStart,AdjCNStart,AdjCNChgStart,ArmEnd,AdjAFEnd,AdjCNEnd,AdjCNChgEnd,Ploidy");

                // cluster info
                mSvFileWriter.write(",ClusterId,SubClusterId,ClusterCount,ClusterReason");
                mSvFileWriter.write(",ClusterDesc,IsResolved,ResolvedType,Consistency,ArmCount");

                // SV info
                mSvFileWriter.write(",Homology,InexactHOStart,InexactHOEnd,InsertSeq,Imprecise,QualScore,RefContextStart,RefContextEnd,InsSeqAlignments");

                mSvFileWriter.write(",FSStart,FSEnd,LEStart,LEEnd"); // ,DupBEStart,DupBEEnd

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

                // gene & replication info
                mSvFileWriter.write(",GeneStart,GeneEnd,RepOriginStart,RepOriginEnd");

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

                final String insSeqAlignments = dbData.insertSequenceAlignments().replaceAll(",", ";");

                writer.write(
                        String.format(",%s,%d,%d,%s,%s,%.0f,%s,%s,%s",
                                dbData.insertSequence().isEmpty() && var.type() != INS ? dbData.homology() : "",
                                dbData.inexactHomologyOffsetStart(), dbData.inexactHomologyOffsetEnd(),
                                dbData.insertSequence(), dbData.imprecise(), dbData.qualityScore(),
                                dbData.startRefContext(), dbData.endRefContext(), insSeqAlignments));

                writer.write(
                        String.format(",%s,%s,%s,%s",
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

                writer.write(String.format(",%s,%s,%.4f,%.4f",
                        var.getGeneInBreakend(true), var.getGeneInBreakend(false),
                        var.getReplicationOrigin(true), var.getReplicationOrigin(false)));

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
                        String.format(",%d,%d,%d,%d,%d,%.2f,%.2f",
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


    public void close()
    {
        closeBufferedWriter(mSvFileWriter);
        closeBufferedWriter(mClusterFileWriter);
        closeBufferedWriter(mLinksFileWriter);
        mVisWriter.close();

        // log perf stats
        mPerfCounter.stop();
        mPerfCounter.logStats(false);
        // mPc1.logStats(false);
        mPc2.logStats(false);
        mPc3.logStats(false);
        // mPc4.logStats(false);
        mPc5.logStats(false);

        mAnalyser.logStats();
    }

}
