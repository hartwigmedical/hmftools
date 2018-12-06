package com.hartwig.hmftools.svanalysis.analysis;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.LinkFinder.NO_DB_MARKER;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.common.variant.structural.annotation.SvPONAnnotator.REGION_DISTANCE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_LOW_QUALITY;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SIMPLE_SV;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator;
import com.hartwig.hmftools.common.variant.structural.annotation.SvPONAnnotator;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.common.variant.structural.annotation.SvPON;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class SvSampleAnalyser {

    private final SvaConfig mConfig;
    private final ClusterAnalyser mAnalyser;

    // data per run (ie sample)
    private String mSampleId;
    private List<SvVarData> mAllVariants; // the original list to analyse

    BufferedWriter mSVFileWriter;
    BufferedWriter mClusterFileWriter;
    BufferedWriter mLinksFileWriter;
    SvPONAnnotator mSvPONAnnotator;
    FragileSiteAnnotator mFragileSiteAnnotator;
    LineElementAnnotator mLineElementAnnotator;
    SvClusteringMethods mClusteringMethods;
    CNAnalyser mCnAnalyser;

    PerformanceCounter mPerfCounter;
    PerformanceCounter mPc1;
    PerformanceCounter mPc2;
    PerformanceCounter mPc3;
    PerformanceCounter mPc4;
    PerformanceCounter mPc5;

    private static final Logger LOGGER = LogManager.getLogger(SvSampleAnalyser.class);

    public SvSampleAnalyser(final SvaConfig config)
    {
        mConfig = config;
        mClusteringMethods = new SvClusteringMethods(mConfig.ProximityDistance);
        mAnalyser = new ClusterAnalyser(config, mClusteringMethods);

        mSVFileWriter = null;
        mLinksFileWriter = null;
        mClusterFileWriter = null;

        mSvPONAnnotator = new SvPONAnnotator();
        mSvPONAnnotator.loadPonFile(mConfig.SvPONFile);

        mFragileSiteAnnotator = new FragileSiteAnnotator();
        mFragileSiteAnnotator.loadFragileSitesFile(mConfig.FragileSiteFile);

        mLineElementAnnotator = new LineElementAnnotator();
        mLineElementAnnotator.loadLineElementsFile(mConfig.LineElementFile);

        mCnAnalyser = null;
        if(!mConfig.LOHDataFile.isEmpty())
        {
            mCnAnalyser = new CNAnalyser(mConfig.OutputCsvPath, null);
            mCnAnalyser.loadLOHFromCSV(mConfig.LOHDataFile, "");
            mClusteringMethods.setSampleLohData(mCnAnalyser.getSampleLohData());
        }

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

        mPc5.start();

        if(!mConfig.OutputCsvPath.isEmpty())
        {
            writeClusterSVOutput();

            if(mConfig.hasMultipleSamples())
            {
                writeClusterLinkData();
                writeClusterData();
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

            setPonOccurenceCount(var);
            var.setFragileSites(mFragileSiteAnnotator.isFragileSite(var, true), mFragileSiteAnnotator.isFragileSite(var, false));
            var.setLineElement(mLineElementAnnotator.isLineElement(var, true), true);
            var.setLineElement(mLineElementAnnotator.isLineElement(var, false), false);

            // exclude PON
            if (var.getPonCount() >= 2)
            {
                LOGGER.info("filtering sv({}) with PonCount({})", var.id(), var.getPonCount());
                mAllVariants.remove(currentIndex);
                continue;
            }

            String startArm = getChromosomalArm(var.chromosome(true), var.position(true));

            String endArm = "";
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
            BufferedWriter writer = null;

            if(mSVFileWriter != null)
            {
                // check if can continue appending to an existing file
                writer = mSVFileWriter;
            }
            else
            {
                String outputFileName = mConfig.OutputCsvPath;

                if(!outputFileName.endsWith("/"))
                    outputFileName += File.separator;

                if(mConfig.hasMultipleSamples())
                    outputFileName += "CLUSTER.csv";
                else
                    outputFileName += mSampleId + ".csv";

                Path outputFile = Paths.get(outputFileName);

                writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);
                mSVFileWriter = writer;

                // definitional fields
                writer.write("SampleId,Id,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd");

                // position and copy number
                writer.write(",ArmStart,AdjAFStart,AdjCNStart,AdjCNChgStart,ArmEnd,AdjAFEnd,AdjCNEnd,AdjCNChgEnd,Ploidy");

                // cluster info
                writer.write(",ClusterId,SubClusterId,ClusterCount,ClusterReason");

                // cluster-level info
                writer.write(",ClusterDesc,IsResolved,ResolvedType,Consistency,ArmCount");

                // SV info
                writer.write(",Homology,InexactHOStart,InexactHOEnd,InsertSeq,Imprecise,PONCount,QualScore,RefContextStart,RefContextEnd");

                // location attributes
                writer.write(",FSStart,FSEnd,LEStart,LEEnd,DupBEStart,DupBEEnd,ArmCountStart,ArmExpStart,ArmCountEnd,ArmExpEnd");

                // linked pair info
                writer.write(",LnkSvStart,LnkTypeStart,LnkLenStart,LnkSvEnd,LnkTypeEnd,LnkLenEnd");

                // GRIDDS caller info
                writer.write(",AsmbStart,AsmbEnd,AsmbMatchStart,AsmbMatchEnd");

                // chain info
                writer.write(",ChainId,ChainCount,ChainIndex");

                // proximity info and other link info
                writer.write(",NearestLen,NearestType,DBLenStart,DBLenEnd,SynDelDupLen,SynDelDupTILen");

                // proximity info and other link info
                writer.write(",FoldbackLnkStart,FoldbackLenStart,FoldbackLinkInfoStart,FoldbackLnkEnd,FoldbackLenEnd,FoldbackLinkInfoEnd");

                // transitive info
                // writer.write(",TransType,TransLen,TransSvLinks");

                writer.newLine();
            }

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
                        String.format(",%s,%d,%d,%s,%s,%d,%.0f,%s,%s",
                                dbData.insertSequence().isEmpty() && var.type() != INS ? dbData.homology() : "",
                                dbData.inexactHomologyOffsetStart(), dbData.inexactHomologyOffsetEnd(),
                                dbData.insertSequence(), dbData.imprecise(), var.getPonCount(), dbData.qualityScore(),
                                dbLenStart > NO_DB_MARKER & dbLenStart < 0 ? dbData.startRefContext() : "",
                                dbLenEnd > NO_DB_MARKER & dbLenEnd < 0 ? dbData.endRefContext() : ""));

                writer.write(
                        String.format(",%s,%s,%s,%s,%s,%s,%s",
                                var.isFragileSite(true), var.isFragileSite(false),
                                var.getLineElement(true), var.getLineElement(false),
                                var.isDupBreakend(true), var.isDupBreakend(false),
                                mClusteringMethods.getChrArmData(var)));

                // linked pair info
                final SvLinkedPair startLP = var.getLinkedPair(true) != null ? var.getLinkedPair(true) : cluster.getLinkedPair(var, true);
                String startLinkStr = "0,,-1";
                String assemblyMatchStart = var.getAssemblyMatchType(true);
                if(startLP != null)
                {
                    startLinkStr = String.format("%s,%s,%d",
                            startLP.first().equals(var, true) ? startLP.second().origId() : startLP.first().origId(),
                            startLP.linkType(), startLP.length());
                }

                final SvLinkedPair endLP = var.getLinkedPair(false) != null ? var.getLinkedPair(false) : cluster.getLinkedPair(var, false);
                String endLinkStr = "0,,-1";
                String assemblyMatchEnd = var.getAssemblyMatchType(false);
                if(endLP != null)
                {
                    endLinkStr = String.format("%s,%s,%d",
                            endLP.first().equals(var, true) ? endLP.second().origId() : endLP.first().origId(),
                            endLP.linkType(), endLP.length());
                }

                // assembly info
                writer.write(String.format(",%s,%s,%s,%s,%s,%s",
                        startLinkStr, endLinkStr, var.getAssemblyData(true), var.getAssemblyData(false), assemblyMatchStart, assemblyMatchEnd));

                // chain info
                final SvChain chain = cluster.findChain(var);
                String chainStr = ",0,0,";

                if(chain != null)
                {
                    chainStr = String.format(",%d,%d,%s", chain.id(), chain.getUniqueSvCount(), chain.getSvIndices(var));
                }

                writer.write(chainStr);

                writer.write(String.format(",%d,%s,%d,%d,%d,%d",
                        var.getNearestSvDistance(), var.getNearestSvRelation(), dbLenStart, dbLenEnd,
                        cluster.getSynDelDupLength(), cluster.getSynDelDupTILength()));

                writer.write(String.format(",%s,%d,%s,%s,%d,%s",
                        var.getFoldbackLink(true), var.getFoldbackLen(true), var.getFoldbackLinkInfo(true),
                        var.getFoldbackLink(false), var.getFoldbackLen(false), var.getFoldbackLinkInfo(false)));

                // writer.write(String.format(",%s,%d,%s", var.getTransType(), var.getTransLength(), var.getTransSvLinks()));

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
            BufferedWriter writer = null;

            if(mClusterFileWriter != null)
            {
                // check if can continue appending to an existing file
                writer = mClusterFileWriter;
            }
            else
            {
                String outputFileName = mConfig.OutputCsvPath;

                if(!outputFileName.endsWith("/"))
                    outputFileName += File.separator;

                outputFileName += "SVA_CLUSTERS.csv";

                Path outputFile = Paths.get(outputFileName);

                writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);
                mClusterFileWriter = writer;

                writer.write("SampleId,ClusterId,ClusterDesc,ClusterCount,ResolvedType,FullyChained,ChainCount");
                writer.write(",DelCount,DupCount,InsCount,InvCount,BndCount,SglCount");
                writer.write(",Consistency,ArmCount,IsLINE,HasReplicated");
                writer.write(",AssemblyLinks,LongDelDups,UnlinkedRemotes,ShortTIRemotes,Annotations,ChainInfo");
                writer.newLine();
            }

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
                        String.format(",%d,%d,%s,%s",
                                cluster.getConsistencyCount(), cluster.getArmCount(),
                                cluster.hasLinkingLineElements(), cluster.hasReplicatedSVs()));

                final String chainInfo = cluster.getChains().stream()
                        .filter(x -> !x.getDetails().isEmpty())
                        .map(x -> x.getDetails())
                        .collect (Collectors.joining (";"));

                writer.write(
                        String.format(",%d,%d,%d,%d,%s,%s",
                                cluster.getAssemblyLinkedPairs().size(), cluster.getLongDelDups().size(),
                                cluster.getUnlinkedRemoteSVs().size(), cluster.getShortTIRemoteSVs().size(),
                                cluster.getAnnotations(), chainInfo));

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
            BufferedWriter writer = null;

            if(mLinksFileWriter != null)
            {
                // check if can continue appending to an existing file
                writer = mLinksFileWriter;
            }
            else
            {
                String outputFileName = mConfig.OutputCsvPath;

                if(!outputFileName.endsWith("/"))
                    outputFileName += File.separator;

                outputFileName += "SVA_LINKS.csv";

                Path outputFile = Paths.get(outputFileName);

                writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);
                mLinksFileWriter = writer;

                writer.write("SampleId,ClusterId,ClusterDesc,ClusterCount,ResolvedType,IsLINE,FullyChained");
                writer.write(",ChainId,ChainCount,ChainConsistent,Id1,Id2,ChrArm");
                writer.write(",IsAssembled,TILength,NextSVLength,DBLenStart,DBLenEnd,OnArmOfOrigin,CopyNumberGain,TraversedSVCount,AssembledCount");
                writer.newLine();
            }

            for(final SvCluster cluster : getClusters())
            {
                int clusterSvCount = cluster.getUniqueSvCount();

                isSpecificCluster(cluster);

                List<SvChain> chains = null;

                if(cluster.hasLinkingLineElements())
                {
                    // line elements currently aren't chained, so manufacture one for the looping below
                    SvChain tempChain = new SvChain(0);
                    tempChain.getLinkedPairs().addAll(cluster.getAssemblyLinkedPairs());

                    chains = Lists.newArrayList();
                    chains.add(tempChain);
                }
                else
                {
                    chains = cluster.getChains();
                }

                for (final SvChain chain : chains)
                {
                    int chainSvCount = chain.getSvCount();
                    boolean chainConsistent = !cluster.hasLinkingLineElements() ? chain.isConsistent() : false;

                    for (final SvLinkedPair pair : chain.getLinkedPairs())
                    {
                        if(pair.linkType() != LINK_TYPE_TI)
                            continue;

                        if(pair.first().type() == SGL || pair.second().type() == SGL)
                            continue;

                        writer.write(
                                String.format("%s,%d,%s,%d,%s,%s,%s",
                                        mSampleId, cluster.id(), cluster.getDesc(), clusterSvCount, cluster.getResolvedType(),
                                        cluster.hasLinkingLineElements(), cluster.isFullyChained()));

                        writer.write(
                                String.format(",%d,%d,%s,%s,%s,%s",
                                        chain.id(), chainSvCount, chainConsistent,
                                        pair.first().origId(), pair.second().origId(),
                                        pair.first().getBreakend(pair.firstLinkOnStart()).getChrArm()));

                        writer.write(
                                String.format(",%s,%d,%d,%d,%d,%s,%s,%d,%d",
                                        pair.isAssembled(), pair.length(), pair.getNearestSVDistance(),
                                        pair.getDBLenFirst(), pair.getDBLenSecond(), pair.onArmOfOrigin(),
                                        pair.hasCopyNumberGain(), pair.getTraversedSVCount(), pair.getAssembledChainCount()));

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

    public void setPonOccurenceCount(SvVarData svData)
    {
        if(mSvPONAnnotator == null || !mSvPONAnnotator.hasEntries())
            return;

        if(mSvPONAnnotator.getPonList().isEmpty())
            return;

        final List<SvPON> ponList = mSvPONAnnotator.getPonList();
        final Map<String, Integer> chrIndexMap = mSvPONAnnotator.getChrIndexMap();

        // use CRMS to find start index
        if(!chrIndexMap.containsKey(svData.chromosome(true)))
            return;

        int index = chrIndexMap.get(svData.chromosome(true));

        for(; index < ponList.size(); ++index)
        {
            final SvPON svPon = ponList.get(index);

            if(!svPon.chrStart().equals(svData.chromosome(true)))
            {
                return;
            }

            if(isMatch(svData, svPon,0))
            {
                LOGGER.debug("var({}) found in PON with count({})", svData.posId(), svPon.count());
                svData.setPonCount(svPon.count());
                svData.setPonRegionCount(svPon.count()); // set to match
                return;
            }
            else if(svData.getPonRegionCount() == 0 && isMatch(svData, svPon, REGION_DISTANCE))
            {
                LOGGER.debug("var({}) found in PON region({}->{}) with count({})",
                        svData.posId(), svPon.posStart(), svPon.posEnd(), svPon.count());

                svData.setPonRegionCount(svPon.count());
            }
        }
    }

    private static boolean isMatch(final SvVarData svData, final SvPON svPon, final int permittedDiff)
    {
        if(!svPon.chrStart().equals(svData.chromosome(true))
        || !svPon.chrEnd().equals(svData.chromosome(false))
        || svPon.orientStart() != svData.orientation(true)
        || svPon.orientEnd() != svData.orientation(false)
        || !svPon.type().equals(svData.type().toString()))
        {
            return false;
        }

        return Math.abs(svPon.posStart() - svData.position(true)) <= permittedDiff
                && Math.abs(svPon.posEnd() - svData.position(false)) <= permittedDiff;
    }

    private void logSampleClusterInfo()
    {
        int simpleClusterCount = 0;
        int complexClusterCount = 0;
        Map<String, Integer> armSimpleClusterCount = new HashMap();
        Map<String, Integer> armComplexClusterCount = new HashMap();

        for(final SvCluster cluster : mAnalyser.getClusters())
        {
            Map<String, Integer> targetMap;

            if(cluster.getResolvedType() == RESOLVED_LOW_QUALITY)
                continue;

            if(cluster.getResolvedType() == RESOLVED_TYPE_SIMPLE_SV || cluster.isSyntheticSimpleType())
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
        try
        {
            if(mSVFileWriter != null)
                mSVFileWriter.close();

            if(mLinksFileWriter != null)
                mLinksFileWriter.close();

            if(mClusterFileWriter != null)
                mClusterFileWriter.close();
        }
        catch (final IOException e)
        {
        }

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
