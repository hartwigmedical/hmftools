package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM_CHR;
import static com.hartwig.hmftools.linx.analysis.ClusterClassification.getClusterCategory;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.asStr;
import static com.hartwig.hmftools.linx.types.ArmCluster.getArmClusterData;
import static com.hartwig.hmftools.linx.types.LinxConstants.NO_DB_MARKER;
import static com.hartwig.hmftools.linx.types.LinxConstants.SHORT_TI_LENGTH;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.linx.analysis.ClusterMetrics;
import com.hartwig.hmftools.linx.annotators.LineElementType;
import com.hartwig.hmftools.linx.chaining.ChainMetrics;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.ArmClusterType;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.ArmCluster;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisDataWriter;

// a class to own all the cohort-level writers for multi-threading synchronicity
public class CohortDataWriter
{
    private final LinxConfig mConfig;

    private final VisDataWriter mVisWriter;

    private final Map<String,BufferedWriter> mWriters;

    public static final String COHORT_WRITER_SV = "SvData";
    public static final String COHORT_WRITER_CLUSTER = "Cluster";
    public static final String COHORT_WRITER_LINK = "Link";

    public CohortDataWriter(final LinxConfig config, final EnsemblDataCache geneDataCache)
    {
        mConfig = config;

        mWriters = Maps.newHashMap();

        mVisWriter = new VisDataWriter(
                config.OutputDataPath, geneDataCache, config.Output.writeVisualisationData(), config.hasMultipleSamples(), config.IsGermline);

        if(mConfig.Output.WriteTypes.contains(WriteType.SV_DATA))
            mWriters.put(COHORT_WRITER_SV, createSvDataFile());

        if(mConfig.Output.WriteTypes.contains(WriteType.CLUSTER))
            mWriters.put(COHORT_WRITER_CLUSTER, createClusterFile());

        if(mConfig.Output.WriteTypes.contains(WriteType.LINK))
            mWriters.put(COHORT_WRITER_LINK, createLinksFile());
    }

    public final VisDataWriter getVisWriter() { return mVisWriter; }

    public void close()
    {
        mWriters.values().forEach(x -> closeBufferedWriter(x));
        mVisWriter.close();
    }

    public boolean writeCohortFiles()
    {
        return mConfig.hasMultipleSamples() || mConfig.Output.WriteCohortFiles;
    }

    public boolean hasWriter(final String fileType) { return mWriters.containsKey(fileType); }

    public synchronized void write(final CohortFileInterface cohortFile, final List<String> lines)
    {
        if(mConfig.OutputDataPath == null)
            return;

        BufferedWriter writer = mWriters.get(cohortFile.fileType());

        if(writer == null)
        {
            writer = cohortFile.createWriter(mConfig.OutputDataPath);
            mWriters.put(cohortFile.fileType(), writer);
        }

        if(writer == null)
            return;

        try
        {
            for(String line : lines)
            {
                writer.write(line);
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to write to {} file: {}", cohortFile.fileType(), e.toString());
        }
    }

    public static String cohortDataFilename(final String outputDir, final String fileId)
    {
        return outputDir + "LNX_" + fileId + TSV_EXTENSION;
    }

    private BufferedWriter createSvDataFile()
    {
        if(!writeCohortFiles())
            return null;

        try
        {
            String outputFileName = cohortDataFilename(mConfig.OutputDataPath, "SVS");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            // definitional fields
            writer.write("SampleId\tId\tType\tClusterId\tClusterCount");
            writer.write("\tChrStart\tPosStart\tOrientStart\tArmStart\tChrEnd\tPosEnd\tOrientEnd\tArmEnd");

            // position and copy number
            writer.write("\tCNStart\tCNChgStart\tCNEnd\tCNChgEnd\tJcn\tJcnMin\tJcnMax");

            // cluster info
            writer.write("\tClusterReason\tClusterDesc\tResolvedType");

            writer.write("\tFSStart\tFSEnd\tLEStart\tLEEnd");

            // linked pair info
            writer.write("\tLnkSvStart\tLnkLenStart\tLnkSvEnd\tLnkLenEnd\tAsmbStart\tAsmbEnd");

            // chain info
            writer.write("\tChainId\tChainCount\tChainIndex");

            // proximity info and other link info
            writer.write("\tNearestLen\tNearestType\tDBLenStart\tDBLenEnd");

            // proximity info and other link info
            writer.write("\tFoldbackLnkStart\tFoldbackLenStart\tFoldbackInfoStart\tFoldbackLnkEnd\tFoldbackLenEnd\tFoldbackInfoEnd");

            // local topology from arm cluster
            writer.write("\tLocTopIdStart\tLocTopTypeStart\tLocTopTIStart\tLocTopIdEnd\tLocTopTypeEnd\tLocTopTIEnd");

            // gene info
            writer.write("\tGeneStart\tGeneEnd");

            if(mConfig.Output.writeSvData())
            {
                // extra copy number info
                writer.write("\tMinorAPStartPrev\tMinorAPStartPost\tMinorAPEndPrev\tMinorAPEndPost\tAFStart\tAFEnd");

                // SV table info
                writer.write("\tHomologyStart\tHomologyEnd\tInsertSeq\tQualScore");
                writer.write("\tInsSeqAlignments");
                writer.write("\tRepeatClass\tRepeatType\tAnchorStart\tAnchorEnd");
            }

            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to open and write output file headers");
            return null;
        }
    }

    private static final int INF_DB_MARKER = -2000;

    public synchronized void writeSvData(final String sampleId, final List<SvVarData> svDataList)
    {
        BufferedWriter writer = mWriters.get(COHORT_WRITER_SV);

        if(writer == null)
            return;

        try
        {
            for(final SvVarData var : svDataList)
            {
                final SvCluster cluster = var.getCluster();

                if(cluster == null)
                {
                    LNX_LOGGER.error("SV({}) not assigned to any cluster", var.posId());
                    continue;
                }

                if(mConfig.IsGermline && var.getGenesList(true).isEmpty() && var.getGenesList(false).isEmpty())
                    continue;

                final StructuralVariantData dbData = var.getSvData();

                final ArmCluster armClusterStart = cluster.findArmCluster(var.getBreakend(true));

                final ArmCluster armClusterEnd = !var.isSglBreakend() ? cluster.findArmCluster(var.getBreakend(false)) : null;

                writer.write(String.format("%s\t%d\t%s\t%d\t%d",
                        sampleId, var.id(), var.typeStr(), cluster.id(), cluster.getSvCount()));

                writer.write(String.format("\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s",
                        var.chromosome(true), var.position(true), var.orientation(true), asStr(var.arm(true)),
                        var.chromosome(false), var.position(false), var.orientation(false), asStr(var.arm(false))));

                writer.write(String.format("\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f",
                        dbData.adjustedStartCopyNumber(), dbData.adjustedStartCopyNumberChange(),
                        dbData.adjustedEndCopyNumber(), dbData.adjustedEndCopyNumberChange(),
                        dbData.junctionCopyNumber(), var.jcnMin(), var.jcnMax()));

                writer.write(String.format("\t%s\t%s\t%s",
                        var.getClusterReason(), cluster.getDesc(), cluster.getResolvedType()));

                writer.write(String.format("\t%s\t%s\t%s\t%s",
                        var.isFragileSite(true), var.isFragileSite(false),
                        LineElementType.toString(var.getLineElement(true)), LineElementType.toString(var.getLineElement(false))));

                // linked pair info
                for(int be = SE_START; be <= SE_END; ++be)
                {
                    boolean isStart = isStart(be);
                    final LinkedPair link = var.getLinkedPair(isStart);
                    if(link != null)
                    {
                        writer.write(String.format("\t%s\t%d",
                                link.first() == var ? link.second().id() : link.first().id(), link.baseLength()));
                    }
                    else
                    {
                        writer.write("\t\t-1");
                    }
                }

                // assembly info
                writer.write(String.format("\t%s\t%s", var.assemblyInfoStr(true), var.assemblyInfoStr(false)));

                // chain info
                final SvChain chain = cluster.findChain(var);
                String chainStr = "";

                if(chain != null)
                {
                    chainStr = String.format("\t%d\t%d\t%s", chain.id(), chain.getSvCount(), chain.getSvIndices(var));
                }
                else
                {
                    chainStr = String.format("\t%d\t0\t", cluster.getChainId(var));
                }

                writer.write(chainStr);

                // only log DB lengths if the partner is in the cluster
                final DbPair dbStart = var.getDBLink(true);
                final DbPair dbEnd = var.getDBLink(false);

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

                writer.write(String.format("\t%d\t%s\t%d\t%d",
                        var.getNearestSvDistance(), var.getNearestSvRelation(), dbLenStart, dbLenEnd));

                writer.write(String.format("\t%d\t%d\t%s\t%d\t%d\t%s",
                        var.getFoldbackId(true), var.getFoldbackLength(true), var.getFoldbackInfo(true),
                        var.getFoldbackId(false), var.getFoldbackLength(false), var.getFoldbackInfo(false)));

                for(int be = SE_START; be <= SE_END; ++be)
                {
                    ArmCluster armCluster = be == SE_START ? armClusterStart : armClusterEnd;

                    if(armCluster != null)
                        writer.write(String.format("\t%d\t%s\t%d", armCluster.id(), armCluster.getTypeStr(), armCluster.getTICount()));
                    else
                        writer.write("\t-1\t\t0");
                }

                writer.write(String.format("\t%s\t%s",
                        var.getGeneInBreakend(true, true), var.getGeneInBreakend(false, true)));

                if(mConfig.Output.writeSvData())
                {
                    writer.write(String.format("\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f",
                            var.getBreakend(true).minorAlleleJcn(true),
                            var.getBreakend(true).minorAlleleJcn(false),
                            !var.isSglBreakend() ? var.getBreakend(false).minorAlleleJcn(true) : 0,
                            !var.isSglBreakend() ? var.getBreakend(false).minorAlleleJcn(false) : 0,
                            dbData.adjustedStartAF(), dbData.adjustedEndAF()));

                    final String insSeqAlignments = dbData.insertSequenceAlignments().replaceAll(",", ";");

                    writer.write(String.format("\t%s\t%s\t%s\t%.0f\t%s",
                            dbData.startHomologySequence(), dbData.endHomologySequence(),
                            dbData.insertSequence(), dbData.qualityScore(), insSeqAlignments));

                    writer.write(String.format("\t%s\t%s\t%d\t%d",
                            dbData.insertSequenceRepeatClass(), dbData.insertSequenceRepeatType(),
                            dbData.startAnchoringSupportDistance(), dbData.endAnchoringSupportDistance()));
                }

                writer.newLine();
            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    private BufferedWriter createClusterFile()
    {
        if(!writeCohortFiles())
            return null;

        try
        {
            String outputFileName = cohortDataFilename(mConfig.OutputDataPath, "CLUSTERS");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("SampleId\tClusterId\tClusterDesc\tClusterCount\tCategory\tResolvedType\tSynthetic\tFullyChained\tChainCount");
            writer.write("\tDelCount\tDupCount\tInsCount\tInvCount\tBndCount\tSglCount\tInfCount");
            writer.write("\tClusterReasons\tConsistency\tIsLINE\tReplication\tMinJcn\tMaxJcn\tFoldbacks");
            writer.write("\tArmCount\tOriginArms\tFragmentArms\tConsistentArms\tComplexArms\tAnnotations\tAlleleValidPerc");

            writer.write("\tTotalTIs\tAssemblyTIs\tShortTIs\tIntTIs\tExtTIs\tIntShortTIs\tExtShortTIs\tIntTIsCnGain");
            writer.write("\tExtTIsCnGain\tOverlapTIs\tChainEndsFace\tChainEndsAway\tUnchainedSVs");

            writer.write("\tDBs\tShortDBs\tTotalDBLength\tTotalDeleted\tTravDelCount\tTravDelLength");
            writer.write("\tTotalRange\tChainedLength\tImpliedTIs");

            writer.write("\tArmClusterCount\tAcTotalTIs\tAcIsolatedBE\tAcTIOnly\tAcDsb\tAcSimpleDup");
            writer.write("\tAcSingleFb\tAcFbDsb\tAcComplexFb\tAcComplexLine\tAcSameOrient\tAcComplexOther");

            writer.newLine();
            return writer;
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing cluster-data to outputFile: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeClusterData(final String sampleId, final List<SvCluster> clusters)
    {
        BufferedWriter writer = mWriters.get(COHORT_WRITER_CLUSTER);

        if(writer == null)
            return;

        try
        {
            for(final SvCluster cluster : clusters)
            {
                int clusterSvCount = cluster.getSvCount();

                if(clusterSvCount == 1 && !mConfig.Output.WriteSingleSVClusters)
                    continue;

                if(mConfig.IsGermline)
                {
                    if(cluster.getSVs().stream()
                            .noneMatch(x -> !x.getGenesList(true).isEmpty() || !x.getGenesList(false).isEmpty()))
                    {
                        continue;
                    }
                }

                ResolvedType resolvedType = cluster.getResolvedType();

                final String category = getClusterCategory(cluster);

                writer.write(String.format("%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%d",
                        sampleId, cluster.id(), cluster.getDesc(), clusterSvCount,
                        category, resolvedType, cluster.isSyntheticType(),
                        cluster.isFullyChained(false), cluster.getChains().size()));

                writer.write(String.format("\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
                        cluster.getTypeCount(DEL), cluster.getTypeCount(DUP), cluster.getTypeCount(INS),
                        cluster.getTypeCount(INV), cluster.getTypeCount(BND), cluster.getTypeCount(SGL), cluster.getTypeCount(INF)));

                double foldbackCount = 0;

                for(final SvVarData var : cluster.getFoldbacks())
                {
                    // avoid double-count chained foldbacks
                    if(var.getFoldbackBreakend(true) != null)
                        foldbackCount += 0.5;
                    if(var.getFoldbackBreakend(false) != null)
                        foldbackCount += 0.5;
                }

                writer.write(String.format("\t%s\t%d\t%s\t%s\t%.1f\t%.1f\t%.0f",
                        cluster.getClusteringReasons(), cluster.getConsistencyCount(), cluster.hasLinkingLineElements(),
                        cluster.requiresReplication(), cluster.getMinJcn(), cluster.getMaxJcn(), foldbackCount));

                final ClusterMetrics metrics = cluster.getMetrics();

                writer.write(String.format("\t%d\t%d\t%d\t%d\t%d\t%s\t%.2f",
                        cluster.getArmCount(), metrics.OriginArms, metrics.FragmentArms, metrics.ConsistentArms,
                        metrics.ComplexArms, cluster.getAnnotations(), metrics.ValidAlleleJcnSegmentPerc));

                long shortTIs = cluster.getLinkedPairs().stream().filter(x -> x.baseLength() <= SHORT_TI_LENGTH).count();

                writer.write(String.format("\t%d\t%d\t%d",
                        cluster.getLinkedPairs().size(), cluster.getAssemblyLinkedPairs().size(), shortTIs));

                final ChainMetrics chainMetrics = cluster.getLinkMetrics();

                writer.write(String.format("\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
                        chainMetrics.InternalTIs, chainMetrics.ExternalTIs, chainMetrics.InternalShortTIs, chainMetrics.ExternalShortTIs,
                        chainMetrics.InternalTICnGain, chainMetrics.ExternalTICnGain, chainMetrics.OverlappingTIs,
                        chainMetrics.ChainEndsFace, chainMetrics.ChainEndsAway, cluster.getUnlinkedSVs().size()));

                writer.write(String.format("\t%d\t%d\t%d\t%d\t%d\t%d",
                        metrics.DBCount, metrics.ShortDBCount, metrics.TotalDBLength,
                        metrics.TotalDeleted, metrics.TraversedDelCount, metrics.TraversedDelLength));

                writer.write(String.format("\t%d\t%d\t%d",
                        metrics.TotalRange, metrics.ChainedLength, metrics.ImpliedTICount));

                final int[] armClusterData = getArmClusterData(cluster);
                long armClusterTIs = cluster.getArmClusters().stream().mapToInt(x -> x.getTICount()).sum();

                writer.write(String.format("\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
                        cluster.getArmClusters().size(), armClusterTIs, armClusterData[ArmClusterType.ISOLATED_BE.ordinal()],
                        armClusterData[ArmClusterType.TI_ONLY.ordinal()], armClusterData[ArmClusterType.DSB.ordinal()],
                        armClusterData[ArmClusterType.SIMPLE_DUP.ordinal()], armClusterData[ArmClusterType.FOLDBACK.ordinal()],
                        armClusterData[ArmClusterType.FOLDBACK_DSB.ordinal()], armClusterData[ArmClusterType.COMPLEX_FOLDBACK.ordinal()],
                        armClusterData[ArmClusterType.COMPLEX_LINE.ordinal()], armClusterData[ArmClusterType.SAME_ORIENT.ordinal()],
                        armClusterData[ArmClusterType.COMPLEX_OTHER.ordinal()]));

                writer.newLine();

            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing cluster-data to outputFile: {}", e.toString());
        }
    }

    private BufferedWriter createLinksFile()
    {
        if(!mConfig.Output.writeLinks() || !writeCohortFiles())
            return null;

        try
        {
            String outputFileName = cohortDataFilename(mConfig.OutputDataPath, "LINKS");;

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("SampleId\tClusterId\tClusterCount\tResolvedType");
            writer.write("\tChainId\tChainCount\tChainConsistent\tLinkReason\tLinkIndex\tChainIndex\tJcn\tJcnUncertainty");
            writer.write("\tIsAssembled\tTILength\tNextSvDist\tNextClusteredSvDist\tTraversedSVCount");
            writer.write("\tLocationType\tOverlapCount\tCopyNumberGain");
            writer.write("\tId1\tId2\tChrArm\tPosStart\tPosEnd\tLocTopTypeStart\tLocTopTypeEnd\tGeneStart\tGeneEnd\tExonMatch\tIsDM");
            writer.newLine();

            return writer;
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing links to outputFile: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeLinksData(final String sampleId, final List<SvCluster> clusters)
    {
        if(!mConfig.Output.writeLinks())
            return;

        BufferedWriter writer = mWriters.get(COHORT_WRITER_LINK);

        if(writer == null)
            return;

        try
        {
            for(final SvCluster cluster : clusters)
            {
                int clusterSvCount = cluster.getSvCount();

                List<SvChain> chains = cluster.getChains();

                for(final SvChain chain : chains)
                {
                    int chainSvCount = chain.getSvCount();
                    boolean chainConsistent = chain.isConsistent();
                    boolean isDoubleMinute = chain.isDoubleMinute();

                    List<LinkedPair> uniquePairs = Lists.newArrayList();
                    final List<LinkedPair> chainLinks = chain.getLinkedPairs();

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
                                chainIndexStr = appendStr(chainIndexStr, String.valueOf(j), ITEM_DELIM_CHR);
                            }
                        }

                        final SvBreakend beStart = pair.getBreakend(true);
                        final SvBreakend beEnd = pair.getBreakend(false);

                        writer.write(String.format("%s\t%d\t%d\t%s",
                                sampleId, cluster.id(), clusterSvCount, cluster.getResolvedType()));

                        writer.write(String.format("\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%.3f",
                                chain.id(), chainSvCount, chainConsistent, pair.getLinkReason(), pair.getLinkIndex(),
                                chainIndexStr, formatJcn(chain.jcn()), chain.jcnUncertainty()));

                        writer.write(String.format("\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%s",
                                pair.isAssembled(), pair.baseLength(),
                                pair.getNextSvDistance(), pair.getNextClusteredSvDistance(), pair.getTraversedSVCount(),
                                pair.locationType(), pair.overlapCount(), pair.hasCopyNumberGain()));

                        ArmCluster acStart = cluster.findArmCluster(beStart);
                        ArmCluster acEnd = cluster.findArmCluster(beEnd);

                        writer.write(String.format("\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
                                beStart.getSV().id(), beEnd.getSV().id(),
                                beStart.getChrArm(), beStart.position(), beEnd.position(),
                                acStart != null ? acStart.getTypeStr() : "", acEnd != null ? acEnd.getTypeStr() : "",
                                beStart.getSV().getGeneInBreakend(beStart.usesStart(), false),
                                beEnd.getSV().getGeneInBreakend(beEnd.usesStart(), false),
                                pair.getExonMatchData(), isDoubleMinute));

                        writer.newLine();
                    }
                }
            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing links to outputFile: {}", e.toString());
        }
    }
}
