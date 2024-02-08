package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
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
                config.OutputDataPath, geneDataCache, config.Output.WriteVisualisationData, config.hasMultipleSamples(), config.IsGermline);

        mWriters.put(COHORT_WRITER_SV, createSvDataFile());
        mWriters.put(COHORT_WRITER_CLUSTER, createClusterFile());
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

    private BufferedWriter createSvDataFile()
    {
        if(!writeCohortFiles())
            return null;

        try
        {
            String outputFileName = mConfig.OutputDataPath + "LNX_SVS.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            // definitional fields
            writer.write("SampleId,Id,Type,ClusterId,ClusterCount");
            writer.write(",ChrStart,PosStart,OrientStart,ArmStart,ChrEnd,PosEnd,OrientEnd,ArmEnd");

            // position and copy number
            writer.write(",CNStart,CNChgStart,CNEnd,CNChgEnd,Jcn,JcnMin,JcnMax");

            // cluster info
            writer.write(",ClusterReason,ClusterDesc,ResolvedType");

            writer.write(",FSStart,FSEnd,LEStart,LEEnd");

            // linked pair info
            writer.write(",LnkSvStart,LnkLenStart,LnkSvEnd,LnkLenEnd,AsmbStart,AsmbEnd");

            // chain info
            writer.write(",ChainId,ChainCount,ChainIndex");

            // proximity info and other link info
            writer.write(",NearestLen,NearestType,DBLenStart,DBLenEnd");

            // proximity info and other link info
            writer.write(",FoldbackLnkStart,FoldbackLenStart,FoldbackInfoStart,FoldbackLnkEnd,FoldbackLenEnd,FoldbackInfoEnd");

            // local topology from arm cluster
            writer.write(",LocTopIdStart,LocTopTypeStart,LocTopTIStart,LocTopIdEnd,LocTopTypeEnd,LocTopTIEnd");

            // gene info
            writer.write(",GeneStart,GeneEnd,Annotations");

            if(mConfig.Output.WriteSvData)
            {
                // extra copy number info
                writer.write(",MinorAPStartPrev,MinorAPStartPost,MinorAPEndPrev,MinorAPEndPost,AFStart,AFEnd");

                // SV table info
                writer.write(",HomologyStart,HomologyEnd,InsertSeq,Imprecise,QualScore");
                writer.write(",RefContextStart,RefContextEnd,InsSeqAlignments");
                writer.write(",Recovered,RepeatClass,RepeatType,AnchorStart,AnchorEnd");
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

                writer.write(String.format("%s,%d,%s,%d,%d",
                        sampleId, var.id(), var.typeStr(), cluster.id(), cluster.getSvCount()));

                writer.write(String.format(",%s,%d,%d,%s,%s,%d,%d,%s",
                        var.chromosome(true), var.position(true), var.orientation(true), asStr(var.arm(true)),
                        var.chromosome(false), var.position(false), var.orientation(false), asStr(var.arm(false))));

                writer.write(String.format(",%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f",
                        dbData.adjustedStartCopyNumber(), dbData.adjustedStartCopyNumberChange(),
                        dbData.adjustedEndCopyNumber(), dbData.adjustedEndCopyNumberChange(),
                        dbData.junctionCopyNumber(), var.jcnMin(), var.jcnMax()));

                writer.write(String.format(",%s,%s,%s",
                        var.getClusterReason(), cluster.getDesc(), cluster.getResolvedType()));

                writer.write(String.format(",%s,%s,%s,%s",
                        var.isFragileSite(true), var.isFragileSite(false),
                        LineElementType.toString(var.getLineElement(true)), LineElementType.toString(var.getLineElement(false))));

                // linked pair info
                for(int be = SE_START; be <= SE_END; ++be)
                {
                    boolean isStart = isStart(be);
                    final LinkedPair link = var.getLinkedPair(isStart);
                    if(link != null)
                    {
                        writer.write(String.format(",%s,%d",
                                link.first() == var ? link.second().id() : link.first().id(), link.baseLength()));
                    }
                    else
                    {
                        writer.write(",,-1");
                    }
                }

                // assembly info
                writer.write(String.format(",%s,%s",
                        var.getAssemblyData(true), var.getAssemblyData(false)));

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

                writer.write(String.format(",%d,%s,%d,%d",
                        var.getNearestSvDistance(), var.getNearestSvRelation(), dbLenStart, dbLenEnd));

                writer.write(String.format(",%d,%d,%s,%d,%d,%s",
                        var.getFoldbackId(true), var.getFoldbackLength(true), var.getFoldbackInfo(true),
                        var.getFoldbackId(false), var.getFoldbackLength(false), var.getFoldbackInfo(false)));

                for(int be = SE_START; be <= SE_END; ++be)
                {
                    ArmCluster armCluster = be == SE_START ? armClusterStart : armClusterEnd;

                    if(armCluster != null)
                        writer.write(String.format(",%d,%s,%d", armCluster.id(), armCluster.getTypeStr(), armCluster.getTICount()));
                    else
                        writer.write(",-1,,0");
                }

                writer.write(String.format(",%s,%s,%s",
                        var.getGeneInBreakend(true, true), var.getGeneInBreakend(false, true),
                        var.getAnnotations()));

                if(mConfig.Output.WriteSvData)
                {
                    writer.write(String.format(",%.2f,%.2f,%.2f,%.2f,%.2f,%.2f",
                            var.getBreakend(true).minorAlleleJcn(true),
                            var.getBreakend(true).minorAlleleJcn(false),
                            !var.isSglBreakend() ? var.getBreakend(false).minorAlleleJcn(true) : 0,
                            !var.isSglBreakend() ? var.getBreakend(false).minorAlleleJcn(false) : 0,
                            dbData.adjustedStartAF(), dbData.adjustedEndAF()));

                    final String insSeqAlignments = dbData.insertSequenceAlignments().replaceAll(",", ";");

                    writer.write(String.format(",%s,%s,%s,%s,%.0f,%s,%s,%s",
                            dbData.startHomologySequence(), dbData.endHomologySequence(),
                            dbData.insertSequence(), dbData.imprecise(), dbData.qualityScore(),
                            dbData.startRefContext(), dbData.endRefContext(), insSeqAlignments));

                    writer.write(String.format(",%s,%s,%s,%d,%d",
                            dbData.recovered(), dbData.insertSequenceRepeatClass(), dbData.insertSequenceRepeatType(),
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
            String outputFileName = mConfig.OutputDataPath + "LNX_CLUSTERS.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("SampleId,ClusterId,ClusterDesc,ClusterCount,Category,ResolvedType,Synthetic,FullyChained,ChainCount");
            writer.write(",DelCount,DupCount,InsCount,InvCount,BndCount,SglCount,InfCount");
            writer.write(",ClusterReasons,Consistency,IsLINE,Replication,MinJcn,MaxJcn,Foldbacks");
            writer.write(",ArmCount,OriginArms,FragmentArms,ConsistentArms,ComplexArms,Annotations,AlleleValidPerc");

            writer.write(",TotalTIs,AssemblyTIs,ShortTIs,IntTIs,ExtTIs,IntShortTIs,ExtShortTIs,IntTIsCnGain");
            writer.write(",ExtTIsCnGain,OverlapTIs,ChainEndsFace,ChainEndsAway,UnchainedSVs");

            writer.write(",DBs,ShortDBs,TotalDBLength,TotalDeleted,TravDelCount,TravDelLength");
            writer.write(",TotalRange,ChainedLength,ImpliedTIs");

            writer.write(",ArmClusterCount,AcTotalTIs,AcIsolatedBE,AcTIOnly,AcDsb,AcSimpleDup");
            writer.write(",AcSingleFb,AcFbDsb,AcComplexFb,AcComplexLine,AcSameOrient,AcComplexOther");

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

                writer.write(String.format("%s,%d,%s,%d,%s,%s,%s,%s,%d",
                        sampleId, cluster.id(), cluster.getDesc(), clusterSvCount,
                        category, resolvedType, cluster.isSyntheticType(),
                        cluster.isFullyChained(false), cluster.getChains().size()));

                writer.write(String.format(",%d,%d,%d,%d,%d,%d,%d",
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

                writer.write(String.format(",%s,%d,%s,%s,%.1f,%.1f,%.0f",
                        cluster.getClusteringReasons(), cluster.getConsistencyCount(), cluster.hasLinkingLineElements(),
                        cluster.requiresReplication(), cluster.getMinJcn(), cluster.getMaxJcn(), foldbackCount));

                final ClusterMetrics metrics = cluster.getMetrics();

                writer.write(String.format(",%d,%d,%d,%d,%d,%s,%.2f",
                        cluster.getArmCount(), metrics.OriginArms, metrics.FragmentArms, metrics.ConsistentArms,
                        metrics.ComplexArms, cluster.getAnnotations(), metrics.ValidAlleleJcnSegmentPerc));

                long shortTIs = cluster.getLinkedPairs().stream().filter(x -> x.baseLength() <= SHORT_TI_LENGTH).count();

                writer.write(String.format(",%d,%d,%d",
                        cluster.getLinkedPairs().size(), cluster.getAssemblyLinkedPairs().size(), shortTIs));

                final ChainMetrics chainMetrics = cluster.getLinkMetrics();

                writer.write(String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
                        chainMetrics.InternalTIs, chainMetrics.ExternalTIs, chainMetrics.InternalShortTIs, chainMetrics.ExternalShortTIs,
                        chainMetrics.InternalTICnGain, chainMetrics.ExternalTICnGain, chainMetrics.OverlappingTIs,
                        chainMetrics.ChainEndsFace, chainMetrics.ChainEndsAway, cluster.getUnlinkedSVs().size()));

                writer.write(String.format(",%d,%d,%d,%d,%d,%d",
                        metrics.DBCount, metrics.ShortDBCount, metrics.TotalDBLength,
                        metrics.TotalDeleted, metrics.TraversedDelCount, metrics.TraversedDelLength));

                writer.write(String.format(",%d,%d,%d",
                        metrics.TotalRange, metrics.ChainedLength, metrics.ImpliedTICount));

                final int[] armClusterData = getArmClusterData(cluster);
                long armClusterTIs = cluster.getArmClusters().stream().mapToInt(x -> x.getTICount()).sum();

                writer.write(String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
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
        if(!mConfig.Output.WriteLinks || !writeCohortFiles())
            return null;

        try
        {
            String outputFileName = mConfig.OutputDataPath + "LNX_LINKS.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("SampleId,ClusterId,ClusterCount,ResolvedType");
            writer.write(",ChainId,ChainCount,ChainConsistent,LinkReason,LinkIndex,ChainIndex,Jcn,JcnUncertainty");
            writer.write(",IsAssembled,TILength,NextSvDist,NextClusteredSvDist,TraversedSVCount");
            writer.write(",LocationType,OverlapCount,CopyNumberGain");
            writer.write(",Id1,Id2,ChrArm,PosStart,PosEnd,LocTopTypeStart,LocTopTypeEnd,GeneStart,GeneEnd,ExonMatch,IsDM");
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
        if(!mConfig.Output.WriteLinks)
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
                                chainIndexStr = appendStr(chainIndexStr, String.valueOf(j), ';');
                            }
                        }

                        final SvBreakend beStart = pair.getBreakend(true);
                        final SvBreakend beEnd = pair.getBreakend(false);

                        writer.write(String.format("%s,%d,%d,%s",
                                sampleId, cluster.id(), clusterSvCount, cluster.getResolvedType()));

                        writer.write(String.format(",%d,%d,%s,%s,%d,%s,%s,%.3f",
                                chain.id(), chainSvCount, chainConsistent, pair.getLinkReason(), pair.getLinkIndex(),
                                chainIndexStr, formatJcn(chain.jcn()), chain.jcnUncertainty()));

                        writer.write(String.format(",%s,%d,%d,%d,%d,%s,%d,%s",
                                pair.isAssembled(), pair.baseLength(),
                                pair.getNextSvDistance(), pair.getNextClusteredSvDistance(), pair.getTraversedSVCount(),
                                pair.locationType(), pair.overlapCount(), pair.hasCopyNumberGain()));

                        ArmCluster acStart = cluster.findArmCluster(beStart);
                        ArmCluster acEnd = cluster.findArmCluster(beEnd);

                        writer.write(String.format(",%d,%d,%s,%d,%d,%s,%s,%s,%s,%s,%s",
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
