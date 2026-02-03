package com.hartwig.hmftools.linx;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StartEndIterator.isStart;
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
import static com.hartwig.hmftools.linx.types.ArmCluster.getArmClusterData;
import static com.hartwig.hmftools.linx.types.LinxConstants.NO_DB_MARKER;
import static com.hartwig.hmftools.linx.types.LinxConstants.SHORT_TI_LENGTH;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

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
            String outputFileName = cohortDataFilename(mConfig.OutputDataPath, mConfig.IsGermline ? "GERMLINE_SVS" : "SVS");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            // definitional fields
            sj.add("SampleId").add("Id").add("Type").add("ClusterId").add("ClusterCount");
            sj.add("ChrStart").add("PosStart").add("OrientStart").add("ArmStart").add("ChrEnd").add("PosEnd").add("OrientEnd").add("ArmEnd");

            if(mConfig.isSomatic())
            {
                // position and copy number
                sj.add("CNStart").add("CNChgStart").add("CNEnd").add("CNChgEnd").add("Jcn").add("JcnMin").add("JcnMax");
            }

            // cluster info
            sj.add("ClusterReason").add("ClusterDesc").add("ResolvedType");

            sj.add("FSStart").add("FSEnd").add("LEStart").add("LEEnd");

            // linked pair info
            sj.add("LnkSvStart").add("LnkLenStart").add("LnkSvEnd").add("LnkLenEnd").add("AsmbStart").add("AsmbEnd");

            // chain info
            sj.add("ChainId").add("ChainCount").add("ChainIndex");

            // proximity info and other link info
            sj.add("NearestLen").add("NearestType").add("DBLenStart").add("DBLenEnd");

            if(mConfig.isSomatic())
            {
                // proximity info and other link info
                sj.add("FoldbackLnkStart").add("FoldbackLenStart").add("FoldbackInfoStart");
                sj.add("FoldbackLnkEnd").add("FoldbackLenEnd").add("FoldbackInfoEnd");

                // local topology from arm cluster
                sj.add("LocTopIdStart").add("LocTopTypeStart").add("LocTopTIStart").add("LocTopIdEnd").add("LocTopTypeEnd").add("LocTopTIEnd");
            }

            // gene info
            sj.add("GeneStart").add("GeneEnd");

            if(mConfig.Output.writeSvData())
            {
                // extra copy number info
                if(mConfig.isSomatic())
                    sj.add("MinorAPStartPrev").add("MinorAPStartPost").add("MinorAPEndPrev").add("MinorAPEndPost");

                // SV table info
                sj.add("HomologyStart").add("HomologyEnd").add("InsertSeq").add("QualScore").add("AFStart").add("AFEnd");
                sj.add("InsSeqAlignments").add("RepeatClass").add("RepeatType").add("AnchorStart").add("AnchorEnd");
            }

            writer.write(sj.toString());

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

                // if(mConfig.IsGermline && var.getGenesList(true).isEmpty() && var.getGenesList(false).isEmpty())
                //    continue;

                final StructuralVariantData dbData = var.getSvData();

                final ArmCluster armClusterStart = cluster.findArmCluster(var.getBreakend(true));

                final ArmCluster armClusterEnd = !var.isSglBreakend() ? cluster.findArmCluster(var.getBreakend(false)) : null;

                StringJoiner sj = new StringJoiner(TSV_DELIM);

                sj.add(sampleId);
                sj.add(String.valueOf(var.id()));
                sj.add(var.typeStr());
                sj.add(String.valueOf(cluster.id()));
                sj.add(String.valueOf(cluster.getSvCount()));

                sj.add(var.chromosome(true));
                sj.add(String.valueOf(var.position(true)));
                sj.add(String.valueOf(var.orientation(true)));
                sj.add(var.arm(true).toString());
                sj.add(var.chromosome(false));
                sj.add(String.valueOf(var.position(false)));
                sj.add(String.valueOf(var.orientation(false)));
                sj.add(var.arm(false).toString());

                if(mConfig.isSomatic())
                {
                    sj.add(format("%.2f", dbData.adjustedStartCopyNumber()));
                    sj.add(format("%.2f", dbData.adjustedStartCopyNumberChange()));
                    sj.add(format("%.2f", dbData.adjustedEndCopyNumber()));
                    sj.add(format("%.2f", dbData.adjustedEndCopyNumberChange()));
                    sj.add(format("%.2f", dbData.junctionCopyNumber()));
                    sj.add(format("%.2f", var.jcnMin()));
                    sj.add(format("%.2f", var.jcnMax()));
                }

                sj.add(var.getClusterReason());
                sj.add(cluster.getDesc());
                sj.add(String.valueOf(cluster.getResolvedType()));

                sj.add(String.valueOf(var.isFragileSite(true)));
                sj.add(String.valueOf(var.isFragileSite(false)));
                sj.add(LineElementType.toString(var.getLineElement(true)));
                sj.add(LineElementType.toString(var.getLineElement(false)));

                // linked pair info
                for(int be = SE_START; be <= SE_END; ++be)
                {
                    boolean isStart = isStart(be);
                    final LinkedPair link = var.getLinkedPair(isStart);
                    if(link != null)
                    {
                        sj.add(String.valueOf(link.first() == var ? link.second().id() : link.first().id()));
                        sj.add(String.valueOf(link.baseLength()));
                    }
                    else
                    {
                        sj.add("");
                        sj.add("-1");
                    }
                }

                // assembly info
                sj.add(var.assemblyInfoStr(true));
                sj.add(var.assemblyInfoStr(false));

                // chain info
                final SvChain chain = cluster.findChain(var);
                String chainStr = "";

                if(chain != null)
                {
                    chainStr = format("%d\t%d\t%s", chain.id(), chain.getSvCount(), chain.getSvIndices(var));
                }
                else
                {
                    chainStr = format("%d\t0\t", cluster.getChainId(var));
                }

                sj.add(chainStr);

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

                sj.add(String.valueOf(var.getNearestSvDistance()));
                sj.add(var.getNearestSvRelation());
                sj.add(String.valueOf(dbLenStart));
                sj.add(String.valueOf(dbLenEnd));

                if(mConfig.isSomatic())
                {
                    sj.add(String.valueOf(var.getFoldbackId(true)));
                    sj.add(String.valueOf(var.getFoldbackLength(true)));
                    sj.add(var.getFoldbackInfo(true));
                    sj.add(String.valueOf(var.getFoldbackId(false)));
                    sj.add(String.valueOf(var.getFoldbackLength(false)));
                    sj.add(var.getFoldbackInfo(false));

                    for(int be = SE_START; be <= SE_END; ++be)
                    {
                        ArmCluster armCluster = be == SE_START ? armClusterStart : armClusterEnd;

                        if(armCluster != null)
                            sj.add(format("%d\t%s\t%d", armCluster.id(), armCluster.getTypeStr(), armCluster.getTICount()));
                        else
                            sj.add("-1\t\t0");
                    }
                }

                sj.add(var.getGeneInBreakend(true, false, true));
                sj.add(var.getGeneInBreakend(false, false, true));

                if(mConfig.Output.writeSvData())
                {
                    if(mConfig.isSomatic())
                    {
                        sj.add(format("%.2f", var.getBreakend(true).minorAlleleJcn(true)));
                        sj.add(format("%.2f", var.getBreakend(true).minorAlleleJcn(false)));
                        sj.add(format("%.2f", !var.isSglBreakend() ? var.getBreakend(false).minorAlleleJcn(true) : 0));
                        sj.add(format("%.2f", !var.isSglBreakend() ? var.getBreakend(false).minorAlleleJcn(false) : 0));
                    }

                    final String insSeqAlignments = dbData.insertSequenceAlignments().replaceAll(",", ";");

                    sj.add(dbData.startHomologySequence());
                    sj.add(dbData.endHomologySequence());
                    sj.add(dbData.insertSequence());
                    sj.add(format("%.0f", dbData.qualityScore()));
                    sj.add(format("%.2f", dbData.adjustedStartAF()));
                    sj.add(format("%.2f", dbData.adjustedEndAF()));

                    sj.add(insSeqAlignments);
                    sj.add(dbData.insertSequenceRepeatClass());
                    sj.add(dbData.insertSequenceRepeatType());
                    sj.add(String.valueOf(dbData.startAnchoringSupportDistance()));
                    sj.add(String.valueOf(dbData.endAnchoringSupportDistance()));
                }

                writer.write(sj.toString());
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
            String outputFileName = cohortDataFilename(mConfig.OutputDataPath, mConfig.IsGermline ? "GERMLINE_CLUSTERS" : "CLUSTERS");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("SampleId").add("ClusterId").add("ClusterDesc").add("ClusterCount").add("Category").add("ResolvedType").add("Synthetic");
            sj.add("FullyChained").add("ChainCount").add("DelCount").add("DupCount").add("InsCount").add("InvCount").add("BndCount");
            sj.add("SglCount").add("InfCount").add("ClusterReasons").add("Consistency").add("IsLINE");

            if(mConfig.isSomatic())
            {
                sj.add("Replication").add("MinJcn").add("MaxJcn").add("Foldbacks").add("ArmCount").add("OriginArms");
                sj.add("FragmentArms").add("ConsistentArms").add("ComplexArms").add("Annotations").add("AlleleValidPerc");

                sj.add("TotalTIs").add("AssemblyTIs").add("ShortTIs").add("IntTIs").add("ExtTIs").add("IntShortTIs").add("ExtShortTIs");
                sj.add("IntTIsCnGain").add("ExtTIsCnGain").add("OverlapTIs").add("ChainEndsFace").add("ChainEndsAway").add("UnchainedSVs");
            }

            sj.add("DBs").add("ShortDBs").add("TotalDBLength").add("TotalDeleted").add("TravDelCount").add("TravDelLength");
            sj.add("TotalRange").add("ChainedLength").add("ImpliedTIs");

            sj.add("ArmClusterCount").add("AcTotalTIs").add("AcIsolatedBE").add("AcTIOnly").add("AcDsb").add("AcSimpleDup");
            sj.add("AcSingleFb").add("AcFbDsb").add("AcComplexFb").add("AcComplexLine").add("AcSameOrient").add("AcComplexOther");

            writer.write(sj.toString());
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

                StringJoiner sj = new StringJoiner(TSV_DELIM);

                sj.add(sampleId);
                sj.add(String.valueOf(cluster.id()));
                sj.add(cluster.getDesc());
                sj.add(String.valueOf(clusterSvCount));
                sj.add(category);
                sj.add(String.valueOf(resolvedType));
                sj.add(String.valueOf(cluster.isSyntheticType()));
                sj.add(String.valueOf(cluster.isFullyChained(false)));
                sj.add(String.valueOf(cluster.getChains().size()));

                sj.add(String.valueOf(cluster.getTypeCount(DEL)));
                sj.add(String.valueOf(cluster.getTypeCount(DUP)));
                sj.add(String.valueOf(cluster.getTypeCount(INS)));
                sj.add(String.valueOf(cluster.getTypeCount(INV)));
                sj.add(String.valueOf(cluster.getTypeCount(BND)));
                sj.add(String.valueOf(cluster.getTypeCount(SGL)));
                sj.add(String.valueOf(cluster.getTypeCount(INF)));

                double foldbackCount = 0;

                for(final SvVarData var : cluster.getFoldbacks())
                {
                    // avoid double-count chained foldbacks
                    if(var.getFoldbackBreakend(true) != null)
                        foldbackCount += 0.5;
                    if(var.getFoldbackBreakend(false) != null)
                        foldbackCount += 0.5;
                }

                sj.add(cluster.getClusteringReasons());
                sj.add(String.valueOf(cluster.getConsistencyCount()));
                sj.add(String.valueOf(cluster.hasLinkingLineElements()));

                final ClusterMetrics metrics = cluster.getMetrics();

                if(mConfig.isSomatic())
                {
                    sj.add(String.valueOf(cluster.requiresReplication()));
                    sj.add(format("%.2f", cluster.getMinJcn()));
                    sj.add(format("%.2f", cluster.getMaxJcn()));
                    sj.add(String.valueOf(foldbackCount));

                    sj.add(String.valueOf(cluster.getArmCount()));
                    sj.add(String.valueOf(metrics.OriginArms));
                    sj.add(String.valueOf(metrics.FragmentArms));
                    sj.add(String.valueOf(metrics.ConsistentArms));
                    sj.add(String.valueOf(metrics.ComplexArms));
                    sj.add(cluster.getAnnotations());
                    sj.add(format("%.2f", metrics.ValidAlleleJcnSegmentPerc));

                    long shortTIs = cluster.getLinkedPairs().stream().filter(x -> x.baseLength() <= SHORT_TI_LENGTH).count();

                    sj.add(String.valueOf(cluster.getLinkedPairs().size()));
                    sj.add(String.valueOf(cluster.getAssemblyLinkedPairs().size()));
                    sj.add(String.valueOf(shortTIs));

                    final ChainMetrics chainMetrics = cluster.getLinkMetrics();

                    sj.add(String.valueOf(chainMetrics.InternalTIs));
                    sj.add(String.valueOf(chainMetrics.ExternalTIs));
                    sj.add(String.valueOf(chainMetrics.InternalShortTIs));
                    sj.add(String.valueOf(chainMetrics.ExternalShortTIs));
                    sj.add(String.valueOf(chainMetrics.InternalTICnGain));
                    sj.add(String.valueOf(chainMetrics.ExternalTICnGain));
                    sj.add(String.valueOf(chainMetrics.OverlappingTIs));
                    sj.add(String.valueOf(chainMetrics.ChainEndsFace));
                    sj.add(String.valueOf(chainMetrics.ChainEndsAway));
                    sj.add(String.valueOf(cluster.getUnlinkedSVs().size()));
                }

                sj.add(String.valueOf(metrics.DBCount));
                sj.add(String.valueOf(metrics.ShortDBCount));
                sj.add(String.valueOf(metrics.TotalDBLength));
                sj.add(String.valueOf(metrics.TotalDeleted));
                sj.add(String.valueOf(metrics.TraversedDelCount));
                sj.add(String.valueOf(metrics.TraversedDelLength));

                sj.add(String.valueOf(metrics.TotalRange));
                sj.add(String.valueOf(metrics.ChainedLength));
                sj.add(String.valueOf(metrics.ImpliedTICount));

                final int[] armClusterData = getArmClusterData(cluster);
                long armClusterTIs = cluster.getArmClusters().stream().mapToInt(x -> x.getTICount()).sum();

                sj.add(String.valueOf(cluster.getArmClusters().size()));
                sj.add(String.valueOf(armClusterTIs));
                sj.add(String.valueOf(armClusterData[ArmClusterType.ISOLATED_BE.ordinal()]));
                sj.add(String.valueOf(armClusterData[ArmClusterType.TI_ONLY.ordinal()]));
                sj.add(String.valueOf(armClusterData[ArmClusterType.DSB.ordinal()]));
                sj.add(String.valueOf(armClusterData[ArmClusterType.SIMPLE_DUP.ordinal()]));
                sj.add(String.valueOf(armClusterData[ArmClusterType.FOLDBACK.ordinal()]));
                sj.add(String.valueOf(armClusterData[ArmClusterType.FOLDBACK_DSB.ordinal()]));
                sj.add(String.valueOf(armClusterData[ArmClusterType.COMPLEX_FOLDBACK.ordinal()]));
                sj.add(String.valueOf(armClusterData[ArmClusterType.COMPLEX_LINE.ordinal()]));
                sj.add(String.valueOf(armClusterData[ArmClusterType.SAME_ORIENT.ordinal()]));
                sj.add(String.valueOf(armClusterData[ArmClusterType.COMPLEX_OTHER.ordinal()]));

                writer.write(sj.toString());
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
            String outputFileName = cohortDataFilename(mConfig.OutputDataPath, mConfig.IsGermline ? "GERMLINE_LINKS" : "LINKS");;

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add("SampleId").add("ClusterId").add("ClusterCount").add("ResolvedType");
            sj.add("ChainId").add("ChainCount").add("ChainConsistent").add("LinkReason").add("LinkIndex").add("ChainIndex").add("Jcn").add("JcnUncertainty");
            sj.add("IsAssembled").add("TILength").add("NextSvDist").add("NextClusteredSvDist");

            if(mConfig.isSomatic())
            {
                sj.add("TraversedSVCount").add("LocationType").add("OverlapCount").add("CopyNumberGain");
            }

            sj.add("Id1").add("Id2").add("ChrArm").add("PosStart").add("PosEnd");
            sj.add("LocTopTypeStart").add("LocTopTypeEnd").add("GeneStart").add("GeneEnd").add("ExonMatch");

            if(mConfig.isSomatic())
            {
                sj.add("IsDM");
            }

            writer.write(sj.toString());
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

                        StringJoiner sj = new StringJoiner(TSV_DELIM);

                        sj.add(sampleId);
                        sj.add(String.valueOf(cluster.id()));
                        sj.add(String.valueOf(clusterSvCount));
                        sj.add(String.valueOf(cluster.getResolvedType()));

                        sj.add(String.valueOf(chain.id()));
                        sj.add(String.valueOf(chainSvCount));
                        sj.add(String.valueOf(chainConsistent));
                        sj.add(pair.getLinkReason());
                        sj.add(String.valueOf(pair.getLinkIndex()));
                        sj.add(chainIndexStr);
                        sj.add(formatJcn(chain.jcn()));
                        sj.add(format("%.3f", chain.jcnUncertainty()));

                        sj.add(String.valueOf(pair.isAssembled()));
                        sj.add(String.valueOf(pair.baseLength()));
                        sj.add(String.valueOf(pair.getNextSvDistance()));
                        sj.add(String.valueOf(pair.getNextClusteredSvDistance()));

                        if(mConfig.isSomatic())
                        {
                            sj.add(String.valueOf(pair.getTraversedSVCount()));
                            sj.add(pair.locationType());
                            sj.add(String.valueOf(pair.overlapCount()));
                            sj.add(String.valueOf(pair.hasCopyNumberGain()));
                        }

                        ArmCluster acStart = cluster.findArmCluster(beStart);
                        ArmCluster acEnd = cluster.findArmCluster(beEnd);

                        sj.add(String.valueOf(beStart.getSV().id()));
                        sj.add(String.valueOf(beEnd.getSV().id()));
                        sj.add(beStart.getChrArm());
                        sj.add(String.valueOf(beStart.position()));
                        sj.add(String.valueOf(beEnd.position()));
                        sj.add(acStart != null ? acStart.getTypeStr() : "");
                        sj.add(acEnd != null ? acEnd.getTypeStr() : "");
                        sj.add(beStart.getSV().getGeneInBreakend(beStart.usesStart(), false));
                        sj.add(beEnd.getSV().getGeneInBreakend(beEnd.usesStart(), false));

                        if(mConfig.IsGermline)
                        {

                        }


                        sj.add(pair.getExonMatchData());

                        if(mConfig.isSomatic())
                            sj.add(String.valueOf(isDoubleMinute));

                        writer.write(sj.toString());
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
