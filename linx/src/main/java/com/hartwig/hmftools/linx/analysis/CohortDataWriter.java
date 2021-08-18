package com.hartwig.hmftools.linx.analysis;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
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
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.asStr;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_COMPLEX_FOLDBACK;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_COMPLEX_LINE;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_COMPLEX_OTHER;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_DSB;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_FOLDBACK;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_FOLDBACK_DSB;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_ISOLATED_BE;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_SAME_ORIENT;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_SIMPLE_DUP;
import static com.hartwig.hmftools.linx.types.ArmCluster.ARM_CL_TI_ONLY;
import static com.hartwig.hmftools.linx.types.ArmCluster.getArmClusterData;
import static com.hartwig.hmftools.linx.types.LinxConstants.NO_DB_MARKER;
import static com.hartwig.hmftools.linx.types.LinxConstants.SHORT_TI_LENGTH;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.annotators.LineElementType;
import com.hartwig.hmftools.linx.chaining.ChainMetrics;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.germline.GermlinePonCache;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.ArmCluster;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisDataWriter;

public class CohortDataWriter
{
    private final LinxConfig mConfig;

    private final BufferedWriter mSvFileWriter;
    private final BufferedWriter mClusterFileWriter;
    private final BufferedWriter mLinksFileWriter;
    private final VisDataWriter mVisWriter;
    private final GermlinePonCache mGermlinePonCache;

    public CohortDataWriter(final LinxConfig config)
    {
        mConfig = config;
        mVisWriter = new VisDataWriter(config.OutputDataPath, config.Output.WriteVisualisationData, config.hasMultipleSamples());

        mGermlinePonCache = config.IsGermline ? new GermlinePonCache(config.CmdLineArgs) : null;

        mSvFileWriter = createSvDataFile();
        mClusterFileWriter = createClusterFile();
        mLinksFileWriter = createLinksFile();
    }

    public final VisDataWriter getVisWriter() { return mVisWriter; }

    public void close()
    {
        closeBufferedWriter(mSvFileWriter);
        closeBufferedWriter(mClusterFileWriter);
        closeBufferedWriter(mLinksFileWriter);
        mVisWriter.close();
    }

    public boolean writeCohortFiles()
    {
        return mConfig.hasMultipleSamples() || mConfig.Output.WriteCohortFiles;
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

            // gene & replication info
            writer.write(",GeneStart,GeneEnd,RepOriginStart,RepOriginEnd,Annotations");

            if(mConfig.Output.WriteSvData)
            {
                // extra copy number info
                writer.write(",MinorAPStartPrev,MinorAPStartPost,MinorAPEndPrev,MinorAPEndPost,AFStart,AFEnd");

                // SV table info
                writer.write(",HomologyStart,HomologyEnd,InsertSeq,Imprecise,QualScore");
                writer.write(",RefContextStart,RefContextEnd,InsSeqAlignments");
                writer.write(",Recovered,RepeatClass,RepeatType,AnchorStart,AnchorEnd");
            }

            if(mConfig.IsGermline)
            {
                writer.write(",Filter,PonCount");
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

    public void writeSvData(final String sampleId, final List<SvVarData> svDataList)
    {
        if(mSvFileWriter == null)
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

                mSvFileWriter.write(String.format("%s,%d,%s,%d,%d",
                        sampleId, var.id(), var.typeStr(), cluster.id(), cluster.getSvCount()));

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
                        LineElementType.toString(var.getLineElement(true)), LineElementType.toString(var.getLineElement(false))));

                // linked pair info
                for (int be = SE_START; be <= SE_END; ++be)
                {
                    boolean isStart = isStart(be);
                    final LinkedPair link = var.getLinkedPair(isStart);
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

                mSvFileWriter.write(String.format(",%d,%s,%d,%d",
                        var.getNearestSvDistance(), var.getNearestSvRelation(), dbLenStart, dbLenEnd));

                mSvFileWriter.write(String.format(",%d,%d,%s,%d,%d,%s",
                        var.getFoldbackId(true), var.getFoldbackLength(true), var.getFoldbackInfo(true),
                        var.getFoldbackId(false), var.getFoldbackLength(false), var.getFoldbackInfo(false)));

                for (int be = SE_START; be <= SE_END; ++be)
                {
                    ArmCluster armCluster = be == SE_START ? armClusterStart : armClusterEnd;

                    if (armCluster != null)
                        mSvFileWriter.write(String.format(",%d,%s,%d", armCluster.id(), armCluster.getTypeStr(), armCluster.getTICount()));
                    else
                        mSvFileWriter.write(",-1,,0");
                }

                mSvFileWriter.write(String.format(",%s,%s,%.4f,%.4f,%s",
                        var.getGeneInBreakend(true, true), var.getGeneInBreakend(false, true),
                        var.getReplicationOrigin(true), var.getReplicationOrigin(false), var.getAnnotations()));

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
                }

                if(mConfig.IsGermline)
                {
                    int ponCount = var.getSvData().filter().equals(PON_FILTER_PON) ? mGermlinePonCache.getPonCount(var) : 0;
                    mSvFileWriter.write(String.format(",%s,%d", var.getSvData().filter(), ponCount));
                }

                mSvFileWriter.newLine();
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

            if(mConfig.IndelAnnotation)
                writer.write(",IndelCount,IndelProb");

            writer.newLine();
            return writer;
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing cluster-data to outputFile: {}", e.toString());
            return null;
        }
    }

    public void writeClusterData(final String sampleId, final List<SvCluster> clusters)
    {
        if(mClusterFileWriter == null)
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

                mClusterFileWriter.write(String.format("%s,%d,%s,%d,%s,%s,%s,%s,%d",
                        sampleId, cluster.id(), cluster.getDesc(), clusterSvCount,
                        category, resolvedType, cluster.isSyntheticType(),
                        cluster.isFullyChained(false), cluster.getChains().size()));

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
                        cluster.getArmCount(), metrics.OriginArms, metrics.FragmentArms, metrics.ConsistentArms,
                        metrics.ComplexArms, cluster.getAnnotations(), metrics.ValidAlleleJcnSegmentPerc));

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

                if(mConfig.IndelAnnotation)
                {
                    mClusterFileWriter.write(String.format(",%d,%f", metrics.IndelCount, metrics.IndelProbability));
                }

                mClusterFileWriter.newLine();

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

    public void writeLinksData(final String sampleId, final List<SvCluster> clusters)
    {
        if(!mConfig.Output.WriteLinks || mLinksFileWriter == null)
            return;

        try
        {
            for(final SvCluster cluster : clusters)
            {
                int clusterSvCount = cluster.getSvCount();

                List<SvChain> chains = cluster.getChains();

                for (final SvChain chain : chains)
                {
                    int chainSvCount = chain.getSvCount();
                    boolean chainConsistent = chain.isConsistent();
                    boolean isDoubleMinute = chain.isDoubleMinute();

                    List<LinkedPair> uniquePairs = Lists.newArrayList();
                    final List<LinkedPair> chainLinks = chain.getLinkedPairs();

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

                        mLinksFileWriter.write(String.format("%s,%d,%d,%s",
                                sampleId, cluster.id(), clusterSvCount, cluster.getResolvedType()));

                        mLinksFileWriter.write(String.format(",%d,%d,%s,%s,%d,%s,%s,%.3f",
                                chain.id(), chainSvCount, chainConsistent, pair.getLinkReason(), pair.getLinkIndex(),
                                chainIndexStr, formatJcn(chain.jcn()), chain.jcnUncertainty()));

                        mLinksFileWriter.write(String.format(",%s,%d,%d,%d,%d,%s,%d,%s",
                                pair.isAssembled(), pair.length(),
                                pair.getNextSvDistance(), pair.getNextClusteredSvDistance(), pair.getTraversedSVCount(),
                                pair.locationType(), pair.overlapCount(), pair.hasCopyNumberGain()));

                        ArmCluster acStart = cluster.findArmCluster(beStart);
                        ArmCluster acEnd = cluster.findArmCluster(beEnd);

                        mLinksFileWriter.write(String.format(",%d,%d,%s,%d,%d,%s,%s,%s,%s,%s,%s",
                                beStart.getSV().id(), beEnd.getSV().id(),
                                beStart.getChrArm(), beStart.position(), beEnd.position(),
                                acStart != null ? acStart.getTypeStr() : "", acEnd != null ? acEnd.getTypeStr() : "",
                                beStart.getSV().getGeneInBreakend(beStart.usesStart(), false),
                                beEnd.getSV().getGeneInBreakend(beEnd.usesStart(), false),
                                pair.getExonMatchData(), isDoubleMinute));

                        mLinksFileWriter.newLine();
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
