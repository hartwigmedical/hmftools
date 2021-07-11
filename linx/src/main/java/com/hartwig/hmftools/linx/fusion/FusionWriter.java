package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.linx.LinxFusion.context;
import static com.hartwig.hmftools.common.sv.linx.LinxFusion.fusionJcn;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.BreakendTransData;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

public class FusionWriter
{
    private final String mOutputDir;
    private BufferedWriter mFusionWriter;

    public FusionWriter(final String outputDir)
    {
        mOutputDir = outputDir;
        mFusionWriter = null;
    }

    public static void convertBreakendsAndFusions(
            final List<GeneFusion> geneFusions, final List<BreakendTransData> transcripts,
            final List<LinxFusion> fusions, final List<LinxBreakend> breakends)
    {
        int breakendId = 0;
        Map<BreakendTransData,Integer> transIdMap = Maps.newHashMap();

        for(final BreakendTransData transcript : transcripts)
        {
            transIdMap.put(transcript, breakendId);

            final BreakendGeneData gene = transcript.gene();

            breakends.add(ImmutableLinxBreakend.builder()
                    .id(breakendId++)
                    .svId(transcript.gene().id())
                    .isStart(transcript.gene().isStart())
                    .type(gene.type().toString())
                    .gene(transcript.geneName())
                    .transcriptId(transcript.transName())
                    .canonical(transcript.isCanonical())
                    .geneOrientation(transcript.isUpstream() ? "Upstream" : "Downstream")
                    .disruptive(transcript.isDisruptive())
                    .reportedDisruption(transcript.reportableDisruption())
                    .junctionCopyNumber(gene.jcn())
                    .undisruptedCopyNumber(transcript.undisruptedCopyNumber())
                    .regionType(transcript.regionType().toString())
                    .codingContext(transcript.codingType().toString())
                    .biotype(transcript.bioType())
                    .exonicBasePhase(transcript.Phase)
                    .nextSpliceExonRank(transcript.nextSpliceExonRank())
                    .nextSpliceExonPhase(transcript.Phase)
                    .nextSpliceDistance(transcript.isUpstream() ? transcript.prevSpliceAcceptorDistance() : transcript.nextSpliceAcceptorDistance())
                    .totalExonCount(transcript.TransData.exons().size())
                    .chromosome(gene.chromosome())
                    .orientation(gene.orientation())
                    .strand(gene.Strand)
                    .chrBand(gene.karyotypeBand())
                    .exonUp(transcript.ExonUpstream)
                    .exonDown(transcript.ExonDownstream)
                    .build());
        }

        for(final GeneFusion geneFusion : geneFusions)
        {
            int upBreakendId = transIdMap.get(geneFusion.upstreamTrans());
            int downBreakendId = transIdMap.get(geneFusion.downstreamTrans());

            fusions.add(ImmutableLinxFusion.builder()
                    .fivePrimeBreakendId(upBreakendId)
                    .threePrimeBreakendId(downBreakendId)
                    .name(geneFusion.name())
                    .reported(geneFusion.reportable())
                    .reportedType(geneFusion.knownTypeStr())
                    .phased(geneFusion.phaseType())
                    .likelihood(geneFusion.likelihoodType())
                    .chainLength(geneFusion.getChainLength())
                    .chainLinks(geneFusion.getChainLinks())
                    .chainTerminated(geneFusion.isTerminated())
                    .domainsKept(geneFusion.downstreamTrans().getProteinFeaturesKept())
                    .domainsLost(geneFusion.downstreamTrans().getProteinFeaturesLost())
                    .skippedExonsUp(geneFusion.getExonsSkipped(true))
                    .skippedExonsDown(geneFusion.getExonsSkipped(false))
                    .fusedExonUp(geneFusion.getFusedExon(true))
                    .fusedExonDown(geneFusion.getFusedExon(false))
                    .geneStart(geneFusion.geneName(FS_UP))
                    .geneTranscriptStart(geneFusion.upstreamTrans().transName())
                    .geneContextStart(context(geneFusion.upstreamTrans(), geneFusion.getFusedExon(true)))
                    .geneEnd(geneFusion.geneName(FS_DOWN))
                    .geneTranscriptEnd(geneFusion.downstreamTrans().transName())
                    .geneContextEnd(context(geneFusion.downstreamTrans(), geneFusion.getFusedExon(false)))
                    .junctionCopyNumber(fusionJcn(geneFusion.upstreamTrans().gene().jcn(), geneFusion.downstreamTrans().gene().jcn()))
                    .build());
        }
    }

    public void writeSampleData(final String sampleId, final List<LinxFusion> fusions, final List<LinxBreakend> breakends)
    {
        try
        {
            // write flat files for database loading
            final String breakendsFile = LinxBreakend.generateFilename(mOutputDir, sampleId);
            LinxBreakend.write(breakendsFile, breakends);

            final String fusionsFile = LinxFusion.generateFilename(mOutputDir, sampleId);
            LinxFusion.write(fusionsFile, fusions);
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to write fusions file: {}", e.toString());
        }
    }

    public void initialiseOutputFiles()
    {
        // initialise the integrated, verbose fusion and breakend output file
        try
        {
            if(mFusionWriter == null)
            {
                mFusionWriter = createBufferedWriter(mOutputDir + "LNX_FUSIONS.csv", false);

                mFusionWriter.write("SampleId,Reportable,ReportableReason,KnownType,Phased,KnownExons,ClusterId,ClusterCount,ResolvedType");

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    String upDown = se == SE_START ? "Up" : "Down";

                    String fieldsStr = ",SvId" + upDown;
                    fieldsStr += ",Chr" + upDown;
                    fieldsStr += ",Pos" + upDown;
                    fieldsStr += ",Orient" + upDown;
                    fieldsStr += ",Type" + upDown;
                    fieldsStr += ",Ploidy" + upDown;
                    fieldsStr += ",GeneId" + upDown;
                    fieldsStr += ",GeneName" + upDown;
                    fieldsStr += ",Transcript" + upDown;
                    fieldsStr += ",Strand" + upDown;
                    fieldsStr += ",RegionType" + upDown;
                    fieldsStr += ",CodingType" + upDown;
                    fieldsStr += ",BreakendExon" + upDown;
                    fieldsStr += ",FusedExon" + upDown;
                    fieldsStr += ",ExonsSkipped" + upDown;
                    fieldsStr += ",Phase" + upDown;
                    fieldsStr += ",ExonMax" + upDown;
                    fieldsStr += ",Disruptive" + upDown;
                    fieldsStr += ",CodingBases" + upDown;
                    fieldsStr += ",TotalCoding" + upDown;
                    fieldsStr += ",CodingStart" + upDown;
                    fieldsStr += ",CodingEnd" + upDown;
                    fieldsStr += ",TransStart" + upDown;
                    fieldsStr += ",TransEnd" + upDown;
                    fieldsStr += ",DistancePrev" + upDown;
                    fieldsStr += ",Canonical" + upDown;
                    fieldsStr += ",Biotype" + upDown;
                    mFusionWriter.write(fieldsStr);
                }

                mFusionWriter.write(",ProteinsKept,ProteinsLost,PriorityScore,FusionId,ChainTerminatedUp,ChainTerminatedDown,ChainInfo");
                mFusionWriter.newLine();
            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing fusions: {}", e.toString());
        }
    }

    public void writeVerboseFusionData(final GeneFusion fusion, final String sampleId)
    {
        if(mFusionWriter == null)
            return;

        try
        {
            BufferedWriter writer = mFusionWriter;

            final FusionAnnotations annotations = fusion.getAnnotations();

            if(annotations == null)
            {
                LNX_LOGGER.error("fusion({}) annotations not set", fusion.name());
                return;
            }

            writer.write(String.format("%s,%s,%s,%s",
                    sampleId, fusion.reportable(), fusion.reportableReason(), fusion.knownTypeStr()));

            writer.write(String.format(",%s,%s,%d,%d,%s",
                    fusion.phaseType(), fusion.knownExons(),
                    annotations.clusterId(), annotations.clusterCount(), annotations.resolvedType()));

            // write upstream SV, transcript and exon info
            for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
            {
                boolean isUpstream = (fs == FS_UP);
                final BreakendTransData trans = fusion.transcripts()[fs];
                final BreakendGeneData gene = trans.gene();

                writer.write(String.format(",%d,%s,%d,%d,%s,%.6f",
                        gene.id(), gene.chromosome(), gene.position(), gene.orientation(),
                        gene.type(), gene.jcn()));

                writer.write(String.format(",%s,%s,%s,%d,%s,%s",
                        gene.StableId, fusion.geneName(fs), trans.transName(),
                        gene.Strand, trans.regionType(), trans.codingType()));

                writer.write(String.format(",%d,%d,%d,%d,%d,%s",
                        isUpstream ? trans.ExonUpstream : trans.ExonDownstream,
                        fusion.getFusedExon(isUpstream), fusion.getExonsSkipped()[fs],
                        trans.Phase, trans.exonCount(), trans.isDisruptive()));

                writer.write(String.format(",%d,%d,%d,%d,%d,%d,%d,%s,%s",
                        trans.CodingBases, trans.TotalCodingBases,
                        trans.codingStart(), trans.codingEnd(), trans.transStart(), trans.transEnd(),
                        trans.prevSpliceAcceptorDistance(), trans.isCanonical(), trans.bioType()));
            }

            writer.write(String.format(",%s,%s,%f,%d",
                        fusion.downstreamTrans().getProteinFeaturesKept(), fusion.downstreamTrans().getProteinFeaturesLost(),
                        fusion.priority(), fusion.id()));

            String chainInfo = String.format(",%s,%s", annotations.terminatedUp(), annotations.terminatedDown());

            if(annotations.chainInfo() != null)
            {
                chainInfo += String.format(",%d;%d;%d;%s;%s",
                        annotations.chainInfo().chainId(), annotations.chainInfo().chainLinks(), annotations.chainInfo().chainLength(),
                        annotations.chainInfo().validTraversal(), annotations.chainInfo().traversalAssembled());
            }
            else
            {
                chainInfo += ",-1;0;0;true;false";
            }

            writer.write(String.format("%s", chainInfo));

            writer.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing fusions: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mFusionWriter);
    }
}
