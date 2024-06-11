package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.linx.LinxBreakend.BREAKEND_ORIENTATION_DOWNSTREAM;
import static com.hartwig.hmftools.common.linx.LinxBreakend.BREAKEND_ORIENTATION_UPSTREAM;
import static com.hartwig.hmftools.common.linx.LinxFusion.context;
import static com.hartwig.hmftools.common.linx.LinxFusion.fusionJcn;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.linx.CohortDataWriter;
import com.hartwig.hmftools.linx.CohortFileInterface;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.linx.gene.BreakendTransData;

public class FusionWriter implements CohortFileInterface
{
    private final String mOutputDir;
    private final CohortDataWriter mCohortDataWriter;

    public FusionWriter(final String outputDir, final CohortDataWriter cohortDataWriter)
    {
        mOutputDir = outputDir;
        mCohortDataWriter = cohortDataWriter;
    }

    public static void convertBreakendsAndFusions(
            final List<GeneFusion> geneFusions, final List<BreakendTransData> transcripts,
            final List<LinxFusion> fusions, final List<LinxBreakend> breakends)
    {
        int breakendId = 0;
        Map<BreakendTransData, Integer> transIdMap = Maps.newHashMap();

        for(final BreakendTransData transcript : transcripts)
        {
            transIdMap.put(transcript, breakendId);

            final BreakendGeneData gene = transcript.gene();

            breakends.add(ImmutableLinxBreakend.builder()
                    .id(breakendId++)
                    .svId(transcript.gene().id())
                    .isStart(transcript.gene().isStart())
                    .type(gene.type())
                    .gene(transcript.geneName())
                    .transcriptId(transcript.transName())
                    .canonical(transcript.isCanonical())
                    .geneOrientation(transcript.isUpstream() ? BREAKEND_ORIENTATION_UPSTREAM : BREAKEND_ORIENTATION_DOWNSTREAM)
                    .disruptive(transcript.isDisruptive())
                    .reportedDisruption(transcript.reportableDisruption())
                    .junctionCopyNumber(gene.jcn())
                    .undisruptedCopyNumber(transcript.undisruptedCopyNumber())
                    .regionType(transcript.regionType())
                    .codingType(transcript.codingType())
                    .biotype(transcript.bioType())
                    .exonicBasePhase(transcript.Phase)
                    .nextSpliceExonRank(transcript.nextSpliceExonRank())
                    .nextSpliceExonPhase(transcript.Phase)
                    .nextSpliceDistance(transcript.isUpstream()
                            ? transcript.prevSpliceAcceptorDistance()
                            : transcript.nextSpliceAcceptorDistance())
                    .totalExonCount(transcript.TransData.exons().size())
                    .chromosome(gene.chromosome())
                    .orientation(gene.orientation())
                    .strand(gene.strand())
                    .chrBand(gene.GeneData.KaryotypeBand)
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
                    .reportableReasons(geneFusion.reportableReasonsStr())
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
                    .geneContextStart(context(geneFusion.upstreamTrans()
                            .regionType(), geneFusion.knownType(), geneFusion.getFusedExon(true)))
                    .geneEnd(geneFusion.geneName(FS_DOWN))
                    .geneTranscriptEnd(geneFusion.downstreamTrans().transName())
                    .geneContextEnd(context(geneFusion.downstreamTrans()
                            .regionType(), geneFusion.knownType(), geneFusion.getFusedExon(false)))
                    .junctionCopyNumber(fusionJcn(geneFusion.upstreamTrans().gene().jcn(), geneFusion.downstreamTrans().gene().jcn()))
                    .build());
        }
    }

    public void writeSampleData(final String sampleId, final List<LinxFusion> fusions, final List<LinxBreakend> breakends)
    {
        if(mOutputDir == null)
        {
            return;
        }

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

    private static final String COHORT_WRITER_FUSION = "Fusion";

    @Override
    public String fileType()
    {
        return COHORT_WRITER_FUSION;
    }

    @Override
    public BufferedWriter createWriter(final String outputDir)
    {
        // initialise the integrated, verbose fusion and breakend output file
        try
        {
            BufferedWriter writer = createBufferedWriter(outputDir + "LNX_FUSIONS.csv", false);

            writer.write("SampleId,Reportable,ReportableReason,KnownType,Phased,KnownExons,ClusterId,ClusterCount,ResolvedType");

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
                writer.write(fieldsStr);
            }

            writer.write(",ProteinsKept,ProteinsLost,PriorityScore,FusionId,ChainTerminatedUp,ChainTerminatedDown,ChainInfo");
            writer.newLine();
            return writer;
        }
        catch(final IOException e)
        {
            LNX_LOGGER.error("error writing fusions: {}", e.toString());
            return null;
        }
    }

    public void writeVerboseFusionData(final GeneFusion fusion, final String sampleId)
    {
        final FusionAnnotations annotations = fusion.getAnnotations();

        if(annotations == null)
        {
            LNX_LOGGER.error("fusion({}) annotations not set", fusion.name());
            return;
        }

        StringBuilder sb = new StringBuilder();

        sb.append(String.format("%s,%s,%s,%s",
                sampleId, fusion.reportable(), fusion.reportableReasonsStr(), fusion.knownTypeStr()));

        sb.append(String.format(",%s,%s,%d,%d,%s",
                fusion.phaseType(), fusion.knownExons(),
                annotations.clusterId(), annotations.clusterCount(), annotations.resolvedType()));

        // write upstream SV, transcript and exon info
        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            boolean isUpstream = (fs == FS_UP);
            final BreakendTransData trans = fusion.transcripts()[fs];
            final BreakendGeneData gene = trans.gene();

            sb.append(String.format(",%d,%s,%d,%d,%s,%.6f",
                    gene.id(), gene.chromosome(), gene.position(), gene.orientation(),
                    gene.type(), gene.jcn()));

            sb.append(String.format(",%s,%s,%s,%d,%s,%s",
                    gene.geneId(), fusion.geneName(fs), trans.transName(),
                    gene.strand(), trans.regionType(), trans.codingType()));

            sb.append(String.format(",%d,%d,%d,%d,%d,%s",
                    isUpstream ? trans.ExonUpstream : trans.ExonDownstream,
                    fusion.getFusedExon(isUpstream), fusion.getExonsSkipped()[fs],
                    trans.Phase, trans.exonCount(), trans.isDisruptive()));

            sb.append(String.format(",%d,%d,%d,%d,%d,%d,%d,%s,%s",
                    trans.CodingBases, trans.TotalCodingBases,
                    trans.codingStart(), trans.codingEnd(), trans.transStart(), trans.transEnd(),
                    trans.prevSpliceAcceptorDistance(), trans.isCanonical(), trans.bioType()));
        }

        sb.append(String.format(",%s,%s,%f,%d",
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

        sb.append(String.format("%s", chainInfo));

        mCohortDataWriter.write(this, Lists.newArrayList(sb.toString()));
    }
}
