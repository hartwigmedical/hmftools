package com.hartwig.hmftools.linx.fusion;

import static java.lang.String.valueOf;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.linx.LinxBreakend.BREAKEND_ORIENTATION_DOWNSTREAM;
import static com.hartwig.hmftools.common.linx.LinxBreakend.BREAKEND_ORIENTATION_UPSTREAM;
import static com.hartwig.hmftools.common.linx.LinxFusion.context;
import static com.hartwig.hmftools.common.linx.LinxFusion.fusionJcn;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.CohortDataWriter.cohortDataFilename;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

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

            BreakendGeneData geneData = transcript.breakendGeneData();

            breakends.add(ImmutableLinxBreakend.builder()
                    .id(breakendId++)
                    .svId(geneData.varId())
                    .vcfId(geneData.vcfId())
                    .coords(geneData.coordsStr())
                    .isStart(geneData.isStart())
                    .gene(transcript.geneName())
                    .transcriptId(transcript.transName())
                    .canonical(transcript.isCanonical())
                    .geneOrientation(transcript.isUpstream() ? BREAKEND_ORIENTATION_UPSTREAM : BREAKEND_ORIENTATION_DOWNSTREAM)
                    .disruptive(transcript.isDisruptive())
                    .reportedDisruption(transcript.reportableDisruption())
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
                    .exonUp(transcript.ExonUpstream)
                    .exonDown(transcript.ExonDownstream)
                    .build());
        }

        for(final GeneFusion geneFusion : geneFusions)
        {
            BreakendGeneData upGeneData = geneFusion.upstreamTrans().breakendGeneData();
            BreakendGeneData downGeneData = geneFusion.downstreamTrans().breakendGeneData();

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
                    .fivePrimeVcfId(upGeneData.vcfId())
                    .threePrimeVcfId(downGeneData.vcfId())
                    .fivePrimeCoords(upGeneData.coordsStr())
                    .threePrimeCoords(downGeneData.coordsStr())
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
                    .junctionCopyNumber(fusionJcn(geneFusion.upstreamTrans().breakendGeneData().jcn(), geneFusion.downstreamTrans().breakendGeneData().jcn()))
                    .build());
        }
    }

    public void writeSampleData(final String sampleId, final List<LinxFusion> fusions, final List<LinxBreakend> breakends)
    {
        if(mOutputDir == null)
            return;

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
            String cohortFilename = cohortDataFilename(outputDir, "FUSIONS");
            BufferedWriter writer = createBufferedWriter(cohortFilename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("SampleId").add("Reportable").add("ReportableReason").add("KnownType").add("Phased").add("KnownExons");
            sj.add("ClusterId").add("ClusterCount").add("ResolvedType");

            for(int se = SE_START; se <= SE_END; ++se)
            {
                String upDown = se == SE_START ? "Up" : "Down";

                sj.add("SvId" + upDown);
                sj.add("Chr" + upDown);
                sj.add("Pos" + upDown);
                sj.add("Orient" + upDown);
                sj.add("Type" + upDown);
                sj.add("Ploidy" + upDown);
                sj.add("GeneId" + upDown);
                sj.add("GeneName" + upDown);
                sj.add("Transcript" + upDown);
                sj.add("Strand" + upDown);
                sj.add("RegionType" + upDown);
                sj.add("CodingType" + upDown);
                sj.add("BreakendExon" + upDown);
                sj.add("FusedExon" + upDown);
                sj.add("ExonsSkipped" + upDown);
                sj.add("Phase" + upDown);
                sj.add("ExonMax" + upDown);
                sj.add("Disruptive" + upDown);
                sj.add("CodingBases" + upDown);
                sj.add("TotalCoding" + upDown);
                sj.add("CodingStart" + upDown);
                sj.add("CodingEnd" + upDown);
                sj.add("TransStart" + upDown);
                sj.add("TransEnd" + upDown);
                sj.add("DistancePrev" + upDown);
                sj.add("Canonical" + upDown);
                sj.add("Biotype" + upDown);
            }

            sj.add("ProteinsKept").add("ProteinsLost").add("PriorityScore").add("FusionId").add("ChainTerminatedUp");
            sj.add("ChainTerminatedDown").add("ChainInfo");

            writer.write(sj.toString());
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

        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(sampleId);
        sj.add(valueOf(fusion.reportable()));
        sj.add(fusion.reportableReasonsStr());
        sj.add(fusion.knownTypeStr());

        sj.add(valueOf(fusion.phaseType()));
        sj.add(valueOf(fusion.knownExons()));
        sj.add(valueOf(annotations.clusterId()));
        sj.add(valueOf(annotations.clusterCount()));
        sj.add(annotations.resolvedType());

        // write upstream SV, transcript and exon info
        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            boolean isUpstream = (fs == FS_UP);
            final BreakendTransData trans = fusion.transcripts()[fs];
            final BreakendGeneData gene = trans.breakendGeneData();

            sj.add(valueOf(gene.varId()));
            sj.add(gene.chromosome());
            sj.add(valueOf(gene.position()));
            sj.add(valueOf(gene.orientation()));
            sj.add(valueOf(gene.svType()));
            sj.add(valueOf(gene.jcn()));

            sj.add(gene.geneId());
            sj.add(fusion.geneName(fs));
            sj.add(trans.transName());
            sj.add(valueOf(gene.strand()));
            sj.add(valueOf(trans.regionType()));
            sj.add(valueOf(trans.codingType()));

            sj.add(valueOf(isUpstream ? trans.ExonUpstream : trans.ExonDownstream));
            sj.add(valueOf(fusion.getFusedExon(isUpstream)));
            sj.add(valueOf(fusion.getExonsSkipped()[fs]));
            sj.add(valueOf(trans.Phase));
            sj.add(valueOf(trans.exonCount()));
            sj.add(valueOf(trans.isDisruptive()));

            sj.add(valueOf(trans.CodingBases));
            sj.add(valueOf(trans.TotalCodingBases));
            sj.add(valueOf(trans.codingStart()));
            sj.add(valueOf(trans.codingEnd()));
            sj.add(valueOf(trans.transStart()));
            sj.add(valueOf(trans.transEnd()));
            sj.add(valueOf(trans.prevSpliceAcceptorDistance()));
            sj.add(valueOf(trans.isCanonical()));
            sj.add(trans.bioType());
        }

        sj.add(fusion.downstreamTrans().getProteinFeaturesKept());
        sj.add(fusion.downstreamTrans().getProteinFeaturesLost());
        sj.add(valueOf(fusion.priority()));
        sj.add(valueOf(fusion.id()));

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

        sj.add(String.format("%s", chainInfo));

        mCohortDataWriter.write(this, Lists.newArrayList(sj.toString()));
    }
}
