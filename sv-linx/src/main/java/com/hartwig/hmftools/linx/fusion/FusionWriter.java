package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile.context;
import static com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile.fusionPloidy;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.FusionAnnotations;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableDisruptionFile;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ReportableGeneFusionFile;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FusionWriter
{
    private final String mOutputDir;
    private BufferedWriter mFusionWriter;

    private static final Logger LOGGER = LogManager.getLogger(FusionWriter.class);

    public FusionWriter(final String outputDir)
    {
        mOutputDir = outputDir;
        mFusionWriter = null;

    }

    public void writeSampleData(final String sampleId, final List<GeneFusion> fusions, final List<Transcript> disruptions)
    {
        // write sample files for patient reporter
        List<ReportableGeneFusion> reportedFusions = Lists.newArrayList();
        List<ReportableDisruption> reportedDisruptions = Lists.newArrayList();
        for(final GeneFusion fusion : fusions)
        {
            if(fusion.reportable())
            {
                reportedFusions.add(ImmutableReportableGeneFusion.builder()
                        .geneStart(fusion.upstreamTrans().geneName())
                        .geneTranscriptStart(fusion.upstreamTrans().StableId)
                        .geneContextStart(context(fusion.upstreamTrans(), fusion.getFusedExon(true)))
                        .geneEnd(fusion.downstreamTrans().geneName())
                        .geneTranscriptEnd(fusion.downstreamTrans().StableId)
                        .geneContextEnd(context(fusion.downstreamTrans(), fusion.getFusedExon(false)))
                        .ploidy(fusionPloidy(fusion.upstreamTrans().parent().ploidy(), fusion.downstreamTrans().parent().ploidy()))
                        .build());
            }
        }

        for(final Transcript transcript : disruptions)
        {
            final GeneAnnotation gene = transcript.parent();

            reportedDisruptions.add(ImmutableReportableDisruption.builder()
                    .svId(gene.id())
                    .chromosome(gene.chromosome())
                    .orientation(gene.orientation())
                    .strand(gene.Strand)
                    .chrBand(gene.karyotypeBand())
                    .gene(transcript.geneName())
                    .type(gene.type().toString())
                    .ploidy(gene.ploidy())
                    .exonUp(transcript.ExonUpstream)
                    .exonDown(transcript.ExonDownstream)
                    .build());
        }

        try
        {
            final String fusionsFile = ReportableGeneFusionFile.generateFilename(mOutputDir, sampleId);
            ReportableGeneFusionFile.write(fusionsFile, reportedFusions);

            final String disruptionsFile = ReportableDisruptionFile.generateFilename(mOutputDir, sampleId);
            ReportableDisruptionFile.write(disruptionsFile, reportedDisruptions);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write fusions file: {}", e.toString());
        }
    }
    public void initialiseOutputFile(final String fileName)
    {
        try
        {
            if(mFusionWriter == null)
            {
                String outputFilename = mOutputDir;

                if (!outputFilename.endsWith(File.separator))
                    outputFilename += File.separator;

                outputFilename += fileName;

                mFusionWriter = createBufferedWriter(outputFilename, false);

                mFusionWriter.write("SampleId,Reportable,KnownType,PhaseMatched,ClusterId,ClusterCount,ResolvedType");

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
                    fieldsStr += ",ExactBase" + upDown;
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

                mFusionWriter.write(",ProteinsKept,ProteinsLost,OverlapUp,OverlapDown,ChainInfo");
                mFusionWriter.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing fusions: {}", e.toString());
        }
    }

    public void writeFusionData(final GeneFusion fusion, final String sampleId)
    {
        if(mFusionWriter == null)
            return;

        try
        {
            BufferedWriter writer = mFusionWriter;

            final FusionAnnotations annotations = fusion.getAnnotations();

            if(annotations == null)
            {
                LOGGER.error("annotations not set");
                return;
            }

            writer.write(String.format("%s,%s,%s",
                    sampleId, fusion.reportable(), fusion.getKnownType()));

            writer.write(String.format(",%s,%d,%d,%s",
                    fusion.phaseMatched(), annotations.clusterId(), annotations.clusterCount(), annotations.resolvedType()));

            // write upstream SV, transcript and exon info
            for(int se = SE_START; se <= SE_END; ++se)
            {
                boolean isUpstream = (se == SE_START);
                final Transcript trans = isUpstream ? fusion.upstreamTrans() : fusion.downstreamTrans();
                final GeneAnnotation gene = trans.parent();

                writer.write(String.format(",%d,%s,%d,%d,%s,%.6f",
                        gene.id(), gene.chromosome(), gene.position(), gene.orientation(),
                        gene.type(), gene.ploidy()));

                writer.write(String.format(",%s,%s,%s,%d,%s,%s",
                        gene.StableId, gene.GeneName, trans.StableId,
                        gene.Strand, trans.regionType(), trans.codingType()));

                writer.write(String.format(",%d,%d,%d,%d,%d,%s",
                        isUpstream ? trans.ExonUpstream : trans.ExonDownstream,
                        fusion.getFusedExon(isUpstream), fusion.getExonsSkipped(isUpstream),
                        isUpstream ? trans.ExonUpstreamPhase : trans.ExonDownstreamPhase,
                        trans.ExonMax, trans.isDisruptive()));

                writer.write(String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%s,%s",
                        trans.exactCodingBase(), trans.calcCodingBases(true), trans.totalCodingBases(),
                        trans.codingStart(), trans.codingEnd(), trans.TranscriptStart, trans.TranscriptEnd,
                        trans.exonDistanceUp(), trans.isCanonical(), trans.bioType()));
            }

            writer.write(String.format(",%s,%s",
                        fusion.downstreamTrans().getProteinFeaturesKept(), fusion.downstreamTrans().getProteinFeaturesLost()));

            String chainInfo = "";
            String defaultValues = ",0:false;0;0;0;false";
            if(annotations.disruptionUp() != null)
            {
                chainInfo += String.format(",%d;%s;%d;%d;%d;%s",
                        annotations.disruptionUp().facingBreakends(), annotations.disruptionUp().allLinksAssembled(),
                        annotations.disruptionUp().totalBreakends(), annotations.disruptionUp().minDistance(),
                        annotations.disruptionUp().disruptedExons(), annotations.disruptionUp().transcriptTerminated());
            }
            else
            {
                chainInfo += defaultValues;

            }

            if(annotations.disruptionDown() != null)
            {
                chainInfo += String.format(",%d;%s;%d;%d;%d;%s",
                        annotations.disruptionDown().facingBreakends(), annotations.disruptionDown().allLinksAssembled(),
                        annotations.disruptionDown().totalBreakends(), annotations.disruptionDown().minDistance(),
                        annotations.disruptionDown().disruptedExons(), annotations.disruptionDown().transcriptTerminated());
            }
            else
            {
                chainInfo += defaultValues;
            }

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
            LOGGER.error("error writing fusions: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mFusionWriter);
    }


}
