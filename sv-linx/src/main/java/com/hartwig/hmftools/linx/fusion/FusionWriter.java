package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

import com.hartwig.hmftools.common.variant.structural.annotation.FusionAnnotations;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
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

                mFusionWriter.write("SampleId,Reportable,KnownType");

                mFusionWriter.write(",PhaseMatched,ClusterId,ClusterCount,ResolvedType");

                mFusionWriter.write(",SvIdUp,ChrUp,PosUp,OrientUp,TypeUp,PloidyUp,GeneIdUp,GeneNameUp,ChrBandUp,TranscriptUp,StrandUp,RegionTypeUp,CodingTypeUp");
                mFusionWriter.write(",ExonUp,PhaseUp,ExonMaxUp,DisruptiveUp,ExactBaseUp,CodingBasesUp,TotalCodingUp");
                mFusionWriter.write(",CodingStartUp,CodingEndUp,TransStartUp,TransEndUp,DistancePrevUp,CanonicalUp,BiotypeUp,ExonsSkippedUp");

                mFusionWriter.write(",SvIdDown,ChrDown,PosDown,OrientDown,TypeDown,PloidyDown,GeneIdDown,GeneNameDown,ChrBandDown,TranscriptDown,StrandDown,RegionTypeDown,CodingTypeDown");
                mFusionWriter.write(",ExonDown,PhaseDown,ExonMaxDown,DisruptiveDown,ExactBaseDown,CodingBasesDown,TotalCodingDown");
                mFusionWriter.write(",CodingStartDown,CodingEndDown,TransStartDown,TransEndDown,DistancePrevDown,CanonicalDown,BiotypeDown,ExonsSkippedDown");

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

            final Transcript upTrans = fusion.upstreamTrans();
            final Transcript downTrans = fusion.downstreamTrans();

            final GeneAnnotation upGene = upTrans.parent();
            final GeneAnnotation downGene = downTrans.parent();

            final FusionAnnotations annotations = fusion.getAnnotations();
            String annotationsStr = ",,";

            if(annotations != null)
            {
                // PhaseMatched,ClusterId,ClusterCount,ResolvedType,OverlapUp,OverlapDown,ChainInfo
                annotationsStr = String.format("%s,%d,%d,%s",
                        fusion.phaseMatched(), annotations.clusterId(), annotations.clusterCount(), annotations.resolvedType());

                String defaultValues = ",0:false;0;0;0;false";
                if(annotations.disruptionUp() != null)
                {
                    annotationsStr += String.format(",%d;%s;%d;%d;%d;%s",
                            annotations.disruptionUp().facingBreakends(), annotations.disruptionUp().allLinksAssembled(),
                            annotations.disruptionUp().totalBreakends(), annotations.disruptionUp().minDistance(),
                            annotations.disruptionUp().disruptedExons(), annotations.disruptionUp().transcriptTerminated());
                }
                else
                {
                    annotationsStr += defaultValues;

                }

                if(annotations.disruptionDown() != null)
                {
                    annotationsStr += String.format(",%d;%s;%d;%d;%d;%s",
                            annotations.disruptionDown().facingBreakends(), annotations.disruptionDown().allLinksAssembled(),
                            annotations.disruptionDown().totalBreakends(), annotations.disruptionDown().minDistance(),
                            annotations.disruptionDown().disruptedExons(), annotations.disruptionDown().transcriptTerminated());
                }
                else
                {
                    annotationsStr += defaultValues;
                }

                if(annotations.chainInfo() != null)
                {
                    annotationsStr += String.format(",%d;%d;%d;%s;%s",
                            annotations.chainInfo().chainId(), annotations.chainInfo().chainLinks(), annotations.chainInfo().chainLength(),
                            annotations.chainInfo().validTraversal(), annotations.chainInfo().traversalAssembled());
                }
                else
                {
                    annotationsStr += ",-1;0;0;true;false";
                }
            }

            writer.write(String.format("%s,%s,%s,%s",
                    sampleId, fusion.reportable(), fusion.getKnownFusionType(), annotationsStr));

            // write upstream SV, transcript and exon info
            writer.write(
                    String.format(",%d,%s,%d,%d,%s,%.6f",
                            upGene.id(), upGene.chromosome(), upGene.position(), upGene.orientation(),
                            upGene.type(), upGene.ploidy()));

            writer.write(
                    String.format(",%s,%s,%s,%s,%d,%s,%s",
                            upGene.StableId, upGene.GeneName, upGene.karyotypeBand(), upTrans.StableId,
                            upGene.Strand, upTrans.regionType(), upTrans.codingType()));

            writer.write(
                    String.format(",%d,%d,%d,%s",
                            upTrans.ExonUpstream, upTrans.ExonUpstreamPhase, upTrans.ExonMax, upTrans.isDisruptive()));
            writer.write(
                    String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%s,%s",
                            upTrans.exactCodingBase(), upTrans.calcCodingBases(true), upTrans.totalCodingBases(),
                            upTrans.codingStart(), upTrans.codingEnd(), upTrans.TranscriptStart, upTrans.TranscriptEnd,
                            upTrans.exonDistanceUp(), upTrans.isCanonical(), upTrans.bioType()));

            writer.write(
                    String.format(",%d,%s,%d,%d,%s,%.6f",
                            downGene.id(), downGene.chromosome(), downGene.position(), downGene.orientation(),
                            downGene.type(), downGene.ploidy()));

            writer.write(
                    String.format(",%s,%s,%s,%s,%d,%s,%s",
                            downGene.StableId, downGene.GeneName, downGene.karyotypeBand(), downTrans.StableId,
                            downGene.Strand, downTrans.regionType(), downTrans.codingType()));

            writer.write(
                    String.format(",%d,%d,%d,%s",
                            downTrans.ExonDownstream, downTrans.ExonDownstreamPhase, downTrans.ExonMax, downTrans.isDisruptive()));

            writer.write(
                    String.format(",%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s,%s",
                            downTrans.exactCodingBase(), downTrans.calcCodingBases(false), downTrans.totalCodingBases(),
                            downTrans.codingStart(), downTrans.codingEnd(), downTrans.TranscriptStart, downTrans.TranscriptEnd,
                            downTrans.exonDistanceUp(), downTrans.isCanonical(), downTrans.bioType(),
                            downTrans.getProteinFeaturesKept(), downTrans.getProteinFeaturesLost()));

            // move to phasing section once move to SVA fusions
            writer.write(String.format(",%d,%d", fusion.getExonsSkipped(true), fusion.getExonsSkipped(false)));

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
