package com.hartwig.hmftools.linx.fusion.rna;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.Transcript;

public class RnaMatchWriter
{
    private BufferedWriter mWriter;

    public RnaMatchWriter(final String outputDir)
    {
        if(outputDir != null)
            initialiseWriter(outputDir);
    }

    private void initialiseWriter(final String outputDir)
    {
        try
        {
            if(mWriter == null)
            {
                String outputFilename = outputDir;

                outputFilename += "LNX_RNA_FUSION_MATCH.csv";

                mWriter = createBufferedWriter(outputFilename, false);

                mWriter.write("SampleId,FusionId,FusionName,ViableFusion");
                mWriter.write(",PhaseMatched,DnaMatchType,DnaMatchInfo,KnownType,RnaPhaseMatched");

                for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
                {
                    String upDown = fs == FS_UPSTREAM ? "Up" : "Down";
                    String fieldsStr = "";

                    fieldsStr += ",GeneId" + upDown;
                    fieldsStr += ",GeneName" + upDown;
                    fieldsStr += ",SvId" + upDown;
                    fieldsStr += ",Chr" + upDown;
                    fieldsStr += ",Pos" + upDown;
                    fieldsStr += ",RnaPos" + upDown;
                    fieldsStr += ",Orient" + upDown;
                    fieldsStr += ",Strand" + upDown;
                    fieldsStr += ",Type" + upDown;
                    fieldsStr += ",ClusterInfo" + upDown;
                    fieldsStr += ",TransViable" + upDown;
                    fieldsStr += ",TransValidLoc" + upDown;
                    fieldsStr += ",TransId" + upDown;
                    fieldsStr += ",ExonsSkipped" + upDown;
                    fieldsStr += ",RegionType" + upDown;
                    fieldsStr += ",CodingType" + upDown;
                    fieldsStr += ",Exon" + upDown;
                    fieldsStr += ",Disruptive" + upDown;
                    fieldsStr += ",DistancePrev" + upDown;
                    fieldsStr += ",RnaExonRank" + upDown;
                    fieldsStr += ",RnaExonPhase" + upDown;
                    fieldsStr += ",RnaExonMatch" + upDown;
                    fieldsStr += ",RnaTransId" + upDown;
                    mWriter.write(fieldsStr);
                }

                mWriter.write(",ChainInfo,JunctionReadCount,SpanningFragCount,JuncTypeUp,JuncTypeDown");

                mWriter.newLine();
            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing RNA match data: {}", e.toString());
        }

    }

    public void writeRnaMatchData(final String sampleId, final RnaFusionData rnaFusion)
    {
        if(mWriter == null)
            return;

        try
        {
            mWriter.write(String.format("%s,%s,%s,%s,%s,%s,%s,%s,%s",
                    sampleId, rnaFusion.FusionId, rnaFusion.name(), rnaFusion.isViableFusion(), rnaFusion.isPhaseMatchedFusion(),
                    rnaFusion.getDnaFusionMatchType(), rnaFusion.getDnaFusionMatchInfo(), rnaFusion.getKnownType(), rnaFusion.hasRnaPhasedFusion()));

            for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
            {
                boolean isUpstream = (fs == SE_START);
                final Transcript trans = rnaFusion.getMatchedfTranscripts()[fs];

                if(trans != null)
                {
                    final GeneAnnotation gene = trans.gene();

                    mWriter.write(String.format(",%s,%s,%d,%s,%d,%d,%d,%d,%s,%s",
                            gene.StableId, gene.GeneName, gene.id(), gene.chromosome(), gene.position(),
                            rnaFusion.Positions[fs], gene.orientation(), gene.Strand, gene.type(),
                            rnaFusion.getClusterInfo(isUpstream)));

                    mWriter.write(String.format(",%s,%s,%s,%d,%s,%s,%d,%s,%d",
                            rnaFusion.isTransViable()[fs], rnaFusion.isTransCorrectLocation()[fs],
                            trans.StableId, rnaFusion.getExonsSkipped()[fs],
                            trans.regionType(), trans.codingType(),
                            isUpstream ? trans.ExonUpstream : trans.ExonDownstream, trans.isDisruptive(), trans.prevSpliceAcceptorDistance()));
                }
                else
                {
                    mWriter.write(String.format(",%s,%s,%d,%s,%d,%d,%d,%d,%s,%s",
                            "", "", 0, rnaFusion.Chromosomes[fs], 0, rnaFusion.Positions[fs], 0, 0, "", ""));

                    mWriter.write(String.format(",%s,%s,,,,,,,",
                            rnaFusion.isTransViable()[fs], rnaFusion.isTransCorrectLocation()[fs]));
                }

                boolean hasExonData = rnaFusion.exonRank()[fs] > 0;

                mWriter.write(String.format(",%s,%s,%d,%s",
                        hasExonData ? rnaFusion.exonRank()[fs] :"", hasExonData ? rnaFusion.exonPhase()[fs] : "",
                        rnaFusion.getExactMatchTransIds()[fs].size(), rnaFusion.getRnaPhasedFusionTransId()[fs]));
            }

            mWriter.write(String.format(",%s,%d,%d,%s,%s",
                    !rnaFusion.getChainInfo().isEmpty() ? rnaFusion.getChainInfo() : "0;0",
                    rnaFusion.JunctionFragments, rnaFusion.DiscordantFragments,
                    rnaFusion.JunctionTypes[FS_UPSTREAM], rnaFusion.JunctionTypes[FS_DOWNSTREAM]));

            mWriter.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing RNA match data: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mWriter);
    }

}
