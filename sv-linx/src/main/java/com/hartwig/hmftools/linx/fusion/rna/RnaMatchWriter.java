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

                mWriter.write("SampleId,FusionName,GeneNameUp,GeneNameDown,ViableFusion");
                mWriter.write(",PhaseMatched,DnaMatchType,DnaMatchInfo,KnownType,RnaPhaseMatched");

                for(int se = SE_START; se <= SE_END; ++se)
                {
                    String upDown = se == SE_START ? "Up" : "Down";

                    String fieldsStr = ",SvId" + upDown;
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
            mWriter.write(String.format("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s",
                    sampleId, rnaFusion.name(), rnaFusion.GeneIds[FS_UPSTREAM], rnaFusion.GeneIds[FS_DOWNSTREAM],
                    rnaFusion.GeneNames[FS_UPSTREAM], rnaFusion.GeneNames[FS_DOWNSTREAM],
                    rnaFusion.isViableFusion(), rnaFusion.isPhaseMatchedFusion(), rnaFusion.getDnaFusionMatchType(),
                    rnaFusion.getDnaFusionMatchInfo(), rnaFusion.getKnownType(), rnaFusion.hasRnaPhasedFusion()));

            for(int se = SE_START; se <= SE_END; ++se)
            {
                boolean isUpstream = (se == SE_START);
                final Transcript trans = rnaFusion.getTrans(isUpstream);

                if(trans != null)
                {
                    final GeneAnnotation gene = trans.gene();

                    mWriter.write(String.format(",%d,%s,%d,%d,%d,%d,%s,%s",
                            gene.id(), gene.chromosome(), gene.position(),
                            isUpstream ? rnaFusion.Positions[FS_UPSTREAM] : rnaFusion.Positions[FS_DOWNSTREAM],
                            gene.orientation(), gene.Strand, gene.type(),
                            rnaFusion.getClusterInfo(isUpstream)));

                    mWriter.write(String.format(",%s,%s,%s,%d,%s,%s,%d,%s,%d",
                            rnaFusion.isTransViable(isUpstream), rnaFusion.isTransCorrectLocation(isUpstream),
                            trans.StableId, rnaFusion.getExonsSkipped(isUpstream),
                            trans.regionType(), trans.codingType(),
                            isUpstream ? trans.ExonUpstream : trans.ExonDownstream, trans.isDisruptive(), trans.prevSpliceAcceptorDistance()));
                }
                else
                {
                    mWriter.write(String.format(",%s,%s,%d,%d,%d,%d,%s,%s",
                            "", isUpstream ? rnaFusion.Chromosomes[FS_UPSTREAM] : rnaFusion.Chromosomes[FS_DOWNSTREAM], 0,
                            isUpstream ? rnaFusion.Positions[FS_UPSTREAM] : rnaFusion.Positions[FS_DOWNSTREAM], 0, 0, "", ""));

                    mWriter.write(String.format(",%s,%s,,,,,,,",
                            rnaFusion.isTransViable(isUpstream), rnaFusion.isTransCorrectLocation(isUpstream)));
                }

                boolean hasExonData = rnaFusion.exonRank(isUpstream) > 0;

                mWriter.write(String.format(",%s,%s,%d,%s",
                        hasExonData ? rnaFusion.exonRank(isUpstream) :"", hasExonData ? rnaFusion.exonPhase(isUpstream) : "",
                        rnaFusion.getExactMatchTransIds(isUpstream).size(), rnaFusion.getRnaPhasedFusionTransId(isUpstream)));
            }

            mWriter.write(String.format(",%s,%d,%d,%s,%s,%s",
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
