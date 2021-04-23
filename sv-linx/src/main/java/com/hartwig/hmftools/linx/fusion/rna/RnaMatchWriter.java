package com.hartwig.hmftools.linx.fusion.rna;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.fusion.rna.RnaFusionData.NO_CLUSTER_INFO;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.BreakendTransData;

public class RnaMatchWriter
{
    private BufferedWriter mWriter;

    public RnaMatchWriter(final String outputDir, final String outputId)
    {
        if(outputDir != null)
            initialiseWriter(outputDir, outputId);
    }

    private void initialiseWriter(final String outputDir, final String outputId)
    {
        try
        {
            if(mWriter == null)
            {
                String outputFilename = outputDir;

                outputFilename += String.format("LNX_RNA_FUSION_MATCH_%s.csv", outputId);

                mWriter = createBufferedWriter(outputFilename, false);

                mWriter.write("SampleId,Source,FusionId,FusionName,ViableFusion");
                mWriter.write(",PhaseMatched,DnaFusionMatchType,Reportable,KnownType,RnaPhaseMatched");
                mWriter.write(",RnaSvType,JunctionFrags,DiscordantFrags");

                for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
                {
                    String upDown = fs == FS_UP ? "Up" : "Down";
                    StringJoiner fieldsStr = new StringJoiner(",");

                    fieldsStr.add("GeneId" + upDown);
                    fieldsStr.add("GeneName" + upDown);
                    fieldsStr.add("Chr" + upDown);
                    fieldsStr.add("RnaPos" + upDown);
                    fieldsStr.add("RnaOrient" + upDown);
                    fieldsStr.add("RnaJunc" + upDown);
                    fieldsStr.add("Strand" + upDown);
                    fieldsStr.add("SvId" + upDown);
                    fieldsStr.add("SvPos" + upDown);
                    fieldsStr.add("SvType" + upDown);
                    fieldsStr.add("ClusterInfo" + upDown);
                    fieldsStr.add("TransViable" + upDown);
                    fieldsStr.add("TransValidLoc" + upDown);
                    fieldsStr.add("TransId" + upDown);
                    fieldsStr.add("ExonsSkipped" + upDown);
                    fieldsStr.add("RegionType" + upDown);
                    fieldsStr.add("CodingType" + upDown);
                    fieldsStr.add("Exon" + upDown);
                    fieldsStr.add("Disruptive" + upDown);
                    fieldsStr.add("DistancePrev" + upDown);
                    fieldsStr.add("RnaTransId" + upDown);
                    fieldsStr.add("RnaExonMatch" + upDown);
                    fieldsStr.add("RnaExonRank" + upDown);
                    fieldsStr.add("RnaExonPhase" + upDown);
                    mWriter.write(String.format(",%s", fieldsStr.toString()));
                }

                mWriter.write(",ChainInfo,RnaCohortCount,RnaOtherData");

                mWriter.newLine();
            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing RNA match data: {}", e.toString());
        }
    }

    public void writeRnaMatchData(final RnaFusionData rnaFusion)
    {
        if(mWriter == null)
            return;

        try
        {
            mWriter.write(String.format("%s,%s,%s,%s",
                    rnaFusion.SampleId, rnaFusion.Source, rnaFusion.FusionId, rnaFusion.name()));

            mWriter.write(String.format(",%s,%s,%s,%s,%s,%s",
                    rnaFusion.isViableFusion(), rnaFusion.isPhaseMatchedFusion(), rnaFusion.getDnaFusionMatchInfo(),
                    rnaFusion.hasDnaReportableFusion(), rnaFusion.getKnownType(), rnaFusion.hasRnaPhasedFusion()));

            mWriter.write(String.format(",%s,%d,%d",
                    rnaFusion.rnaSvType(), rnaFusion.JunctionFragments, rnaFusion.DiscordantFragments));

            for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
            {
                boolean isUpstream = (fs == SE_START);

                mWriter.write(String.format(",%s,%s,%s,%d,%d,%s,%d",
                        rnaFusion.GeneIds[fs], rnaFusion.GeneNames[fs], rnaFusion.Chromosomes[fs],
                        rnaFusion.Positions[fs], rnaFusion.Orientations[fs], rnaFusion.JunctionTypes[fs], rnaFusion.Strands[fs]));

                final BreakendTransData trans = rnaFusion.getMatchedTranscripts()[fs];

                if(trans != null)
                {
                    final BreakendGeneData gene = trans.gene();

                    // SV breakend info
                    mWriter.write(String.format(",%d,%d,%s,%s",
                            gene.id(), gene.position(), gene.type(), rnaFusion.getClusterInfo(isUpstream)));

                    // transcript info and viability
                    mWriter.write(String.format(",%s,%s,%s,%d,%s,%s,%d,%s,%d",
                            rnaFusion.isTransViable()[fs], rnaFusion.isTransCorrectLocation()[fs],
                            trans.transName(), rnaFusion.getExonsSkipped()[fs],
                            trans.regionType(), trans.codingType(),
                            isUpstream ? trans.ExonUpstream : trans.ExonDownstream,
                            trans.isDisruptive(), trans.prevSpliceAcceptorDistance()));
                }
                else
                {
                    mWriter.write(String.format(",%d,%d,%s,%s",
                            0, 0, "", NO_CLUSTER_INFO));

                    mWriter.write(String.format(",%s,%s,,,,,,,",
                            rnaFusion.isTransViable()[fs], rnaFusion.isTransCorrectLocation()[fs]));
                }

                RnaExonMatchData rnaExonData = rnaFusion.getBestExonMatch(fs);

                if(rnaExonData != null)
                {
                    mWriter.write(String.format(",%s,%s,%d,%d",
                            rnaExonData.TransName, rnaExonData.BoundaryMatch, rnaExonData.ExonRank, rnaExonData.ExonPhase));
                }
                else
                {
                    mWriter.write(",NONE,false,0,0");
                }
            }

            mWriter.write(String.format(",%s,%d,%s", rnaFusion.getChainInfo(), rnaFusion.CohortCount, rnaFusion.OtherData));
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
