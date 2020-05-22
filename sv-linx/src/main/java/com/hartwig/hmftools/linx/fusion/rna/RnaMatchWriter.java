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
import java.util.StringJoiner;

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

                mWriter.write("SampleId,Source,FusionId,FusionName,ViableFusion");
                mWriter.write(",PhaseMatched,DnaMatchType,DnaMatchInfo,KnownType,RnaPhaseMatched");
                mWriter.write(",JunctionFrags,DiscordantFrags");

                for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
                {
                    String upDown = fs == FS_UPSTREAM ? "Up" : "Down";
                    StringJoiner fieldsStr = new StringJoiner(",");

                    fieldsStr.add("GeneId" + upDown);
                    fieldsStr.add("GeneName" + upDown);
                    fieldsStr.add("Chr" + upDown);
                    fieldsStr.add("RnaPos" + upDown);
                    fieldsStr.add("RnaJunc" + upDown);
                    fieldsStr.add("Strand" + upDown);
                    fieldsStr.add("SvId" + upDown);
                    fieldsStr.add("SvPos" + upDown);
                    fieldsStr.add("Orient" + upDown);
                    fieldsStr.add("Type" + upDown);
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

                mWriter.write(",ChainInfo");

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
                    rnaFusion.isViableFusion(), rnaFusion.isPhaseMatchedFusion(), rnaFusion.getDnaFusionMatchType(),
                    rnaFusion.getDnaFusionMatchInfo(), rnaFusion.getKnownType(), rnaFusion.hasRnaPhasedFusion()));

            mWriter.write(String.format(",%d,%d",
                    rnaFusion.JunctionFragments, rnaFusion.DiscordantFragments));

            for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
            {
                boolean isUpstream = (fs == SE_START);

                mWriter.write(String.format(",%s,%s,%s,%d,%s,%d",
                        rnaFusion.GeneIds[fs], rnaFusion.GeneNames[fs], rnaFusion.Chromosomes[fs],
                        rnaFusion.Positions[fs], rnaFusion.JunctionTypes[fs], rnaFusion.Strands[fs]));

                final Transcript trans = rnaFusion.getMatchedfTranscripts()[fs];

                if(trans != null)
                {
                    final GeneAnnotation gene = trans.gene();

                    // SV breakend info
                    mWriter.write(String.format(",%d,%d,%d,%s,%s",
                            gene.id(), gene.position(), gene.orientation(), gene.type(), rnaFusion.getClusterInfo(isUpstream)));

                    // transcript info and viability
                    mWriter.write(String.format(",%s,%s,%s,%d,%s,%s,%d,%s,%d",
                            rnaFusion.isTransViable()[fs], rnaFusion.isTransCorrectLocation()[fs],
                            trans.StableId, rnaFusion.getExonsSkipped()[fs],
                            trans.regionType(), trans.codingType(),
                            isUpstream ? trans.ExonUpstream : trans.ExonDownstream,
                            trans.isDisruptive(), trans.prevSpliceAcceptorDistance()));
                }
                else
                {
                    mWriter.write(String.format(",%d,%d,%d,%s,%s",
                            0, 0, 0, "", ""));

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

            mWriter.write(String.format(",%s", !rnaFusion.getChainInfo().isEmpty() ? rnaFusion.getChainInfo() : "0;0"));
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
