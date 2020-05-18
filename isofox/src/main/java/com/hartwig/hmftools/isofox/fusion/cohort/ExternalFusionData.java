package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class ExternalFusionData
{

    /*
    mutate(breakpoint1=paste(ChrUp,PosUp,sep=":"),breakpoint2=paste(ChrDown,PosDown,sep=":")),by=c('breakpoint1','breakpoint2'),all.x=T),
           passFusions %>% mutate(breakpoint1=paste(ChrUp,PosUp,sep=":"),breakpoint2=paste(ChrDown,PosDown,sep=":"),revOrientation=TRUE),by.x=c('breakpoint1','breakpoint2'),by.y=c('breakpoint2','breakpoint1'),all.x=T))# %>%
       filter(is.na(FusionId.x),is.na(FusionId.y),!grepl("read-through",type)))
     */

    public static final String FUSION_SOURCE_ARRIBA = "ARRIBA";

    public final String SourceTool;
    public final String[] Chromosomes;
    public final int[] JunctionPositions;
    public final byte[] JunctionOrientations;
    public final String SvType;
    public final String[] GeneNames;

    public final int SplitFragments;
    public final int DiscordantFragments;
    public final int[] Coverage;

    public final Map<String,String> OtherData;

    public ExternalFusionData(final String sourceTool, final String[] chromosomes, final int[] junctionPositions,
            final byte[] junctionOrientations, final String svType, final String[] geneNames, final int splitFragments,
            final int discordantFragments, final int[] coverage, final Map<String, String> otherData)
    {
        SourceTool = sourceTool;
        Chromosomes = chromosomes;
        JunctionPositions = junctionPositions;
        JunctionOrientations = junctionOrientations;
        SvType = svType;
        GeneNames = geneNames;
        SplitFragments = splitFragments;
        DiscordantFragments = discordantFragments;
        Coverage = coverage;
        OtherData = otherData;
    }

    public static ExternalFusionData loadArribaFusion(final String data)
    {
        final String[] items = data.split("\t", -1);

        // 0       1       2                       3                       4               5               6       7       8       9               10
        // #gene1  gene2   strand1(gene/fusion)    strand2(gene/fusion)    breakpoint1     breakpoint2     site1   site2   type    direction1      direction2
        // 11              12              13                      14              15
        // split_reads1    split_reads2    discordant_mates        coverage1       coverage2
        // confidence      closest_genomic_breakpoint1     closest_genomic_breakpoint2     filters fusion_transcript
        // reading_frame   peptide_sequence        read_identifiers

        final String[] geneNames = new String[] { items[0], items[1] };
        final String[] chromosomes = { items[4].split(":")[0], items[5].split(":")[0] };

        if(!HumanChromosome.contains(chromosomes[SE_START]) || !HumanChromosome.contains(chromosomes[SE_END]))
            return null;

        final int[] positions = { Integer.parseInt(items[4].split(":")[1]), Integer.parseInt(items[5].split(":")[1]) };
        final byte[] orientations = {0, 0};
        String svType = items[8];

        int splitFragments = Integer.parseInt(items[11]) + Integer.parseInt(items[12]);
        int discordantFragments = Integer.parseInt(items[13]);

        int coverage1 = items[14].equals(".") ? 0 : Integer.parseInt(items[14]);
        int coverage2 = items[15].equals(".") ? 0 : Integer.parseInt(items[15]);
        final int[] coverage = new int[] { coverage1, coverage2 };
        final Map<String, String> otherData = Maps.newHashMap();

        return new ExternalFusionData(FUSION_SOURCE_ARRIBA, chromosomes, positions, orientations, svType, geneNames,
                splitFragments, discordantFragments, coverage, otherData);
    }
}
