package com.hartwig.hmftools.isofox.fusion.cohort;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_ORIENT;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.isofox.fusion.FusionJunctionType;

public class ExternalFusionData
{
    public static final String FUSION_SOURCE_ARRIBA = "ARRIBA";

    public final String SourceTool;
    public final String[] Chromosomes;
    public final int[] JunctionPositions;
    public final byte[] JunctionOrientations;
    public final FusionJunctionType[] JunctionTypes;

    public final String SvType;
    public final String[] GeneNames;

    public final int SplitFragments;
    public final int DiscordantFragments;
    public final int[] Coverage;

    public final Map<String,String> OtherData;

    public boolean IsFiltered;

    public ExternalFusionData(final String sourceTool, final String[] chromosomes, final int[] junctionPositions,
            final byte[] junctionOrientations, final FusionJunctionType[] junctionTypes, final String svType, final String[] geneNames,
            final int splitFragments,final int discordantFragments, final int[] coverage, final Map<String, String> otherData)
    {
        SourceTool = sourceTool;
        Chromosomes = chromosomes;
        JunctionPositions = junctionPositions;
        JunctionOrientations = junctionOrientations;
        JunctionTypes = junctionTypes;
        SvType = svType;
        GeneNames = geneNames;
        SplitFragments = splitFragments;
        DiscordantFragments = discordantFragments;
        Coverage = coverage;
        IsFiltered = true;
        OtherData = otherData;
    }

    public static ExternalFusionData loadArribaFusion(final String data)
    {
        final String[] items = data.split("\t", -1);

        // 0       1       2                       3                       4               5               6       7       8       9               10
        // #gene1  gene2   strand1(gene/fusion)    strand2(gene/fusion)    breakpoint1     breakpoint2     site1   site2   type    direction1      direction2
        // 11              12              13                      14              15
        // split_reads1    split_reads2    discordant_mates        coverage1       coverage2
        // 16              17                              18                               19
        // confidence      closest_genomic_breakpoint1     closest_genomic_breakpoint2     filters fusion_transcript
        // reading_frame   peptide_sequence        read_identifiers

        final String[] geneNames = new String[] { items[0].replaceAll(",",";"), items[1].replaceAll(",",";") };
        final String[] chromosomes = { items[4].split(":")[0], items[5].split(":")[0] };

        if(!HumanChromosome.contains(chromosomes[SE_START]) || !HumanChromosome.contains(chromosomes[SE_END]))
            return null;

        final int[] positions = { Integer.parseInt(items[4].split(":")[1]), Integer.parseInt(items[5].split(":")[1]) };

        final byte[] orientations = {0, 0};

        String svType = items[8];

        if(svType.contains("deletion"))
        {
            orientations[0] = positions[1] > positions[0] ? POS_ORIENT : NEG_ORIENT;
            orientations[1] = positions[1] > positions[0] ? NEG_ORIENT : POS_ORIENT;
        }
        else if(svType.contains("duplication"))
        {
            orientations[0] = positions[1] > positions[0] ? NEG_ORIENT : POS_ORIENT;
            orientations[1] = positions[1] > positions[0] ? POS_ORIENT : NEG_ORIENT;
        }
        else if(svType.contains("inversion") || svType.contains("translocation"))
        {
            final String[] geneStrands = { items[2], items[3] };
            orientations[0] = geneStrands[0].charAt(0) == '+' ? POS_ORIENT : NEG_ORIENT;
            orientations[1] = geneStrands[1].charAt(0) == '+' ? POS_ORIENT : NEG_ORIENT;
        }

        FusionJunctionType[] junctionTypes = new FusionJunctionType[] { FusionJunctionType.UNKNOWN, FusionJunctionType.UNKNOWN };

        if(items[6].equals("splice-site"))
            junctionTypes[0] = FusionJunctionType.KNOWN;

        if(items[7].equals("splice-site"))
            junctionTypes[1] = FusionJunctionType.KNOWN;

        int splitFragments = Integer.parseInt(items[11]) + Integer.parseInt(items[12]);
        int discordantFragments = Integer.parseInt(items[13]);

        int coverage1 = items[14].equals(".") ? 0 : Integer.parseInt(items[14]);
        int coverage2 = items[15].equals(".") ? 0 : Integer.parseInt(items[15]);
        final int[] coverage = new int[] { coverage1, coverage2 };

        final Map<String, String> otherData = Maps.newHashMap();

        otherData.put("Conf", items[16]);
        otherData.put("Filters", items[19].replaceAll(",",";"));

        return new ExternalFusionData(FUSION_SOURCE_ARRIBA, chromosomes, positions, orientations, junctionTypes, svType, geneNames,
                splitFragments, discordantFragments, coverage, otherData);
    }

    public boolean matches(final FusionData other)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(Chromosomes[se].equals(other.Chromosomes[SE_START]) && Chromosomes[switchIndex(se)].equals(other.Chromosomes[SE_END])
            && JunctionPositions[se] == other.JunctionPositions[SE_START] && JunctionPositions[switchIndex(se)] == other.JunctionPositions[SE_END])
            {
                return true;
            }
        }

        return false;
    }

    public String otherDataStr()
    {
        if(OtherData.isEmpty())
            return "";

        StringJoiner data = new StringJoiner(" ");
        for(Map.Entry<String,String> entry : OtherData.entrySet())
        {
            data.add(entry.getKey() + "=" + entry.getValue());
        }

        return data.toString();
    }


}
