package com.hartwig.hmftools.isofox.fusion.cohort;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.Map;

import com.hartwig.hmftools.isofox.fusion.FusionJunctionType;

public class FusionData
{
    public final String[] Chromosomes;
    public final int[] JunctionPositions;
    public final byte[] JunctionOrientations;
    public final FusionJunctionType[] JunctionTypes;
    public final String SvType;
    public final String[] GeneIds;
    public final String[] GeneNames;
    public final int[] Coverage;
    public final int SplitFrags;
    public final int RealignedFrags;
    public final int DiscordantFrags;

    public FusionData(final String[] chromosomes, final int[] junctionPositions, final byte[] junctionOrientations,
            final FusionJunctionType[] junctionTypes, final String svType, final String[] geneIds, final String[] geneNames,
            int splitFrags, int realignedFrags, int discordantFrags, final int[] coverage)
    {
        Chromosomes = chromosomes;
        JunctionPositions = junctionPositions;
        JunctionOrientations = junctionOrientations;
        JunctionTypes = junctionTypes;
        SvType = svType;
        GeneIds = geneIds;
        GeneNames = geneNames;
        SplitFrags = splitFrags;
        RealignedFrags = realignedFrags;
        DiscordantFrags = discordantFrags;
        Coverage = coverage;
    }

    public static FusionData fromCsv(final String data, final Map<String,Integer> fieldIndexMap)
    {
        final String[] items = data.split(DELIMITER);

        final String[] chromosomes = new String[] { items[fieldIndexMap.get("ChrUp")], items[fieldIndexMap.get("ChrDown")] };

        final int[] junctionPositions =
                new int[] { Integer.parseInt(items[fieldIndexMap.get("PosUp")]), Integer.parseInt(items[fieldIndexMap.get("PosDown")]) };

        final byte[] junctionOrientations =
                new byte[] { Byte.parseByte(items[fieldIndexMap.get("OrientUp")]), Byte.parseByte(items[fieldIndexMap.get("OrientDown")]) };

        final FusionJunctionType[] junctionTypes = new FusionJunctionType[] {
                FusionJunctionType.valueOf(items[fieldIndexMap.get("JuncTypeUp")]),
                FusionJunctionType.valueOf(items[fieldIndexMap.get("JuncTypeDown")]) };

        final String[] geneIds = new String[] { items[fieldIndexMap.get("GeneIdUp")], items[fieldIndexMap.get("GeneIdDown")] };
        final String[] geneNames = new String[] { items[fieldIndexMap.get("GeneNameUp")], items[fieldIndexMap.get("GeneNameDown")] };

        final String svType = items[fieldIndexMap.get("SVType")];

        final int[] coverage = new int[] {
                Integer.parseInt(items[fieldIndexMap.get("CoverageUp")]), Integer.parseInt(items[fieldIndexMap.get("CoverageDown")]) };

        FusionData fusion = new FusionData(
                chromosomes, junctionPositions, junctionOrientations, junctionTypes, svType, geneIds, geneNames,
                Integer.parseInt(items[fieldIndexMap.get("SplitFrags")]),
                Integer.parseInt(items[fieldIndexMap.get("RealignedFrags")]),
                Integer.parseInt(items[fieldIndexMap.get("DiscordantFrags")]),
                coverage);

        return fusion;
    }

    public int length()
    {
        return SvType.equals(BND.toString()) ? 0 : JunctionPositions[SE_END] - JunctionPositions[SE_START];
    }

    public double alleleFrequency()
    {
        int coverage = max(Coverage[SE_START], Coverage[SE_END]);
        return coverage> 0 ? (SplitFrags + RealignedFrags) / coverage : 0;
    }
}
