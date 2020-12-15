package com.hartwig.hmftools.linx.fusion.rna;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.fusion.rna.RnaJunctionType.KNOWN;

import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public class RnaDataLoader
{
    public static final String RNA_FUSION_SOURCE_ISOFOX = "ISOFOX";
    public static final String RNA_FUSION_SOURCE_MATCHED = "MATCHED";
    public static final String RNA_FUSION_SOURCE_ARRIBA = "ARRIBA";
    public static final String RNA_FUSION_SOURCE_STARFUSION = "STARFUSION";

    public static String getRnaSourceDelimiter(final String source)
    {
        return source.equals(RNA_FUSION_SOURCE_ARRIBA) ? "," : ","; // arriba now reformatted in R
    }

    public static boolean isIsofoxFusion(final String source)
    {
        return source.equals(RNA_FUSION_SOURCE_ISOFOX) || source.equals("MATCH_EXT_UNFILTERED") || source.equals("MATCH");
    }

    public static RnaFusionData loadRnaFusion(final String source, final int index, final String data, final Map<String,Integer> fieldIndexMap)
    {
        if(source.equals(RNA_FUSION_SOURCE_ISOFOX))
            return fromIsofox(data, fieldIndexMap);
        if(source.equals(RNA_FUSION_SOURCE_MATCHED))
            return fromMatched(data, fieldIndexMap);
        else if(source.equals(RNA_FUSION_SOURCE_ARRIBA))
            return fromArriba(index, data, fieldIndexMap);
        else if(source.equals(RNA_FUSION_SOURCE_STARFUSION))
            return fromStarFusion(index, data, fieldIndexMap);
        else
            return null;
    }

    private static RnaFusionData fromIsofox(final String data, final Map<String,Integer> fieldIndexMap)
    {
        final String[] items = data.split(getRnaSourceDelimiter(RNA_FUSION_SOURCE_ISOFOX), -1);

        final String sampleId = items[fieldIndexMap.get("SampleId")];

        String fusionId = items[fieldIndexMap.get("FusionId")];

        final String[] chromosomes = new String[] { items[fieldIndexMap.get("ChrUp")], items[fieldIndexMap.get("ChrDown")] };

        final int[] positions =
                new int[] { Integer.parseInt(items[fieldIndexMap.get("PosUp")]), Integer.parseInt(items[fieldIndexMap.get("PosDown")]) };

        final byte[] orientations =
                new byte[] { Byte.parseByte(items[fieldIndexMap.get("OrientUp")]), Byte.parseByte(items[fieldIndexMap.get("OrientDown")]) };

        final RnaJunctionType[] junctionTypes = new RnaJunctionType[] {
                RnaJunctionType.valueOf(items[fieldIndexMap.get("JuncTypeUp")]),
                RnaJunctionType.valueOf(items[fieldIndexMap.get("JuncTypeDown")]) };

        final String[] geneIds = new String[] { items[fieldIndexMap.get("GeneIdUp")], items[fieldIndexMap.get("GeneIdDown")] };
        final String[] geneNames = new String[] { items[fieldIndexMap.get("GeneNameUp")], items[fieldIndexMap.get("GeneNameDown")] };

        int junctionFrags = Integer.parseInt(items[fieldIndexMap.get("SplitFrags")]) + Integer.parseInt(items[fieldIndexMap.get("RealignedFrags")]);
        int discordantFrags = Integer.parseInt(items[fieldIndexMap.get("DiscordantFrags")]);

        final int cohortCount = fieldIndexMap.containsKey("CohortCount") ? Integer.parseInt(items[fieldIndexMap.get("CohortCount")]) : 0;
        final String otherData = "";

        return new RnaFusionData(
                sampleId, RNA_FUSION_SOURCE_ISOFOX, fusionId, geneIds, geneNames, chromosomes, positions, orientations,
                junctionFrags, discordantFrags, junctionTypes, cohortCount, otherData);
    }

    private static RnaFusionData fromMatched(final String data, final Map<String,Integer> fieldIndexMap)
    {
        final String[] items = data.split(getRnaSourceDelimiter(RNA_FUSION_SOURCE_MATCHED), -1);

        final String matchType = items[fieldIndexMap.get("MatchType")];

        final String sampleId = items[fieldIndexMap.get("SampleId")];

        String fusionId = items[fieldIndexMap.get("FusionId")];

        final String[] chromosomes = new String[] { items[fieldIndexMap.get("ChrUp")], items[fieldIndexMap.get("ChrDown")] };

        final int[] positions =
                new int[] { Integer.parseInt(items[fieldIndexMap.get("PosUp")]), Integer.parseInt(items[fieldIndexMap.get("PosDown")]) };

        final byte[] orientations =
                new byte[] { Byte.parseByte(items[fieldIndexMap.get("OrientUp")]), Byte.parseByte(items[fieldIndexMap.get("OrientDown")]) };

        final RnaJunctionType[] junctionTypes = new RnaJunctionType[] {
                RnaJunctionType.valueOf(items[fieldIndexMap.get("JuncTypeUp")]),
                RnaJunctionType.valueOf(items[fieldIndexMap.get("JuncTypeDown")]) };

        final String[] geneIds = new String[] { items[fieldIndexMap.get("GeneIdUp")], items[fieldIndexMap.get("GeneIdDown")] };

        final String[] geneNames = new String[] { items[fieldIndexMap.get("GeneNameUp")], items[fieldIndexMap.get("GeneNameDown")] };

        int junctionFrags = 0;
        int discordantFrags = 0;

        if(matchType.equals("MATCH") || matchType.equals("FILTERED_IN_ARRIBA") || matchType.equals("ISOFOX_ONLY"))
        {
            junctionFrags = Integer.parseInt(items[fieldIndexMap.get("JuncFrags")]);
            discordantFrags = Integer.parseInt(items[fieldIndexMap.get("DiscFrags")]);
        }
        else
        {
            junctionFrags = Integer.parseInt(items[fieldIndexMap.get("ExtJuncFrags")]);
            discordantFrags = Integer.parseInt(items[fieldIndexMap.get("ExtDiscFrags")]);
        }

        final int cohortCount = fieldIndexMap.containsKey("CohortCount") ? Integer.parseInt(items[fieldIndexMap.get("CohortCount")]) : 0;
        final String otherData = items[fieldIndexMap.get("OtherData")];

        return new RnaFusionData(
                sampleId, matchType, fusionId, geneIds, geneNames, chromosomes, positions, orientations,
                junctionFrags, discordantFrags, junctionTypes, cohortCount, otherData);
    }

    private static RnaFusionData fromArriba(int index, final String data, final Map<String,Integer> fieldIndexMap)
    {
        // SampleId,GeneNameUp,GeneNameDown,ChrUp,ChrDown,PosUp,PosDown,JuncTypeUp,JuncTypeDown,JunctionFrags,DiscordantFrags,CoverageUp,CoverageDown

        final String[] items = data.split(getRnaSourceDelimiter(RNA_FUSION_SOURCE_ARRIBA), -1);

        final String sampleId = items[fieldIndexMap.get("SampleId")];

        final String[] geneIds = new String[] {"", ""};

        final String[] geneNames = new String[] { items[fieldIndexMap.get("GeneNameUp")], items[fieldIndexMap.get("GeneNameDown")] };

        final String[] chromosomes = new String[] { items[fieldIndexMap.get("ChrUp")], items[fieldIndexMap.get("ChrDown")] };

        if(!HumanChromosome.contains(chromosomes[SE_START]) || !HumanChromosome.contains(chromosomes[SE_END]))
            return null;

        final int[] positions = new int[] {
                Integer.parseInt(items[fieldIndexMap.get("PosUp")]), Integer.parseInt(items[fieldIndexMap.get("PosDown")]) };

        final byte[] orientations =
                new byte[] { Byte.parseByte(items[fieldIndexMap.get("OrientUp")]), Byte.parseByte(items[fieldIndexMap.get("OrientDown")]) };

        int splitFragments = Integer.parseInt(items[fieldIndexMap.get("JunctionFrags")]);

        int discordantFragments = Integer.parseInt(items[fieldIndexMap.get("DiscordantFrags")]);

        final String[] spliceSites = new String[] { items[fieldIndexMap.get("JuncTypeUp")], items[fieldIndexMap.get("JuncTypeDown")] };

        final RnaJunctionType[] junctionTypes = new RnaJunctionType[] { RnaJunctionType.NOT_SET, RnaJunctionType.NOT_SET };

        if(spliceSites[FS_UP].equals("splice-site"))
            junctionTypes[FS_UP] = KNOWN;

        if(spliceSites[FS_DOWN].equals("splice-site"))
            junctionTypes[FS_DOWN] = KNOWN;

        final String otherData = "";

        return new RnaFusionData(
                sampleId, RNA_FUSION_SOURCE_ARRIBA, String.valueOf(index), geneIds, geneNames, chromosomes, positions, orientations,
                splitFragments, discordantFragments, junctionTypes, 0, otherData);
    }

    private static RnaFusionData fromStarFusion(int index, final String data, final Map<String,Integer> fieldIndexMap)
    {
        // SampleId,FusionName,JunctionReadCount,SpanningFragCount,SpliceType,GeneNameUp,GeneIdUp,ChrUp,PosUp,OrientUp,GeneNameDown,GeneIdDown,
        // ChrDown,PosDown,OrientDown,JunctionReads,SpanningFrags,
        // LargeAnchorSupport,FFPM,LeftBreakDinuc,LeftBreakEntropy,RightBreakDinuc,RightBreakEntropy,annots
        final String[] items = data.split(getRnaSourceDelimiter(RNA_FUSION_SOURCE_STARFUSION), -1);

        final String sampleId = items[fieldIndexMap.get("SampleId")];

        final String[] geneIds = new String[] {"", ""};

        final String[] geneNames = new String[] { items[fieldIndexMap.get("GeneNameUp")], items[fieldIndexMap.get("GeneNameDown")] };

        final String[] chromosomes = new String[] { items[fieldIndexMap.get("ChrUp")], items[fieldIndexMap.get("ChrDown")] };

        final int[] positions = new int[] {
                Integer.parseInt(items[fieldIndexMap.get("PosUp")]), Integer.parseInt(items[fieldIndexMap.get("PosDown")]) };

        final byte[] orientations =
                new byte[] { Byte.parseByte(items[fieldIndexMap.get("OrientUp")]), Byte.parseByte(items[fieldIndexMap.get("OrientDown")]) };

        final RnaJunctionType[] junctionTypes = new RnaJunctionType[] { RnaJunctionType.NOT_SET, RnaJunctionType.NOT_SET };

        if(items[fieldIndexMap.get("SpliceType")].equals("ONLY_REF_SPLICE"))
        {
            junctionTypes[FS_UP] = KNOWN;
            junctionTypes[FS_DOWN] = KNOWN;
        }

        int junctionCount = Integer.parseInt(items[fieldIndexMap.get("JunctionReadCount")]);
        int spanningCount = Integer.parseInt(items[fieldIndexMap.get("SpanningFragCount")]);

        return new RnaFusionData(
                sampleId, RNA_FUSION_SOURCE_STARFUSION, String.valueOf(index), geneIds, geneNames, chromosomes, positions, orientations,
                junctionCount, spanningCount, junctionTypes, 0, "");
    }
}
