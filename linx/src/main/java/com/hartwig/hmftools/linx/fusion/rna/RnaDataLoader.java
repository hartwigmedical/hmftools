package com.hartwig.hmftools.linx.fusion.rna;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;

import java.util.Map;

public class RnaDataLoader
{
    public static RnaFusionData loadRnaFusion(final String data, final Map<String,Integer> fieldIndexMap)
    {
        return fromIsofox(data, fieldIndexMap);
    }

    private static RnaFusionData fromIsofox(final String data, final Map<String,Integer> fieldIndexMap)
    {
        final String[] items = data.split(CSV_DELIM, -1);

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
                sampleId, fusionId, geneIds, geneNames, chromosomes, positions, orientations,
                junctionFrags, discordantFrags, junctionTypes, cohortCount, otherData);
    }
}
