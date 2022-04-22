package com.hartwig.hmftools.isofox.fusion;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.fusion.FusionJunctionType.KNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.FUSION_ID_PREFIX;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.FUSION_NONE;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.fusionId;
import static com.hartwig.hmftools.isofox.fusion.FusionFilterType.NOT_SET;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public class FusionData
{
    public final int Id;
    public final boolean Valid;
    public final String[] Chromosomes;
    public final int[] JunctionPositions;
    public final byte[] JunctionOrientations;
    public final FusionJunctionType[] JunctionTypes;
    public final String SvType;
    public final boolean NonSupp;
    public final String[] GeneIds;
    public final String[] GeneNames;
    public final byte[] Strands;
    public final int[] Coverage;
    public final int[] AnchorDistance;

    public final int TotalFrags;
    public final int SplitFrags;
    public final int RealignedFrags;
    public final int DiscordantFrags;
    public final String RelatedSplicedIds;
    public final String InitReadId;
    public final String[] TransData;

    private final List<Integer> mRelatedFusionIds;

    // annotations from filtering

    private KnownGeneType mKnownFusionType;
    private boolean mHasRelatedKnownSpliceSites;
    private int mCohortFrequency;
    private FusionFilterType mFilter;

    public FusionData(int id, boolean isValid, final String[] chromosomes, final int[] junctionPositions, final byte[] junctionOrientations,
            final FusionJunctionType[] junctionTypes, final String svType, boolean nonSupp,
            final String[] geneIds, final String[] geneNames, final byte[] strands,
            int totalFrags, int splitFrags, int realignedFrags, int discordantFrags, final int[] coverage, final int[] anchorDistance,
            final String[] transData, final String relatedFusionIds, final String initReadId)
    {
        Id = id;
        Valid = isValid;

        Chromosomes = chromosomes;
        JunctionPositions = junctionPositions;
        JunctionOrientations = junctionOrientations;
        JunctionTypes = junctionTypes;

        SvType = svType;
        NonSupp = nonSupp;

        GeneIds = geneIds;
        GeneNames = geneNames;
        Strands = strands;

        TotalFrags = totalFrags;
        SplitFrags = splitFrags;
        RealignedFrags = realignedFrags;
        DiscordantFrags = discordantFrags;
        Coverage = coverage;
        AnchorDistance = anchorDistance;
        InitReadId = initReadId;
        TransData = transData;
        RelatedSplicedIds = relatedFusionIds;

        // mRawData = null;
        mKnownFusionType = KnownGeneType.OTHER;
        mRelatedFusionIds = Lists.newArrayList();

        if(!relatedFusionIds.isEmpty() && !relatedFusionIds.equals(FUSION_NONE))
        {
            Arrays.stream(relatedFusionIds.split(ITEM_DELIM, -1))
                    .mapToInt(x -> Integer.parseInt(x.replaceAll(FUSION_ID_PREFIX, "")))
                    .forEach(x -> mRelatedFusionIds.add(x));
        }

        mHasRelatedKnownSpliceSites = false;
        mCohortFrequency = 0;
        mFilter = NOT_SET;
    }

    public String name() { return String.format("%s_%s", GeneNames[SE_START], GeneNames[SE_END]); }

    public List<Integer> relatedFusionIds() { return mRelatedFusionIds; }

    public final KnownGeneType getKnownFusionType() { return mKnownFusionType; }
    public void setKnownFusionType(KnownGeneType type) { mKnownFusionType = type; }

    public boolean isRelated(final FusionData other)
    {
        if(!mRelatedFusionIds.contains(other.Id))
            return false;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(GeneNames[se].isEmpty() || !other.GeneNames[se].equals(GeneNames[se]))
                return false;
        }

        return true;
    }

    public boolean hasKnownSpliceSites() { return JunctionTypes[SE_START] == KNOWN && JunctionTypes[SE_END] == KNOWN; }

    public void setHasRelatedKnownSpliceSites() { mHasRelatedKnownSpliceSites = true; }
    public boolean hasRelatedKnownSpliceSites() { return mHasRelatedKnownSpliceSites; }

    public void setCohortFrequency(int count) { mCohortFrequency = count; }
    public int cohortFrequency() { return mCohortFrequency; }

    public void setFilter(FusionFilterType type) { mFilter = type; }
    public FusionFilterType getFilter() { return mFilter; }

    public String toString()
    {
        return String.format("%d: chr(%s-%s) junc(%d-%d %s) genes(%s-%s) frags(%d)",
                Id, Chromosomes[SE_START], Chromosomes[SE_END], JunctionPositions[SE_START], JunctionPositions[SE_END],
                SvType, GeneNames[SE_START], GeneNames[SE_END], TotalFrags);
    }

    public static final String FLD_FUSION_ID = "FusionId";
    private static final String FLD_VALID = "Valid";
    public static final String FLD_CHR = "Chr";
    public static final String FLD_POS = "Pos";
    public static final String FLD_ORIENT = "Orient";
    public static final String FLD_JUNC_TYPE = "JuncType";
    public static final String FLD_STRAND = "Strand";
    public static final String FLD_SV_TYPE = "SVType";
    public static final String FLD_COVERAGE = "Coverage";
    public static final String FLD_MAX_ANCHOR = "MaxAnchorLength";
    public static final String FLD_TOTAL_FRAGS = "TotalFragments";
    public static final String FLD_SPLIT_FRAGS = "SplitFrags";
    public static final String FLD_REALIGN_FLAGS = "RealignedFrags";
    public static final String FLD_DISCORD_FRAGS = "DiscordantFrags";
    public static final String FLD_NON_SUPP = "NonSupp";
    public static final String FLD_TRANS_DATA = "TransData";
    public static final String FLD_REL_SPLICED_IDS = "RelatedSplicedIds";
    public static final String FLD_INIT_READ_ID = "InitReadId";

    // post-filtering fields
    public static final String FLD_FILTER = "Filter";
    public static final String FLD_COHORT_COUNT = "CohortCount";
    public static final String FLD_KNOWN_TYPE = "KnownFusionType";

    // more verbose discovery fields
    private static final String FLD_HOM_OFFSET = "HomologyOffset";

    public static String csvHeader(boolean isFiltered)
    {
        StringJoiner sj = new StringJoiner(DELIMITER);
        sj.add(FLD_FUSION_ID);
        sj.add(FLD_VALID);

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            sj.add(formStreamField(FLD_CHR, fs));
            sj.add(formStreamField(FLD_POS, fs));
            sj.add(formStreamField(FLD_ORIENT, fs));
            sj.add(formStreamField(FLD_GENE_ID, fs));
            sj.add(formStreamField(FLD_GENE_NAME, fs));
            sj.add(formStreamField(FLD_STRAND, fs));
            sj.add(formStreamField(FLD_JUNC_TYPE, fs));
            sj.add(formStreamField(FLD_COVERAGE, fs));
            sj.add(formStreamField(FLD_MAX_ANCHOR, fs));

            if(!isFiltered)
            {
                sj.add(formStreamField(FLD_TRANS_DATA, fs));
            }
        }

        sj.add(FLD_SV_TYPE);
        sj.add(FLD_TOTAL_FRAGS);
        sj.add(FLD_SPLIT_FRAGS);
        sj.add(FLD_REALIGN_FLAGS);
        sj.add(FLD_DISCORD_FRAGS);
        sj.add(FLD_REL_SPLICED_IDS);
        sj.add(FLD_COHORT_COUNT);

        if(isFiltered)
        {
            sj.add(FLD_FILTER);
            sj.add(FLD_KNOWN_TYPE);
        }
        else
        {
            sj.add(FLD_NON_SUPP);
            sj.add(FLD_INIT_READ_ID);
        }

        return sj.toString();
    }

    public String toCsv(boolean isFiltered)
    {
        StringJoiner sj = new StringJoiner(DELIMITER);

        sj.add(fusionId(Id));
        sj.add(String.valueOf(Valid));

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            sj.add(Chromosomes[fs]);
            sj.add(String.valueOf(JunctionPositions[fs]));
            sj.add(String.valueOf(JunctionOrientations[fs]));
            sj.add(GeneIds[fs]);
            sj.add(GeneNames[fs]);
            sj.add(String.valueOf(Strands[fs]));
            sj.add(String.valueOf(JunctionTypes[fs]));
            sj.add(String.valueOf(Coverage[fs]));
            sj.add(String.valueOf(AnchorDistance[fs]));

            if(!isFiltered)
            {
                sj.add(TransData[fs]);
            }
        }

        sj.add(SvType);
        sj.add(String.valueOf(TotalFrags));
        sj.add(String.valueOf(SplitFrags));
        sj.add(String.valueOf(RealignedFrags));
        sj.add(String.valueOf(DiscordantFrags));
        sj.add(RelatedSplicedIds);
        sj.add(String.valueOf(mCohortFrequency));

        // could write an additional fields - eg RelatedProxIds, HomologyOffset were previously written

        if(isFiltered)
        {
            sj.add(mFilter.toString());
            sj.add(getKnownFusionType().toString());
        }
        else
        {
            sj.add(String.valueOf(NonSupp));
            sj.add(InitReadId);
        }

        return sj.toString();
    }

    public static String formStreamField(final String field, final int stream)
    {
        return field + (stream == FS_UP ? "Up" : "Down");
    }

    public static List<FusionData> loadFromFile(final Path filename)
    {
        List<FusionData> fusions = Lists.newArrayList();

        try
        {
            final List<String> lines = Files.readAllLines(filename);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int fusionIdIndex = fieldsIndexMap.get(FLD_FUSION_ID);
            int validIndex = fieldsIndexMap.get(FLD_VALID);
            int chrUpIndex = fieldsIndexMap.get(formStreamField(FLD_CHR, FS_UP));
            int chrDownIndex = fieldsIndexMap.get(formStreamField(FLD_CHR, FS_DOWN));
            int posUpIndex = fieldsIndexMap.get(formStreamField(FLD_POS, FS_UP));
            int posDownIndex = fieldsIndexMap.get(formStreamField(FLD_POS, FS_DOWN));
            int orientUpIndex = fieldsIndexMap.get(formStreamField(FLD_ORIENT, FS_UP));
            int orientDownIndex = fieldsIndexMap.get(formStreamField(FLD_ORIENT, FS_DOWN));
            int geneIdUpIndex = fieldsIndexMap.get(formStreamField(FLD_GENE_ID, FS_UP));
            int geneIdDownIndex = fieldsIndexMap.get(formStreamField(FLD_GENE_ID, FS_DOWN));
            int geneNameUpIndex = fieldsIndexMap.get(formStreamField(FLD_GENE_NAME, FS_UP));
            int geneNameDownIndex = fieldsIndexMap.get(formStreamField(FLD_GENE_NAME, FS_DOWN));
            int svTypeIndex = fieldsIndexMap.get(FLD_SV_TYPE);
            int strandUpIndex = fieldsIndexMap.get(formStreamField(FLD_STRAND, FS_UP));
            int strandDownIndex = fieldsIndexMap.get(formStreamField(FLD_STRAND, FS_DOWN));
            int juncTypeUpIndex = fieldsIndexMap.get(formStreamField(FLD_JUNC_TYPE, FS_UP));
            int juncTypeDownIndex = fieldsIndexMap.get(formStreamField(FLD_JUNC_TYPE, FS_DOWN));

            int totalFragsIndex = fieldsIndexMap.get(FLD_TOTAL_FRAGS);
            int splitFragsIndex = fieldsIndexMap.get(FLD_SPLIT_FRAGS);
            int realignedFragsIndex = fieldsIndexMap.get(FLD_REALIGN_FLAGS);
            int discordantFragsIndex = fieldsIndexMap.get(FLD_DISCORD_FRAGS);

            int ancDistUpIndex = fieldsIndexMap.get(formStreamField(FLD_MAX_ANCHOR, FS_UP));
            int ancDistDownIndex = fieldsIndexMap.get(formStreamField(FLD_MAX_ANCHOR, FS_DOWN));
            int covUpIndex = fieldsIndexMap.get(formStreamField(FLD_COVERAGE, FS_UP));
            int covDownIndex = fieldsIndexMap.get(formStreamField(FLD_COVERAGE, FS_DOWN));
            int relatedFusionIdsIndex = fieldsIndexMap.get(FLD_REL_SPLICED_IDS);

            // optional fields - raw verbose
            Integer transUpIndex = fieldsIndexMap.get(formStreamField(FLD_TRANS_DATA, FS_UP));
            Integer transDownIndex = fieldsIndexMap.get(formStreamField(FLD_TRANS_DATA, FS_DOWN));
            Integer initReadIndex = fieldsIndexMap.get(FLD_INIT_READ_ID);
            Integer nonSuppIndex = fieldsIndexMap.get(FLD_NON_SUPP);

            // optional fields - filtered
            Integer filterIndex = fieldsIndexMap.get(FLD_FILTER);
            Integer cohortFreqIndex = fieldsIndexMap.get(FLD_COHORT_COUNT);

            for(String data : lines)
            {
                final String[] values = data.split(DELIMITER, -1);

                int fusionId = Integer.parseInt(values[fusionIdIndex].replaceAll(FUSION_ID_PREFIX, ""));
                boolean isValid = Boolean.parseBoolean(values[validIndex]);

                final String[] chromosomes = new String[] { values[chrUpIndex], values[chrDownIndex] };

                final int[] junctionPositions =
                        new int[] { Integer.parseInt(values[posUpIndex]), Integer.parseInt(values[posDownIndex]) };

                final byte[] junctionOrientations =
                        new byte[] { Byte.parseByte(values[orientUpIndex]), Byte.parseByte(values[orientDownIndex]) };

                final FusionJunctionType[] junctionTypes = new FusionJunctionType[] {
                        FusionJunctionType.valueOf(values[juncTypeUpIndex]), FusionJunctionType.valueOf(values[juncTypeDownIndex]) };

                final String[] geneIds = new String[] { values[geneIdUpIndex], values[geneIdDownIndex] };
                final String[] geneNames = new String[] { values[geneNameUpIndex], values[geneNameDownIndex] };

                final byte[] strands =
                        new byte[] { Byte.parseByte(values[strandUpIndex]), Byte.parseByte(values[strandDownIndex]) };

                final String svType = values[svTypeIndex];
                boolean nonSupp = nonSuppIndex != null ? Boolean.parseBoolean(values[nonSuppIndex]) : false;

                final int[] coverage = new int[] { Integer.parseInt(values[covUpIndex]), Integer.parseInt(values[covDownIndex]) };

                final int[] anchorDistance = new int[] { Integer.parseInt(values[ancDistUpIndex]), Integer.parseInt(values[ancDistDownIndex]) };

                final String[] transData = transUpIndex != null && transDownIndex != null ?
                        new String[] { values[transUpIndex], values[transDownIndex] } : new String[] { "", "" };

                String initReadId = initReadIndex != null ? values[initReadIndex] : "";

                FusionData fusion = new FusionData(
                        fusionId, isValid, chromosomes, junctionPositions, junctionOrientations, junctionTypes,
                        svType, nonSupp, geneIds, geneNames, strands,
                        Integer.parseInt(values[totalFragsIndex]), Integer.parseInt(values[splitFragsIndex]),
                        Integer.parseInt(values[realignedFragsIndex]), Integer.parseInt(values[discordantFragsIndex]),
                        coverage, anchorDistance, transData, values[relatedFusionIdsIndex], initReadId);

                // optional fields
                if(filterIndex != null)
                    fusion.setFilter(FusionFilterType.valueOf(values[filterIndex]));

                if(cohortFreqIndex != null)
                    fusion.setCohortFrequency(Integer.parseInt(values[cohortFreqIndex]));

                fusions.add(fusion);
            }

            return fusions;
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load fusion file({}): {}", filename.toString(), e.toString());
        }

        return fusions;
    }

    public int length()
    {
        return SvType.equals(BND.toString()) ? 0 : JunctionPositions[SE_END] - JunctionPositions[SE_START];
    }

    public double alleleFrequency()
    {
        double coverage = max(Coverage[SE_START], Coverage[SE_END]);
        return coverage> 0 ? (SplitFrags + RealignedFrags) / coverage : 0;
    }

    public int totalFragments() { return SplitFrags + RealignedFrags + DiscordantFrags; }
    public int supportingFragments() { return SplitFrags + RealignedFrags; }
}
