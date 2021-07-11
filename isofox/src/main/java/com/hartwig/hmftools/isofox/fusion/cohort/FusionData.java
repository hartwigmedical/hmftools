package com.hartwig.hmftools.isofox.fusion.cohort;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.isofox.fusion.FusionJunctionType.KNOWN;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.FUSION_ID_PREFIX;
import static com.hartwig.hmftools.isofox.fusion.FusionReadData.FUSION_NONE;
import static com.hartwig.hmftools.isofox.fusion.cohort.FusionFilterType.NOT_SET;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.isofox.fusion.FusionJunctionType;

import org.apache.commons.compress.utils.Lists;

public class FusionData
{
    public final int Id;
    public final String[] Chromosomes;
    public final int[] JunctionPositions;
    public final byte[] JunctionOrientations;
    public final FusionJunctionType[] JunctionTypes;
    public final String SvType;
    public final String[] GeneIds;
    public final String[] GeneNames;
    public final int[] Coverage;
    public final int[] AnchorDistance;
    public final int SplitFrags;
    public final int RealignedFrags;
    public final int DiscordantFrags;
    public final int CohortCount;

    private String mRawData;
    private KnownGeneType mKnownFusionType;
    private final List<Integer> mRelatedFusionIds;
    private boolean mHasRelatedKnownSpliceSites;
    private int mCohortFrequency;
    private FusionFilterType mFilter;

    public FusionData(int id, final String[] chromosomes, final int[] junctionPositions, final byte[] junctionOrientations,
            final FusionJunctionType[] junctionTypes, final String svType, final String[] geneIds, final String[] geneNames,
            int splitFrags, int realignedFrags, int discordantFrags, final int[] coverage, final int[] anchorDistance, int cohortCount)
    {
        Id = id;
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
        AnchorDistance = anchorDistance;
        CohortCount = cohortCount;

        mRawData = null;
        mKnownFusionType = KnownGeneType.OTHER;
        mRelatedFusionIds = Lists.newArrayList();
        mHasRelatedKnownSpliceSites = false;
        mCohortFrequency = 0;
        mFilter = NOT_SET;
    }

    public String name() { return String.format("%s_%s", GeneNames[SE_START], GeneNames[SE_END]); }

    public void cacheCsvData(final String data) { mRawData = data; }
    public String rawData() { return mRawData; }
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

    public static FusionData fromCsv(final String data, final Map<String,Integer> fieldIndexMap)
    {
        final String[] items = data.split(DELIMITER);

        int fusionId = Integer.parseInt(items[fieldIndexMap.get("FusionId")].replaceAll(FUSION_ID_PREFIX, ""));

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

        final int[] anchorDistance = new int[] {
                Integer.parseInt(items[fieldIndexMap.get("MaxAnchorLengthUp")]), Integer.parseInt(items[fieldIndexMap.get("MaxAnchorLengthDown")]) };

        final int cohortCount = fieldIndexMap.containsKey("CohortCount") ?
                Integer.parseInt(items[fieldIndexMap.get("CohortCount")]) : 0;

        FusionData fusion = new FusionData(
                fusionId, chromosomes, junctionPositions, junctionOrientations, junctionTypes, svType, geneIds, geneNames,
                Integer.parseInt(items[fieldIndexMap.get("SplitFrags")]),
                Integer.parseInt(items[fieldIndexMap.get("RealignedFrags")]),
                Integer.parseInt(items[fieldIndexMap.get("DiscordantFrags")]),
                coverage, anchorDistance, cohortCount);

        if(fieldIndexMap.containsKey("Filter"))
            fusion.setFilter(FusionFilterType.valueOf(items[fieldIndexMap.get("Filter")]));

        final String[] relatedFusionIdStr = items[fieldIndexMap.get("RelatedSplicedIds")].split(ITEM_DELIM, -1);

        Arrays.stream(relatedFusionIdStr)
                .filter(x -> !x.equals(FUSION_NONE))
                .mapToInt(x -> Integer.parseInt(x.replaceAll(FUSION_ID_PREFIX, "")))
                .forEach(x -> fusion.relatedFusionIds().add(x));

        return fusion;
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
