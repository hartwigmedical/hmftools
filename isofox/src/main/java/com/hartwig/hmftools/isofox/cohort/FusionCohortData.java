package com.hartwig.hmftools.isofox.cohort;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.isofox.fusion.FusionJunctionType;

import org.apache.commons.compress.utils.Lists;

public class FusionCohortData
{
    public final String[] Chromosomes;
    public final int[] JunctionPositions;
    public final byte[] JunctionOrientations;
    public final FusionJunctionType[] JunctionTypes;
    public final String SvType;

    public final String[] GeneIds;
    public final String[] GeneNames;

    private final List<String> mSampleIds;
    private int mTotalFragments;
    private int mMaxFragments;

    public FusionCohortData(final String[] chromosomes, final int[] junctionPositions, final byte[] junctionOrientations,
            final FusionJunctionType[] junctionTypes, final String svType, final String[] geneIds, final String[] geneNames)
    {
        mSampleIds = Lists.newArrayList();
        Chromosomes = chromosomes;
        JunctionPositions = junctionPositions;
        JunctionOrientations = junctionOrientations;
        JunctionTypes = junctionTypes;
        SvType = svType;
        GeneIds = geneIds;
        GeneNames = geneNames;

        mTotalFragments = 0;
        mMaxFragments = 0;
    }

    public boolean matches(final FusionCohortData other)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(!Chromosomes[se].equals(other.Chromosomes[se]))
                return false;

            if(JunctionPositions[se] != other.JunctionPositions[se])
                return false;
        }

        return true;
    }

    public void addSample(final String sampleId, int fragmentCount)
    {
        if(!mSampleIds.contains(sampleId))
            mSampleIds.add(sampleId);

        mMaxFragments = max(mMaxFragments, fragmentCount);
        mTotalFragments += fragmentCount;
    }

    public int sampleCount() { return mSampleIds.size(); }
    public int fragmentCount() { return mTotalFragments; }
    public int maxFragmentCount() { return mMaxFragments; }
    public final List<String> sampleIds() { return mSampleIds; }

    public static FusionCohortData fromCsv(final String data, final Map<String,Integer> fieldIndexMap, final String sampleId)
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

        FusionCohortData fusion = new FusionCohortData(
                chromosomes, junctionPositions, junctionOrientations, junctionTypes, svType, geneIds, geneNames);

        fusion.addSample(sampleId, Integer.parseInt(items[fieldIndexMap.get("SplitFrags")]));

        return fusion;
    }
}
