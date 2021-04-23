package com.hartwig.hmftools.isofox.fusion.cohort;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

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

    // consolidated values across samples
    private final List<String> mSampleIds;
    private int mSampleCount;
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

        mSampleCount = 0;
        mTotalFragments = 0;
        mMaxFragments = 0;
    }

    public static FusionCohortData from(final FusionData fusion)
    {
        return new FusionCohortData(fusion.Chromosomes, fusion.JunctionPositions, fusion.JunctionOrientations, fusion.JunctionTypes,
                fusion.SvType, fusion.GeneIds, fusion.GeneNames);
    }

    public boolean matches(final FusionData other)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(!Chromosomes[se].equals(other.Chromosomes[se]))
                return false;

            if(JunctionPositions[se] != other.JunctionPositions[se])
                return false;

            if(JunctionOrientations[se] != other.JunctionOrientations[se])
                return false;
        }

        return true;
    }

    public void addSample(final String sampleId, int fragmentCount)
    {
        if(!mSampleIds.contains(sampleId))
        {
            mSampleIds.add(sampleId);
            ++mSampleCount;
        }

        mMaxFragments = max(mMaxFragments, fragmentCount);
        mTotalFragments += fragmentCount;
    }

    public int sampleCount() { return mSampleCount; }
    public int fragmentCount() { return mTotalFragments; }
    public int maxFragmentCount() { return mMaxFragments; }
    public final List<String> sampleIds() { return mSampleIds; }

    public void setCohortStats(int totalFrags, int maxFrags, int sampleCount)
    {
        mTotalFragments = totalFrags;
        mMaxFragments = maxFrags;
        mSampleCount = sampleCount;
    }

    public static String header()
    {
        StringJoiner header = new StringJoiner(DELIMITER);

        for(int se = SE_START; se <= SE_END; ++se)
        {
            String prefix = se == SE_START ? "Up" : "Down";

            header.add(String.format("GeneId%s", prefix));
            header.add(String.format("GeneName%s", prefix));
            header.add(String.format("Chr%s", prefix));
            header.add(String.format("Pos%s", prefix));
            header.add(String.format("Orient%s", prefix));
            header.add(String.format("JuncType%s", prefix));
        }

        header.add("SVType");
        header.add("SampleCount");
        header.add("TotalFragments");
        header.add("MaxFragments");
        header.add("Samples");
        return header.toString();
    }

    public static String toCsv(final FusionCohortData fusion)
    {
        StringJoiner output = new StringJoiner(DELIMITER);

        for (int se = SE_START; se <= SE_END; ++se)
        {
            output.add(fusion.GeneIds[se]);
            output.add(fusion.GeneNames[se]);
            output.add(fusion.Chromosomes[se]);
            output.add(String.valueOf(fusion.JunctionPositions[se]));
            output.add(String.valueOf(fusion.JunctionOrientations[se]));
            output.add(String.valueOf(fusion.JunctionTypes[se]));
        }

        output.add(fusion.SvType);
        output.add(String.valueOf(fusion.sampleCount()));
        output.add(String.valueOf(fusion.fragmentCount()));
        output.add(String.valueOf(fusion.maxFragmentCount()));
        output.add(appendStrList(fusion.sampleIds(), ';'));

        return output.toString();
    }

    public static FusionCohortData fromCsv(final String data, final Map<String,Integer> fieldIndexMap)
    {
        final String[] items = data.split(DELIMITER);

        final String[] geneIds = new String[] { items[fieldIndexMap.get("GeneIdUp")], items[fieldIndexMap.get("GeneIdDown")] };
        final String[] geneNames = new String[] { items[fieldIndexMap.get("GeneNameUp")], items[fieldIndexMap.get("GeneNameDown")] };

        final String[] chromosomes = new String[] { items[fieldIndexMap.get("ChrUp")], items[fieldIndexMap.get("ChrDown")] };

        final int[] junctionPositions =
                new int[] { Integer.parseInt(items[fieldIndexMap.get("PosUp")]), Integer.parseInt(items[fieldIndexMap.get("PosDown")]) };

        final byte[] junctionOrientations =
                new byte[] { Byte.parseByte(items[fieldIndexMap.get("OrientUp")]), Byte.parseByte(items[fieldIndexMap.get("OrientDown")]) };

        final FusionJunctionType[] junctionTypes = new FusionJunctionType[] {
                FusionJunctionType.valueOf(items[fieldIndexMap.get("JuncTypeUp")]),
                FusionJunctionType.valueOf(items[fieldIndexMap.get("JuncTypeDown")]) };


        final String svType = items[fieldIndexMap.get("SVType")];

        FusionCohortData fusion = new FusionCohortData(
                chromosomes, junctionPositions, junctionOrientations, junctionTypes, svType, geneIds, geneNames);

        int sampleCount = Integer.parseInt(items[fieldIndexMap.get("SampleCount")]);
        int totalFragmentCount = Integer.parseInt(items[fieldIndexMap.get("TotalFragments")]);
        int maxFragmentCount = Integer.parseInt(items[fieldIndexMap.get("MaxFragments")]);

        fusion.setCohortStats(totalFragmentCount, maxFragmentCount, sampleCount);

        return fusion;
    }

}
