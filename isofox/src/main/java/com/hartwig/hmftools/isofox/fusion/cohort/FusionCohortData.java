package com.hartwig.hmftools.isofox.fusion.cohort;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_CHR;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_COHORT_COUNT;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_JUNC_TYPE;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_ORIENT;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_POS;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_SV_TYPE;
import static com.hartwig.hmftools.isofox.fusion.FusionData.formStreamField;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.isofox.fusion.FusionData;
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

    public static String header()
    {
        StringJoiner header = new StringJoiner(DELIMITER);

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            header.add(formStreamField(FLD_GENE_ID, fs));
            header.add(formStreamField(FLD_GENE_NAME, fs));
            header.add(formStreamField(FLD_CHR, fs));
            header.add(formStreamField(FLD_POS, fs));
            header.add(formStreamField(FLD_ORIENT, fs));
            header.add(formStreamField(FLD_JUNC_TYPE, fs));
        }

        header.add(FLD_SV_TYPE);
        header.add(FLD_COHORT_COUNT);
        header.add("TotalFragments");
        header.add("MaxFragments");
        header.add("Samples");
        return header.toString();
    }

    public static String toCsv(final FusionCohortData fusion)
    {
        StringJoiner output = new StringJoiner(DELIMITER);

        for (int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            output.add(fusion.GeneIds[fs]);
            output.add(fusion.GeneNames[fs]);
            output.add(fusion.Chromosomes[fs]);
            output.add(String.valueOf(fusion.JunctionPositions[fs]));
            output.add(String.valueOf(fusion.JunctionOrientations[fs]));
            output.add(String.valueOf(fusion.JunctionTypes[fs]));
        }

        output.add(fusion.SvType);
        output.add(String.valueOf(fusion.sampleCount()));
        output.add(String.valueOf(fusion.fragmentCount()));
        output.add(String.valueOf(fusion.maxFragmentCount()));

        StringJoiner sj = new StringJoiner(";");
        fusion.sampleIds().forEach(x -> sj.add(x));
        output.add(sj.toString());

        return output.toString();
    }
}
