package com.hartwig.hmftools.cup.sample;

import com.hartwig.hmftools.common.purple.gender.Gender;

public class SampleTraitsData
{
    public final String SampleId;
    public final Gender GenderType;
    public final boolean HasWGD;
    public final double Purity;
    public final double Ploidy;
    public final double SnvMbPerMb;
    public final double MsiScore;

    public SampleTraitsData(final String sampleId, final Gender genderType, final boolean hasWGD, final double purity, final double ploidy,
            final double snvMbPerMb, final double msiScore)
    {
        SampleId = sampleId;
        GenderType = genderType;
        HasWGD = hasWGD;
        Purity = purity;
        Ploidy = ploidy;
        SnvMbPerMb = snvMbPerMb;
        MsiScore = msiScore;
    }
}
