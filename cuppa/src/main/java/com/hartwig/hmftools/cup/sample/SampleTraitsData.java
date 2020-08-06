package com.hartwig.hmftools.cup.sample;

import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;

import java.util.Map;

import com.hartwig.hmftools.common.purple.gender.Gender;

public class SampleTraitsData
{
    public final String SampleId;
    public final Gender GenderType;
    public final boolean HasWGD;
    public final double Purity;
    public final double Ploidy;
    public final int SnvCount;
    public final double IndelsMbPerMb;

    public SampleTraitsData(final String sampleId, final Gender genderType, final boolean hasWGD, final double purity, final double ploidy,
            final int snvCount, final double indelsMbPerMb)
    {
        SampleId = sampleId;
        GenderType = genderType;
        HasWGD = hasWGD;
        Purity = purity;
        Ploidy = ploidy;
        SnvCount = snvCount;
        IndelsMbPerMb = indelsMbPerMb;
    }

    public double getDoubleValue(final SampleTraitType type)
    {
        switch(type)
        {
            case PURITY: return Purity;
            case PLOIDY: return Ploidy;
            case SNV_COUNT: return SnvCount;
            case MS_INDELS_TMB: return IndelsMbPerMb;
            default: return -1;
        }
    }

    public String getStrValue(final SampleTraitType type)
    {
        switch(type)
        {
            case GENDER: return GenderType.toString();
            case WGD: return String.valueOf(HasWGD);
            case PURITY: return String.valueOf(Purity);
            case PLOIDY: return String.valueOf(Ploidy);
            case SNV_COUNT: return String.valueOf(SnvCount);
            case MS_INDELS_TMB: return String.valueOf(IndelsMbPerMb);
            default: return "";
        }
    }

    public static SampleTraitsData from(final Map<String,Integer> fieldsIndexMap, final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);

        return new SampleTraitsData(
                items[fieldsIndexMap.get("SampleId")],
                Gender.valueOf(items[fieldsIndexMap.get("Gender")]),
                items[fieldsIndexMap.get("WholeGenomeDuplication")].equals("1"),
                Double.parseDouble(items[fieldsIndexMap.get("Purity")]),
                Double.parseDouble(items[fieldsIndexMap.get("Ploidy")]),
                Integer.parseInt(items[fieldsIndexMap.get("SnvCount")]),
                Double.parseDouble(items[fieldsIndexMap.get("MsIndelsPerMb")]));
    }
}
