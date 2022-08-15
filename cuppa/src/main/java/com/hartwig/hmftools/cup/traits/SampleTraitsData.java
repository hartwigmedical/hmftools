package com.hartwig.hmftools.cup.traits;

import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaConfig.FLD_SAMPLE_ID;

import java.util.Map;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurityContext;

public class SampleTraitsData
{
    public final String SampleId;
    public final Gender GenderType;
    public final boolean HasWGD;
    public final double Purity;
    public final double Ploidy;
    public final double IndelsMbPerMb;
    public final double ChordHrd;

    public SampleTraitsData(final String sampleId, final Gender genderType, final boolean hasWGD, final double purity, final double ploidy,
            final double indelsMbPerMb, final double chordHrd)
    {
        SampleId = sampleId;
        GenderType = genderType;
        HasWGD = hasWGD;
        Purity = purity;
        Ploidy = ploidy;
        IndelsMbPerMb = indelsMbPerMb;
        ChordHrd = chordHrd;
    }

    public double getDoubleValue(final SampleTraitType type)
    {
        switch(type)
        {
            case PURITY: return Purity;
            case PLOIDY: return Ploidy;
            case MS_INDELS_TMB: return IndelsMbPerMb;
            case CHORD_HRD: return ChordHrd;
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
            case MS_INDELS_TMB: return String.valueOf(IndelsMbPerMb);
            case CHORD_HRD: return String.valueOf(ChordHrd);
            default: return "";
        }
    }

    public static SampleTraitsData from(final String sampleId, final PurityContext purityContext, double chordHrd)
    {
        return new SampleTraitsData(
                sampleId, purityContext.gender(), purityContext.wholeGenomeDuplication(),
                purityContext.bestFit().purity(), purityContext.bestFit().ploidy(),
                purityContext.microsatelliteIndelsPerMb(), chordHrd);
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(FLD_SAMPLE_ID);
        sj.add("Gender");
        sj.add("WholeGenomeDuplication");
        sj.add("Purity");
        sj.add("Ploidy");
        sj.add("MsIndelsPerMb");
        sj.add("ChordHrd");
        return sj.toString();
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(GenderType.toString());
        sj.add(String.valueOf(HasWGD));
        sj.add(String.valueOf(Purity));
        sj.add(String.valueOf(Ploidy));
        sj.add(String.valueOf(IndelsMbPerMb));
        sj.add(String.valueOf(ChordHrd));
        return sj.toString();
    }

    public static SampleTraitsData from(final Map<String,Integer> fieldsIndexMap, final String data)
    {
        final String[] items = data.split(DATA_DELIM, -1);

        String wgdValue = items[fieldsIndexMap.get("WholeGenomeDuplication")];

        return new SampleTraitsData(
                items[fieldsIndexMap.get("SampleId")],
                Gender.valueOf(items[fieldsIndexMap.get("Gender")]),
                wgdValue.equals("1") || wgdValue.equalsIgnoreCase("true"),
                Double.parseDouble(items[fieldsIndexMap.get("Purity")]),
                Double.parseDouble(items[fieldsIndexMap.get("Ploidy")]),
                Double.parseDouble(items[fieldsIndexMap.get("MsIndelsPerMb")]),
                Double.parseDouble(items[fieldsIndexMap.get("ChordHrd")]));
    }
}
