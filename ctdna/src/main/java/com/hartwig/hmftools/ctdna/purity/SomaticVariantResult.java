package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.ctdna.common.CommonUtils.DELIMETER;

import java.util.StringJoiner;

public class SomaticVariantResult
{
    public final int VariantCount;
    public final int AlleleFragmentTotal;
    public final double DepthMedian;
    public final double QualPerAdTotal;
    public final double TumorVaf;
    public final double AdjustedTumorVaf;
    public final double SampleVaf;
    public final double SomaticPurity;

    private boolean mValid;

    public static final SomaticVariantResult INVALID_RESULT = new SomaticVariantResult(false);

    public SomaticVariantResult(
            final boolean valid, final int variantCount, final int alleleFragmentTotal, final double qualPerAdTotal, final double depthMedian,
            final double tumorVaf, final double adjustedTumorVaf, final double sampleVaf, final double somaticPurity)
    {
        VariantCount = variantCount;
        AlleleFragmentTotal = alleleFragmentTotal;
        QualPerAdTotal = qualPerAdTotal;
        DepthMedian = depthMedian;
        TumorVaf = tumorVaf;
        AdjustedTumorVaf = adjustedTumorVaf;
        SampleVaf = sampleVaf;
        SomaticPurity = somaticPurity;
        mValid = valid;
    }

    public SomaticVariantResult(final boolean valid)
    {
        mValid = valid;
        VariantCount = 0;
        AlleleFragmentTotal = 0;
        QualPerAdTotal = 0;
        DepthMedian = 0;
        TumorVaf = 0;
        AdjustedTumorVaf = 0;
        SampleVaf = 0;
        SomaticPurity = 0;
    }

    public boolean valid() { return mValid; }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(DELIMETER);
        sj.add("SomaticPurity");
        sj.add("VariantCount");
        sj.add("DepthMedian");
        sj.add("AlleleFragmentTotal");
        sj.add("QualPerAdTotal");
        sj.add("TumorVaf");
        sj.add("AdjustedTumorVaf");
        sj.add("SampleVaf");
        return sj.toString();
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DELIMETER);
        sj.add(format("%4.3e", SomaticPurity));
        sj.add(format("%d", VariantCount));
        sj.add(format("%.1f", DepthMedian));
        sj.add(format("%d", AlleleFragmentTotal));
        sj.add(format("%.1f", QualPerAdTotal));
        sj.add(format("%4.3e", TumorVaf));
        sj.add(format("%4.3e", AdjustedTumorVaf));
        sj.add(format("%4.3e", SampleVaf));

        return sj.toString();
    }
}
