package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.formatPurityValue;

import java.util.StringJoiner;

public class SomaticVariantResult
{
    public final int TotalVariants;
    public final int CalcVariants;
    public final int AlleleFragmentTotal;
    public final double DepthMedian;
    public final int NonZeroDepthCount;
    public final double NonZeroDepthMedian;
    public final double QualPerAdTotal;
    public final double TumorVaf;
    public final double AdjustedTumorVaf;
    public final double SampleVaf;
    public final double SomaticPurity;
    public final double Probability;

    private boolean mValid;

    public static final SomaticVariantResult INVALID_RESULT = new SomaticVariantResult(false);

    public SomaticVariantResult(
            boolean valid, int totalVariants, int calcVariants, int alleleFragmentTotal, double qualPerAdTotal,
            double depthMedian, int nonZeroDepthCount, double nonZeroDepthMedian,
            double tumorVaf, double adjustedTumorVaf, double sampleVaf, double somaticPurity, double probability)
    {
        TotalVariants = totalVariants;
        CalcVariants = calcVariants;
        AlleleFragmentTotal = alleleFragmentTotal;
        QualPerAdTotal = qualPerAdTotal;
        DepthMedian = depthMedian;
        NonZeroDepthCount = nonZeroDepthCount;
        NonZeroDepthMedian = nonZeroDepthMedian;
        TumorVaf = tumorVaf;
        AdjustedTumorVaf = adjustedTumorVaf;
        SampleVaf = sampleVaf;
        SomaticPurity = somaticPurity;
        Probability = probability;
        mValid = valid;
    }

    public SomaticVariantResult(final boolean valid)
    {
        mValid = valid;
        TotalVariants = 0;
        CalcVariants = 0;
        AlleleFragmentTotal = 0;
        QualPerAdTotal = 0;
        DepthMedian = 0;
        NonZeroDepthCount = 0;
        NonZeroDepthMedian = 0;
        TumorVaf = 0;
        AdjustedTumorVaf = 0;
        SampleVaf = 0;
        SomaticPurity = 0;
        Probability = 0;
    }

    public boolean valid() { return mValid; }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("SomaticPurity");
        sj.add("SomaticProbability");
        sj.add("TotalVariants");
        sj.add("CalcVariants");
        sj.add("DepthMedian");
        sj.add("NonZeroDepthCount");
        sj.add("NonZeroDepthMedian");
        sj.add("AlleleFragmentTotal");
        sj.add("QualPerAdTotal");
        sj.add("TumorVaf");
        sj.add("AdjustedTumorVaf");
        sj.add("SampleVaf");
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(formatPurityValue(SomaticPurity));
        sj.add(formatPurityValue(Probability));
        sj.add(format("%d", TotalVariants));
        sj.add(format("%d", CalcVariants));
        sj.add(format("%.1f", DepthMedian));
        sj.add(format("%d", NonZeroDepthCount));
        sj.add(format("%.1f", NonZeroDepthMedian));
        sj.add(format("%d", AlleleFragmentTotal));
        sj.add(format("%.1f", QualPerAdTotal));
        sj.add(formatPurityValue(TumorVaf));
        sj.add(formatPurityValue(AdjustedTumorVaf));
        sj.add(formatPurityValue(SampleVaf));

        return sj.toString();
    }
}
