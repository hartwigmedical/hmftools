package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.formatPurityValue;

import java.util.StringJoiner;

public class SomaticVariantResult
{
    public final int TotalVariants;
    public final int CalcVariants; // used in purity fit - passes filters and either avg qual > threshold or has no allele fragments

    public final int TotalFragments;
    public final int UmiRefNone;
    public final int UmiRefSingle;
    public final int UmiRefDual;
    public final int AlleleFragments; //
    public final int UmiAlleleNone;
    public final int UmiAlleleSingle;
    public final int UmiAlleleDual;
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
            boolean valid, int totalVariants, int calcVariants,
            int totalFragments, int umiRefNone, int umiRefSingle, int umiRefDual,
            int alleleFragments, int umiAlleleNone, int umiAlleleSingle, int umiAlleleDual,
            double qualPerAdTotal, double depthMedian, int nonZeroDepthCount, double nonZeroDepthMedian,
            double tumorVaf, double adjustedTumorVaf, double sampleVaf, double somaticPurity, double probability)
    {
        TotalVariants = totalVariants;
        CalcVariants = calcVariants;
        TotalFragments = totalFragments;
        UmiRefNone = umiRefNone;
        UmiRefSingle = umiRefSingle;
        UmiRefDual = umiRefDual;
        AlleleFragments = alleleFragments;
        UmiAlleleNone = umiAlleleNone;
        UmiAlleleSingle = umiAlleleSingle;
        UmiAlleleDual = umiAlleleDual;
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
        TotalFragments = 0;
        UmiRefNone = 0;
        UmiRefSingle = 0;
        UmiRefDual = 0;
        UmiAlleleNone = 0;
        UmiAlleleSingle = 0;
        UmiAlleleDual = 0;
        AlleleFragments = 0;
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
        sj.add("TotalFragments");
        sj.add("UmiRefNone");
        sj.add("UmiRefSingle");
        sj.add("UmiRefDual");
        sj.add("AlleleFragments");
        sj.add("UmiAlleleNone");
        sj.add("UmiAlleleSingle");
        sj.add("UmiAlleleDual");
        sj.add("QualPerAdTotal");
        sj.add("TumorVaf");
        sj.add("AdjustedTumorVaf");
        sj.add("SampleVaf");
        sj.add("DepthMedian");
        sj.add("NonZeroDepthCount");
        sj.add("NonZeroDepthMedian");
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(formatPurityValue(SomaticPurity));
        sj.add(formatPurityValue(Probability));
        sj.add(format("%d", TotalVariants));
        sj.add(format("%d", CalcVariants));
        sj.add(format("%d", TotalFragments));
        sj.add(format("%d", UmiRefNone));
        sj.add(format("%d", UmiRefSingle));
        sj.add(format("%d", UmiRefDual));
        sj.add(format("%d", AlleleFragments));
        sj.add(format("%d", UmiAlleleNone));
        sj.add(format("%d", UmiAlleleSingle));
        sj.add(format("%d", UmiAlleleDual));
        sj.add(format("%.1f", QualPerAdTotal));
        sj.add(formatPurityValue(TumorVaf));
        sj.add(formatPurityValue(AdjustedTumorVaf));
        sj.add(formatPurityValue(SampleVaf));
        sj.add(format("%.1f", DepthMedian));
        sj.add(format("%d", NonZeroDepthCount));
        sj.add(format("%.1f", NonZeroDepthMedian));

        return sj.toString();
    }
}
