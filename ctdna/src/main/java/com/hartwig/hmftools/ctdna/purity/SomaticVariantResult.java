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
    public final int UmiRefNonDual;
    public final int UmiRefDual;
    public final int AlleleFragments; //
    public final int UmiAlleleNonDual;
    public final int UmiAlleleDual;
    public final double DepthMedian;
    public final int NonZeroDepthCount;
    public final double NonZeroDepthMedian;

    public final double QualPerAdTotal;

    public final double TumorVaf;
    public final double AdjustedTumorVaf;
    public final FragmentCalcResult AllFragsResult;
    public final FragmentCalcResult DualFragsResult;
    public final FragmentCalcResult LimitOfDetectionResult;

    private boolean mValid;

    public static final SomaticVariantResult INVALID_RESULT = new SomaticVariantResult(false);

    public SomaticVariantResult(
            boolean valid, int totalVariants, int calcVariants,
            final SomaticVariantCounts sampleCounts, final UmiTypeCounts umiTypeCounts,
            double qualPerAdTotal, double tumorVaf, double adjustedTumorVaf,
            final FragmentCalcResult allFragsResult, final FragmentCalcResult dualFragsResult, final FragmentCalcResult lodFragsResult)
    {
        TotalVariants = totalVariants;
        CalcVariants = calcVariants;
        TotalFragments = sampleCounts.totalFragments();
        UmiRefNonDual = umiTypeCounts.RefNone + umiTypeCounts.RefSingle;
        UmiRefDual = umiTypeCounts.RefDual;
        AlleleFragments = sampleCounts.alleleFragments();
        UmiAlleleNonDual = umiTypeCounts.AlleleNone + umiTypeCounts.AlleleSingle;
        UmiAlleleDual = umiTypeCounts.AlleleDual;
        QualPerAdTotal = qualPerAdTotal;
        DepthMedian = sampleCounts.medianDepth(false);
        NonZeroDepthCount = sampleCounts.nonZeroDepth();
        NonZeroDepthMedian = sampleCounts.medianDepth(true);
        TumorVaf = tumorVaf;
        AdjustedTumorVaf = adjustedTumorVaf;
        AllFragsResult = allFragsResult;
        DualFragsResult = dualFragsResult;
        LimitOfDetectionResult = lodFragsResult;
        mValid = valid;
    }

    public SomaticVariantResult(final boolean valid)
    {
        mValid = valid;
        TotalVariants = 0;
        CalcVariants = 0;
        TotalFragments = 0;
        UmiRefNonDual = 0;
        UmiRefDual = 0;
        UmiAlleleNonDual = 0;
        UmiAlleleDual = 0;
        AlleleFragments = 0;
        QualPerAdTotal = 0;
        DepthMedian = 0;
        NonZeroDepthCount = 0;
        NonZeroDepthMedian = 0;
        TumorVaf = 0;
        AdjustedTumorVaf = 0;
        AllFragsResult = FragmentCalcResult.INVALID;
        DualFragsResult = FragmentCalcResult.INVALID;
        LimitOfDetectionResult = FragmentCalcResult.INVALID;
    }

    public boolean valid() { return mValid; }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(FragmentCalcResult.header("")); // was All
        sj.add(FragmentCalcResult.header("Dual"));
        sj.add("LodPurity");
        sj.add("TotalVariants");
        sj.add("CalcVariants");
        sj.add("TotalFragments");
        sj.add("UmiRefNonDual");
        sj.add("UmiRefDual");
        sj.add("AlleleFragments");
        sj.add("UmiAlleleNonDual");
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

        // results
        sj.add(AllFragsResult.toTsv());
        sj.add(DualFragsResult.toTsv());
        sj.add(formatPurityValue(LimitOfDetectionResult.EstimatedPurity));

        // inputs
        sj.add(format("%d", TotalVariants));
        sj.add(format("%d", CalcVariants));
        sj.add(format("%d", TotalFragments));
        sj.add(format("%d", UmiRefNonDual));
        sj.add(format("%d", UmiRefDual));
        sj.add(format("%d", AlleleFragments));
        sj.add(format("%d", UmiAlleleNonDual));
        sj.add(format("%d", UmiAlleleDual));
        sj.add(format("%.1f", QualPerAdTotal));
        sj.add(formatPurityValue(TumorVaf));
        sj.add(formatPurityValue(AdjustedTumorVaf));
        sj.add(format("%.1f", DepthMedian));
        sj.add(format("%d", NonZeroDepthCount));
        sj.add(format("%.1f", NonZeroDepthMedian));

        return sj.toString();
    }
}
