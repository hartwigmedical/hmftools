package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatProbabilityValue;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatPurityValue;

import java.util.StringJoiner;

public class SomaticVariantResult
{
    public final int TotalVariants;
    public final int CalcVariants; // used in purity fit - passes filters and either avg qual > threshold or has no allele fragments
    public final int Frag1Variants;
    public final int Frag2PlusVariants;
    public final ClonalityMethod Method;
    public final int ClonalVariants;
    public final double ClonalDropoutRate;
    public final double WeightedAvgDepth;
    public final int TotalFragments;
    public final int TotalDualFragments;
    public final int AlleleFragments;
    public final int AlleleDualFragments;

    public final double TumorVaf;
    public final double AdjustedTumorVaf;
    public final double RawSomaticPurity;
    public final FragmentCalcResult AllFragsResult;
    public final FragmentCalcResult DualFragsResult;
    public final FragmentCalcResult LimitOfDetectionResult;

    private boolean mValid;

    public static final SomaticVariantResult INVALID_RESULT = new SomaticVariantResult(false);

    public SomaticVariantResult(
            boolean valid, int totalVariants, int calcVariants, int frag1Variants, int frag2PlusVariants, int clonalVariants,
            final ClonalityMethod method, final double dropoutRate, double weightedAvgDepth,
            final SomaticVariantCounts sampleCounts, final UmiTypeCounts umiTypeCounts,
            double tumorVaf, double adjustedTumorVaf, double rawSomaticPurity,
            final FragmentCalcResult allFragsResult, final FragmentCalcResult dualFragsResult, final FragmentCalcResult lodFragsResult)
    {
        TotalVariants = totalVariants;
        CalcVariants = calcVariants;
        Frag1Variants = frag1Variants;
        Frag2PlusVariants = frag2PlusVariants;
        ClonalVariants = clonalVariants;
        ClonalDropoutRate = dropoutRate;
        Method = method;
        WeightedAvgDepth = weightedAvgDepth;
        TotalFragments = sampleCounts.totalFragments();
        TotalDualFragments = umiTypeCounts.TotalDual;
        AlleleFragments = sampleCounts.alleleFragments();
        AlleleDualFragments = umiTypeCounts.AlleleDual;
        TumorVaf = tumorVaf;
        AdjustedTumorVaf = adjustedTumorVaf;
        RawSomaticPurity = rawSomaticPurity;
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
        Frag1Variants = 0;
        Frag2PlusVariants = 0;
        ClonalVariants = 0;
        Method = ClonalityMethod.NONE;
        ClonalDropoutRate = 0;
        WeightedAvgDepth = 0;
        TotalFragments = 0;
        TotalDualFragments = 0;
        AlleleDualFragments = 0;
        AlleleFragments = 0;
        TumorVaf = 0;
        AdjustedTumorVaf = 0;
        RawSomaticPurity = 0;
        AllFragsResult = FragmentCalcResult.INVALID;
        DualFragsResult = FragmentCalcResult.INVALID;
        LimitOfDetectionResult = FragmentCalcResult.INVALID;
    }

    public boolean valid() { return mValid; }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("TotalVariants");
        sj.add("CalcVariants");
        sj.add("Frag1Variants");
        sj.add("Frag2PlusVariants");
        sj.add("ClonalPeakVariants");
        sj.add("ClonalMethod");
        sj.add("ClonalDropoutRate");
        sj.add("RawSomaticPurity");
        sj.add(FragmentCalcResult.header(""));
        sj.add("DualSNVProbability");
        sj.add("LodPurity");
        sj.add("TotalFragments");
        sj.add("TotalDual");
        sj.add("AlleleFragments");
        sj.add("AlleleDual");
        sj.add("TumorVaf");
        sj.add("AdjustedTumorVaf");
        sj.add("WeightedAvgDepth");
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(format("%d", TotalVariants));
        sj.add(format("%d", CalcVariants));
        sj.add(format("%d", Frag1Variants));
        sj.add(format("%d", Frag2PlusVariants));
        sj.add(format("%d", ClonalVariants));
        sj.add(String.valueOf(Method));
        sj.add(format("%.2f", ClonalDropoutRate));
        sj.add(formatPurityValue(RawSomaticPurity));
        sj.add(AllFragsResult.toTsv());
        sj.add(formatProbabilityValue(DualFragsResult.PurityProbability));
        sj.add(formatPurityValue(LimitOfDetectionResult.EstimatedPurity));

        // inputs
        sj.add(format("%d", TotalFragments));
        sj.add(format("%d", TotalDualFragments));
        sj.add(format("%d", AlleleFragments));
        sj.add(format("%d", AlleleDualFragments));
        sj.add(formatPurityValue(TumorVaf));
        sj.add(formatPurityValue(AdjustedTumorVaf));
        sj.add(format("%.0f", WeightedAvgDepth));

        return sj.toString();
    }
}
