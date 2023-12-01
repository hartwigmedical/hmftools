package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatProbabilityValue;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatPurityValue;

import java.util.StringJoiner;

public class SomaticPurityResult
{
    public final int TotalVariants;
    public final FragmentTotals FragTotals; // used in purity fit - passes filters and either avg qual > threshold or has no allele fragments
    public final UmiTypeCounts UmiCounts;

    /*
    public final int TotalFragments;
    public final int TotalDualFragments;
    public final int AlleleFragments;
    public final int AlleleDualFragments;

    public final double TumorVaf;
    public final double AdjustedTumorVaf;
    */
    public final PurityCalcData PurityCalcs;

    /*
    public final FragmentCalcResult DualFragsResult;
    public final FragmentCalcResult LimitOfDetectionResult;
    */

    private boolean mValid;

    public static final SomaticPurityResult INVALID_RESULT = new SomaticPurityResult(false);

    public SomaticPurityResult(
            boolean valid, int totalVariants, final FragmentTotals fragmentTotals,
            final UmiTypeCounts umiTypeCounts, final PurityCalcData purityCalcData)
    {
        TotalVariants = totalVariants;
        FragTotals = fragmentTotals;
        UmiCounts = umiTypeCounts;
        PurityCalcs = purityCalcData;

        /*
        TotalFragments = sampleCounts.totalFragments();
        AlleleFragments = sampleCounts.alleleFragments();
        TotalDualFragments = umiTypeCounts.TotalDual;
        AlleleDualFragments = umiTypeCounts.AlleleDual;
        TumorVaf = tumorVaf;
        AdjustedTumorVaf = adjustedTumorVaf;
        FragCalcResult = fragCalcResult;
        DualFragsResult = dualFragsResult;
        LimitOfDetectionResult = lodFragsResult;
        */
        mValid = valid;
    }

    public SomaticPurityResult(final boolean valid)
    {
        mValid = valid;
        TotalVariants = 0;
        PurityCalcs = new PurityCalcData();
        FragTotals = new FragmentTotals();
        UmiCounts = UmiTypeCounts.NO_UMI_COUNTS;
        // DualFragsResult = FragmentCalcResult.INVALID;
        // LimitOfDetectionResult = FragmentCalcResult.INVALID;
    }

    public boolean valid() { return mValid; }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("TotalVariants");
        sj.add("CalcVariants");
        sj.add("SNVPurity");
        sj.add("RawSomaticPurity");
        sj.add("ClonalMethod");
        sj.add("TumorVaf");
        sj.add("AdjSampleVaf");
        sj.add("Frag1Variants");
        sj.add("Frag2PlusVariants");
        sj.add("ClonalPeakVariants");
        sj.add("ClonalDropoutRate");
        sj.add("DualSNVProbability");
        sj.add("LodPurity");
        sj.add("TotalFragments");
        sj.add("TotalDual");
        sj.add("AlleleFragments");
        sj.add("AlleleDual");
        sj.add("WeightedAvgDepth");
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(format("%d", TotalVariants));
        sj.add(format("%d", FragTotals.variantCount()));
        sj.add(formatPurityValue(PurityCalcs.PurityEstimate));
        sj.add(formatPurityValue(PurityCalcs.RawPurityEstimate));
        sj.add(String.valueOf(PurityCalcs.Clonality.Method));
        sj.add(formatPurityValue(FragTotals.rawTumorVaf()));
        sj.add(formatPurityValue(FragTotals.adjSampleVaf()));

        sj.add(format("%d", FragTotals.sampleOneFragmentCount()));
        sj.add(format("%d", FragTotals.sampleTwoPlusCount()));
        sj.add(format("%d", PurityCalcs.Clonality.VariantCount));
        sj.add(format("%.2f", PurityCalcs.Clonality.DropoutRate));
        sj.add(formatProbabilityValue(PurityCalcs.DualProbability));
        sj.add(formatPurityValue(PurityCalcs.LodPurityEstimate));

        sj.add(format("%d", FragTotals.sampleDepthTotal()));
        sj.add(format("%d", UmiCounts.TotalDual));
        sj.add(format("%d", FragTotals.sampleAdTotal()));
        sj.add(format("%d", UmiCounts.AlleleDual));
        sj.add(format("%.1f", FragTotals.weightedSampleDepth()));

        return sj.toString();
    }
}
