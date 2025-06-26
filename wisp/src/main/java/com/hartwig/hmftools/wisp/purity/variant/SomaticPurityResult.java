package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.wisp.purity.DetectionResult.FALSE;
import static com.hartwig.hmftools.wisp.purity.DetectionResult.NA;
import static com.hartwig.hmftools.wisp.purity.DetectionResult.TRUE;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatDetectionResult;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatProbabilityValue;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatPurityValue;
import static com.hartwig.hmftools.wisp.purity.variant.PurityCalcData.CALC_NO_SET;

import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.wisp.purity.DetectionResult;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;

public class SomaticPurityResult
{
    public final int TotalVariants;
    public final String OutlierVariantInfo;
    public final FragmentTotals FragTotals; // used in purity fit - passes filters and either avg qual > threshold or has no allele fragments
    public final UmiTypeCounts UmiCounts;

    public final PurityCalcData PurityCalcs;

    private boolean mValid;

    public static final SomaticPurityResult INVALID_RESULT = new SomaticPurityResult(false);

    public SomaticPurityResult(
            boolean valid, int totalVariants, final String outlierVariantInfo, final FragmentTotals fragmentTotals,
            final UmiTypeCounts umiTypeCounts, final PurityCalcData purityCalcData)
    {
        TotalVariants = totalVariants;
        OutlierVariantInfo = outlierVariantInfo;
        FragTotals = fragmentTotals;
        UmiCounts = umiTypeCounts;
        PurityCalcs = purityCalcData;

        mValid = valid;
    }

    public SomaticPurityResult(final boolean valid)
    {
        mValid = valid;
        TotalVariants = 0;
        OutlierVariantInfo = "";
        PurityCalcs = new PurityCalcData();
        FragTotals = new FragmentTotals();
        UmiCounts = UmiTypeCounts.NO_UMI_COUNTS;
    }

    private DetectionResult formatDetectionResult()
    {
        if(ResultsWriter.formatDetectionResult(PurityCalcs.PurityEstimate, PurityCalcs.LodPurityEstimate) == NA)
            return NA;

        return PurityCalcs.Probability < LOW_PROBABILITY ? TRUE : FALSE;
    }

    public boolean valid() { return mValid; }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("SNV_MRD");
        sj.add("TotalVariants");
        sj.add("ChipVariants");
        sj.add("SNVPurity");
        sj.add("RawSNVPurity");
        sj.add("SNVPValue");
        sj.add("SNVPurityLow");
        sj.add("SNVPurityHigh");
        sj.add("ClonalMethod");
        sj.add("Frag1Variants");
        sj.add("Frag2PlusVariants");
        sj.add("ClonalPeakVariants");
        sj.add("ClonalDropoutRate");
        sj.add("DualSNVPValue");
        sj.add("SNVLod");
        sj.add("DualSNVLod");
        sj.add("TotalFragments");
        sj.add("TotalDual");
        sj.add("AlleleFragments");
        sj.add("AlleleDual");
        sj.add("WeightedAvgDepth");
        sj.add("WeightedAvgVCN");
        sj.add("WeightedAvgCN");
        sj.add("PeakBandwidth");
        sj.add("PeakBandwidthLow");
        sj.add("PeakBandwidthHigh");
        sj.add("OutlierVariants");
        sj.add("ErrorRate");
        sj.add("RawBqrErrorRate");
        sj.add("BqrThreshold");
        sj.add("BqrExtraInfo");
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(formatDetectionResult().toString());
        sj.add(format("%d", TotalVariants));
        sj.add(format("%d", FragTotals.variantCount()));
        sj.add(formatPurityValue(PurityCalcs.PurityEstimate));
        sj.add(formatPurityValue(PurityCalcs.RawPurityEstimate));
        sj.add(formatProbabilityValue(PurityCalcs.Probability));
        sj.add(formatPurityValue(PurityCalcs.PurityRangeLow));
        sj.add(formatPurityValue(PurityCalcs.PurityRangeHigh));
        sj.add(String.valueOf(PurityCalcs.Clonality.Method));

        sj.add(format("%d", FragTotals.sampleOneFragmentCount()));
        sj.add(format("%d", FragTotals.sampleTwoPlusCount()));
        sj.add(format("%d", PurityCalcs.Clonality.VariantCount));
        sj.add(format("%.2f", PurityCalcs.Clonality.DropoutRate));
        sj.add(formatProbabilityValue(PurityCalcs.DualProbability));
        sj.add(formatPurityValue(PurityCalcs.LodPurityEstimate));
        sj.add(formatPurityValue(PurityCalcs.DualLodPurityEstimate));

        sj.add(format("%d", FragTotals.sampleDepthTotal()));
        sj.add(format("%d", UmiCounts.TotalDual));
        sj.add(format("%d", FragTotals.sampleAdTotal()));
        sj.add(format("%d", UmiCounts.AlleleDual));
        sj.add(format("%.1f", FragTotals.weightedSampleDepth()));
        sj.add(format("%.1f", FragTotals.weightedVariantCopyNumber()));
        sj.add(format("%.1f", FragTotals.weightedCopyNumber()));
        sj.add(format("%.4f", PurityCalcs.Clonality.PeakBandwidth));
        sj.add(format("%.4f", PurityCalcs.Clonality.PeakBandwidthLow));
        sj.add(format("%.4f", PurityCalcs.Clonality.PeakBandwidthHigh));


        sj.add(OutlierVariantInfo);

        sj.add(format("%.6f", PurityCalcs.ErrorRate));
        sj.add(format("%.6f", PurityCalcs.RawBqrErrorRate));
        sj.add(format("%d", PurityCalcs.BqrQualThreshold));

        String bqrDataStr = PurityCalcs.BqrExtraInfo.stream().collect(Collectors.joining(";"));
        sj.add(bqrDataStr);

        return sj.toString();
    }
}
