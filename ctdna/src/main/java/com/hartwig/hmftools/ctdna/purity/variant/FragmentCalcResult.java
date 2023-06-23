package com.hartwig.hmftools.ctdna.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.formatProbabilityValue;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.formatPurityValue;

import java.util.StringJoiner;

public class FragmentCalcResult
{
    public final double VAF;
    public final double EstimatedPurity;
    public final double PurityProbability;
    public final double PurityRangeLow;
    public final double PurityRangeHigh;

    public static final FragmentCalcResult INVALID = new FragmentCalcResult(
            0, 0, 0, 0, 0);

    public FragmentCalcResult(
            final double vaf, final double estimatedPurity, final double purityProbability,
            final double purityRangeLow, final double purityRangeHigh)
    {
        VAF = vaf;
        EstimatedPurity = estimatedPurity;
        PurityProbability = purityProbability;
        PurityRangeLow = purityRangeLow;
        PurityRangeHigh = purityRangeHigh;
    }

    public static String header(final String context)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(format("%sSomaticPurity", context));
        sj.add(format("%sSomaticProbability", context));
        sj.add(format("%sSampleVAF", context));
        sj.add(format("%sSomaticPurityLow", context));
        sj.add(format("%sSomaticPurityHigh", context));
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(formatPurityValue(EstimatedPurity));
        sj.add(formatProbabilityValue(PurityProbability));
        sj.add(formatPurityValue(VAF));
        sj.add(formatPurityValue(PurityRangeLow));
        sj.add(formatPurityValue(PurityRangeHigh));
        return sj.toString();
    }

}
