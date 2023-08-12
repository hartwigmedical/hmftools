package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.util.StringJoiner;

import com.hartwig.hmftools.wisp.purity.ResultsWriter;

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
        sj.add(format("%sSNVPurity", context));
        sj.add(format("%sSNVProbability", context));
        sj.add(format("%sSampleVAF", context));
        sj.add(format("%sSNVPurityLow", context));
        sj.add(format("%sSNVPurityHigh", context));
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(ResultsWriter.formatPurityValue(EstimatedPurity));
        sj.add(ResultsWriter.formatProbabilityValue(PurityProbability));
        sj.add(ResultsWriter.formatPurityValue(VAF));
        sj.add(ResultsWriter.formatPurityValue(PurityRangeLow));
        sj.add(ResultsWriter.formatPurityValue(PurityRangeHigh));
        return sj.toString();
    }

}
