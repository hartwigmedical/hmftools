package com.hartwig.hmftools.wisp.purity.cn;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatProbabilityValue;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatPurityValue;

import java.util.StringJoiner;

public class AmberLohResult
{
    public final int SiteCount;
    public final double LohPercent;
    public final double AvgCopyNumber;
    public final double AvgAF;
    public final double EstimatedPurity;
    public final double PValue;
    public final int TotalFragments;

    public static final AmberLohResult INVALID_RESULT = new AmberLohResult(
            0, 0, 0, 0, 0, 0, 0);

    public AmberLohResult(
            final int siteCount, final double estimatedPurity, final double lohPercent,
            final double avgCopyNumber, final double avgAF, final double pValue, final int totalFragments)
    {
        SiteCount = siteCount;
        EstimatedPurity = estimatedPurity;
        LohPercent = lohPercent;
        AvgCopyNumber = avgCopyNumber;
        AvgAF = avgAF;
        PValue = pValue;
        TotalFragments = totalFragments;
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("LohEstimatedPurity");
        sj.add("LohSiteCount");
        sj.add("LohPercent");
        sj.add("LohMeanCopyNumber");
        sj.add("LohMeanAF");
        sj.add("LohPValue");
        sj.add("LohFragments");
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(formatPurityValue(EstimatedPurity));
        sj.add(String.valueOf(SiteCount));
        sj.add(format("%.3f", LohPercent));
        sj.add(format("%.2f", AvgCopyNumber));
        sj.add(format("%.6f", AvgAF));
        sj.add(formatProbabilityValue(PValue));
        sj.add(String.valueOf(TotalFragments));

        return sj.toString();
    }
}
