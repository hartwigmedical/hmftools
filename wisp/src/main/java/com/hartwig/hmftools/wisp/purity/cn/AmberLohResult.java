package com.hartwig.hmftools.wisp.purity.cn;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatPurityValue;

import java.util.StringJoiner;

public class AmberLohResult
{
    public final int SiteCount;
    public final int LohDistance;
    public final double AvgCopyNumber;
    public final double AvgAF;
    public final double EstimtedPurity;

    public static final AmberLohResult INVALID_RESULT = new AmberLohResult(
            0, 0, 0, 0, 0);

    public AmberLohResult(
            final int siteCount, final int lohDistance, final double avgCopyNumber, final double avgAF, final double estimtedPurity)
    {
        SiteCount = siteCount;
        LohDistance = lohDistance;
        AvgCopyNumber = avgCopyNumber;
        AvgAF = avgAF;
        EstimtedPurity = estimtedPurity;
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("LohSiteCount");
        sj.add("LohMb");
        sj.add("LohAvgCopyNumber");
        sj.add("LohAvgAF");
        sj.add("LohEstimatedPurity");
        return sj.toString();
    }

    private static final double LOH_MB = 1_000_000;

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(String.valueOf(SiteCount));
        sj.add(format("%.1f", LohDistance / LOH_MB));
        sj.add(format("%.2f", AvgCopyNumber));
        sj.add(format("%.4f", AvgAF));
        sj.add(formatPurityValue(EstimtedPurity));

        return sj.toString();
    }
}
