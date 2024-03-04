package com.hartwig.hmftools.wisp.purity.cn;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates.GENOME_LENGTH_V38;
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
        sj.add("LohEstimatedPurity");
        sj.add("LohSiteCount");
        sj.add("LohPercent");
        sj.add("LohAvgCopyNumber");
        sj.add("LohAvgAF");
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(formatPurityValue(EstimtedPurity));
        sj.add(String.valueOf(SiteCount));
        sj.add(format("%.3f", LohDistance / GENOME_LENGTH_V38));
        sj.add(format("%.2f", AvgCopyNumber));
        sj.add(format("%.4f", AvgAF));

        return sj.toString();
    }
}
