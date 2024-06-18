package com.hartwig.hmftools.wisp.purity.cn;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatPurityValue;

import java.util.StringJoiner;

public class CnPurityResult
{
    public final double Residuals;
    public final double EstimatedPurity;
    public final double EstimatedPurityLow;
    public final double EstimatedPurityHigh;
    public final boolean Valid;
    public final int CopyNumberSegments;
    public final int GcRatioSegments;
    public final double MedianGcRatiosPerSegments;
    public final double AnueploidyScore;
    public final double ClonalPercent;

    public static final CnPurityResult INVALID_RESULT = new CnPurityResult(
            false, 0, 0, 0, 0,
            0, 0, 0, 0, 0);

    public CnPurityResult(
            boolean valid, double residuals,
            double estimatedPurity, double estimatedPurityLow, double estimatedPurityHigh, double anueploidyScore, double clonalPercent,
            int copyNumberSegments, int gcRatioSegments, double medianGcRatiosPerSegments)
    {
        Valid = valid;
        Residuals = residuals;
        EstimatedPurity = estimatedPurity;
        EstimatedPurityLow = estimatedPurityLow;
        EstimatedPurityHigh = estimatedPurityHigh;
        CopyNumberSegments = copyNumberSegments;
        GcRatioSegments = gcRatioSegments;
        MedianGcRatiosPerSegments = medianGcRatiosPerSegments;
        AnueploidyScore = anueploidyScore;
        ClonalPercent = clonalPercent;
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("CopyNumberPurity");
        sj.add("CopyNumberPurityLow");
        sj.add("CopyNumberPurityHigh");
        sj.add("AnueploidyScore");
        sj.add("ClonalPercent");
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(formatPurityValue(EstimatedPurity));
        sj.add(formatPurityValue(EstimatedPurityLow));
        sj.add(formatPurityValue(EstimatedPurityHigh));
        sj.add(format("%.4f", AnueploidyScore));
        sj.add(format("%.3f", ClonalPercent));

        return sj.toString();
    }
}
