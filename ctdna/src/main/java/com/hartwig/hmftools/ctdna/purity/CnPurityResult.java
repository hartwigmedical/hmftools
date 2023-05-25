package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.ctdna.purity.ResultsWriter.formatPurityValue;

import java.util.StringJoiner;

public class CnPurityResult
{
    public final double FitCoefficient;
    public final double FitIntercept;
    public final double Residuals;
    public final double EstimatedPurity;
    public final boolean Valid;
    public final int CopyNumberSegments;
    public final int GcRatioSegments;
    public final double MedianGcRatiosPerSegments;

    public static final CnPurityResult INVALID_RESULT = new CnPurityResult(
            false, 0, 0, 0, 0,
            0, 0, 0);

    public CnPurityResult(
            final boolean valid, final double fitCoefficient, final double fitIntercept, final double residuals,
            final double estimatedPurity, final int copyNumberSegments, final int gcRatioSegments, final double medianGcRatiosPerSegments)
    {
        Valid = valid;
        FitCoefficient = fitCoefficient;
        FitIntercept = fitIntercept;
        Residuals = residuals;
        EstimatedPurity = estimatedPurity;
        CopyNumberSegments = copyNumberSegments;
        GcRatioSegments = gcRatioSegments;
        MedianGcRatiosPerSegments = medianGcRatiosPerSegments;
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("CopyNumberPurity");
        sj.add("CnFitCoeff");
        sj.add("CnFitIntercept");
        sj.add("CnFitResidualsPerc");
        sj.add("CnSegments");
        sj.add("GcRatioSegments");
        sj.add("GcRatioMedianPerSegment");
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(formatPurityValue(EstimatedPurity));
        sj.add(format("%.4f", FitCoefficient));
        sj.add(format("%.4f", FitIntercept));
        sj.add(format("%.4f", Residuals));
        sj.add(format("%d", CopyNumberSegments));
        sj.add(format("%d", GcRatioSegments));
        sj.add(format("%.1f", MedianGcRatiosPerSegments));

        return sj.toString();
    }
}
