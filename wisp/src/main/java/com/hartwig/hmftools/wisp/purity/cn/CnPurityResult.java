package com.hartwig.hmftools.wisp.purity.cn;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.COPY_NUMBER_LOD_CLONAL_FACTOR;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.COPY_NUMBER_LOD_FACTOR;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatDetectionResult;
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

    // LOD_CNA = 0.4% / sqrt(anueploidy score) / clonalProportion^3

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
        sj.add("CNV_MRD");
        sj.add("CopyNumberPurity");
        sj.add("CopyNumberLod");
        sj.add("CopyNumberPurityLow");
        sj.add("CopyNumberPurityHigh");
        sj.add("AnueploidyScore");
        sj.add("ClonalPercent");
        return sj.toString();
    }

    public double calcLod()
    {
        if(AnueploidyScore < 0 || ClonalPercent <= 0)
            return 1;

        return COPY_NUMBER_LOD_FACTOR / sqrt(AnueploidyScore) / pow(ClonalPercent, COPY_NUMBER_LOD_CLONAL_FACTOR);
    }

    public String toTsv()
    {
        double lod = calcLod();
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(formatDetectionResult(EstimatedPurity, lod));
        sj.add(formatPurityValue(EstimatedPurity));
        sj.add(formatPurityValue(lod));
        sj.add(formatPurityValue(EstimatedPurityLow));
        sj.add(formatPurityValue(EstimatedPurityHigh));
        sj.add(format("%.4f", AnueploidyScore));
        sj.add(format("%.3f", ClonalPercent));

        return sj.toString();
    }
}
