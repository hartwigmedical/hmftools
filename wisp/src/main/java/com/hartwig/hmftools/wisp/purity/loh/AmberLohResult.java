package com.hartwig.hmftools.wisp.purity.loh;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.wisp.purity.DetectionResult.FALSE;
import static com.hartwig.hmftools.wisp.purity.DetectionResult.NA;
import static com.hartwig.hmftools.wisp.purity.DetectionResult.TRUE;
import static com.hartwig.hmftools.wisp.purity.PurityConstants.LOW_PROBABILITY;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatProbabilityValue;
import static com.hartwig.hmftools.wisp.purity.ResultsWriter.formatPurityValue;

import java.util.StringJoiner;

import com.hartwig.hmftools.wisp.purity.DetectionResult;
import com.hartwig.hmftools.wisp.purity.ResultsWriter;

public class AmberLohResult
{
    public final int RegionCount;
    public final int SiteCount;
    public final double LohPercent;
    public final double AvgCopyNumber;
    public final double MedianAF;
    public final double AvgAF;
    public final double RawEstimatedPurity;
    public final double EstimatedPurity;
    public final double Probability;
    public final int TotalFragments;
    public final double LOD;

    public static final AmberLohResult INVALID_RESULT = new AmberLohResult(
            0, 0, 0, 0, 0, 1,
            0, 0, 1, 1, 0);

    public AmberLohResult(
            final int regionCount, final int siteCount, final double rawEstimatedPurity, double estimatedPurity, final double lohPercent,
            final double lod, final double avgCopyNumber, final double medianAF, final double avgAF, final double pValue, final int totalFragments)
    {
        RegionCount = regionCount;
        SiteCount = siteCount;
        RawEstimatedPurity = rawEstimatedPurity;
        EstimatedPurity = estimatedPurity;
        LohPercent = lohPercent;
        LOD = lod;
        AvgCopyNumber = avgCopyNumber;
        MedianAF = medianAF;
        AvgAF = avgAF;
        Probability = pValue;
        TotalFragments = totalFragments;
    }

    private DetectionResult formatDetectionResult()
    {
        if(ResultsWriter.formatDetectionResult(EstimatedPurity, LOD) == NA)
            return NA;

        if(EstimatedPurity > RawEstimatedPurity || Probability < LOW_PROBABILITY)
            return TRUE;
        else
            return FALSE;
    }

    public static String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add("LOH_MRD");
        sj.add("LOHPurity");
        sj.add("LOHRawPurity");
        sj.add("LOHRegionCount");
        sj.add("LOHSiteCount");
        sj.add("LOHPercent");
        sj.add("LOHLod");
        sj.add("LOHAvgCN");
        sj.add("LOHMedianAF");
        sj.add("LOHMeanAF");
        sj.add("LOHPValue");
        sj.add("LOHFragments");
        return sj.toString();
    }

    public String toTsv()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(formatDetectionResult().toString());
        sj.add(formatPurityValue(EstimatedPurity));
        sj.add(formatPurityValue(RawEstimatedPurity));
        sj.add(String.valueOf(RegionCount));
        sj.add(String.valueOf(SiteCount));
        sj.add(format("%.3f", LohPercent));
        sj.add(formatPurityValue(LOD));
        sj.add(format("%.2f", AvgCopyNumber));
        sj.add(format("%.6f", MedianAF));
        sj.add(format("%.6f", AvgAF));
        sj.add(formatProbabilityValue(Probability));
        sj.add(String.valueOf(TotalFragments));

        return sj.toString();
    }
}
