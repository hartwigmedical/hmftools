package com.hartwig.hmftools.wisp.purity.cn;

public class CnFitResult
{
    public final double FitCoefficient;
    public final double FitIntercept;
    public final double EstimatedPurity;
    public final double Residuals;

    public CnFitResult(final double fitCoefficient, final double fitIntercept, final double estimatedPurity, final double residuals)
    {
        FitCoefficient = fitCoefficient;
        FitIntercept = fitIntercept;
        EstimatedPurity = estimatedPurity;
        Residuals = residuals;
    }
}
