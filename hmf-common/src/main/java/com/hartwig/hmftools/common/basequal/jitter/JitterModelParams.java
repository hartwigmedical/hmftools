package com.hartwig.hmftools.common.basequal.jitter;

public class JitterModelParams
{
    public final String RepeatUnit;
    public double OptimalScaleRepeat4;
    public double OptimalScaleRepeat5;
    public double OptimalScaleRepeat6;
    public double ScaleFitGradient;
    public double ScaleFitIntercept;
    public double MicrosatelliteSkew;

    public JitterModelParams(
            final String repeatUnit,
            final double optimalScaleRepeat4, final double optimalScaleRepeat5, final double optimalScaleRepeat6,
            final double scaleFitGradient, final double scaleFitIntercept, final double microsatelliteSkew)
    {
        RepeatUnit = repeatUnit;
        OptimalScaleRepeat4 = optimalScaleRepeat4;
        OptimalScaleRepeat5 = optimalScaleRepeat5;
        OptimalScaleRepeat6 = optimalScaleRepeat6;
        ScaleFitGradient = scaleFitGradient;
        ScaleFitIntercept = scaleFitIntercept;
        MicrosatelliteSkew = microsatelliteSkew;
    }
}
