package com.hartwig.hmftools.common.basequal.jitter;

public class JitterModelCalc
{
    public static double modifiedAsymmetricLaplace(double x, double scale, double rawSkew)
    {
        double skew = 1 + (rawSkew - 1) * Math.min(scale, 1);
        double pdfSum = 0.0;
        for(double x1 = -5.0; x1 <= 5.0; x1 += 1.0)
        {
            pdfSum += rawAsymmetricLaplace(x1, scale, skew);
        }
        return rawAsymmetricLaplace(x, scale, skew) / pdfSum;
    }

    public static double rawAsymmetricLaplace(double x, double scale, double skew)
    {
        return 1 / (scale * (skew + 1 / skew)) * Math.exp((-x * sgn(x) * Math.pow(skew, sgn(x))) / scale);
    }

    private static double sgn(double x)
    {
        if (x >= 0)
        {
            return 1.0;
        } else return -1.0;
    }
}
