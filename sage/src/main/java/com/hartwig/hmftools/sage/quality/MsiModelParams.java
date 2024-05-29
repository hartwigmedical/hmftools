package com.hartwig.hmftools.sage.quality;

import static java.lang.Math.exp;
import static java.lang.Math.min;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.sage.SageConstants.MSI_JITTER_MAX_REPEAT_CHANGE;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.basequal.jitter.JitterModelParams;

public class MsiModelParams
{
    private final JitterModelParams mParams;

    private final Map<Integer,Double> mRepeatCountPdfSum;

    public MsiModelParams(final JitterModelParams params)
    {
        mParams = params;
        mRepeatCountPdfSum = Maps.newHashMap();
    }

    public JitterModelParams params() { return mParams; }

    public double calcSkew(int repeatCount, int repeatCountChange)
    {
        /*
         skew = microsatelliteSkew = 1.0164
         scale = scaleFitIntercept + scaleFitGradient * repeat_length (number of units) = -0.1835 + 0.0444 * 8 = 0.1717

        modified_asymmetric_laplace(x, scale, raw_skew):
            skew = 1 + (raw_skew - 1) * min(scale, 1)
            return raw_asymmetric_laplace(x, scale, skew) / pdf_sum

        Where the raw asymmetric laplace distribution is as follows:

        raw_asymmetric_laplace(x, scale, skew)
            return 1 / (scale * (skew+ 1/skew))) * e ^ ((-x * sgn(x) * skew ^ sgn(x)) / scale)
        */

        double scale = mParams.ScaleFitIntercept + mParams.ScaleFitGradient * repeatCount;

        double pdfSum = getOrCalcPdfSum(repeatCount);

        if(pdfSum == 0)
            return 0;

        double modifiedSkew = modifiedSkew(scale);
        double calc = calcAsymmetricLaplace(scale, modifiedSkew, repeatCountChange);

        return calc / pdfSum;
    }

    private double modifiedSkew(double scale)
    {
        return 1 + (mParams.MicrosatelliteSkew - 1) * min(scale, 1);
    }

    private double getOrCalcPdfSum(int repeatCount)
    {
        Double existingPdfSum = mRepeatCountPdfSum.get(repeatCount);

        if(existingPdfSum != null)
            return existingPdfSum;

        // the pdf_sum is a normalisation factor obtained by summing the values of asymmetric_laplace(x, scale, skew) for all x from –5 to +5
        double pdfSum = 0;

        for(int i = -MSI_JITTER_MAX_REPEAT_CHANGE; i <= MSI_JITTER_MAX_REPEAT_CHANGE; ++i)
        {
            double scale = mParams.ScaleFitIntercept + mParams.ScaleFitGradient * repeatCount;
            double modifiedSkew = modifiedSkew(scale);

            double calc = calcAsymmetricLaplace(scale, modifiedSkew, i);
            pdfSum += calc;
        }

        mRepeatCountPdfSum.put(repeatCount, pdfSum);
        return pdfSum;
    }

    private double calcAsymmetricLaplace(double scale, double skew, int repeatCountChange)
    {
        int sign = repeatCountChange > 0 ? 1 : -1;
        return 1 / (scale * (skew + 1 / skew)) * exp((-repeatCountChange * sign * pow(skew, sign)) / scale);
    }
}
