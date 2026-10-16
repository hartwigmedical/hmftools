package com.hartwig.hmftools.amber;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConstants.LOH_MOVING_ITEM_COUNT;
import static com.hartwig.hmftools.amber.AmberConstants.LOH_REF_ALT_MIN_COUNT;
import static com.hartwig.hmftools.amber.AmberConstants.LOH_REF_ALT_THRESHOLD_HIGH;
import static com.hartwig.hmftools.amber.AmberConstants.LOH_REF_ALT_THRESHOLD_LOW;
import static com.hartwig.hmftools.amber.AmberConstants.LOH_REF_ALT_THRESHOLD_MED;
import static com.hartwig.hmftools.amber.AmberConstants.VALID_BAF_RATE_HIGH;
import static com.hartwig.hmftools.amber.AmberConstants.VALID_BAF_RATE_LOW;
import static com.hartwig.hmftools.amber.AmberConstants.VALID_BAF_RATE_MED;

public class BafLohData
{
    public boolean ValidBafHigh;
    public boolean ValidBafMedium;
    public boolean ValidBafLow;

    public double ValidRateHigh;
    public double ValidRateMedium;
    public double ValidRateLow;

    public boolean LohCandidate;
    public boolean LohRecovered;

    public BafLohData()
    {
        ValidBafHigh = false;
        ValidBafMedium = false;
        ValidBafLow = false;
        ValidRateHigh = 0;
        ValidRateMedium = 0;
        ValidRateLow = 0;
        LohCandidate = false;
        LohRecovered = false;
    }

    public void setBafStatus(int refSupport, int altSupport)
    {
        if(min(refSupport, altSupport) >= LOH_REF_ALT_MIN_COUNT)
        {
            double totalCount = refSupport + altSupport;
            double baf = max(refSupport, altSupport) / totalCount;

            if(baf < LOH_REF_ALT_THRESHOLD_HIGH)
            {
                ValidBafHigh = true;

                if(baf < LOH_REF_ALT_THRESHOLD_MED)
                {
                    ValidBafMedium = true;

                    if(baf < LOH_REF_ALT_THRESHOLD_LOW)
                    {
                        ValidBafLow = true;
                    }
                }
            }
        }
    }

    public void setValidBafData(int validCountHigh, int validCountMed, int validCountLow)
    {
        ValidRateHigh = validCountHigh / (double)LOH_MOVING_ITEM_COUNT;
        ValidRateMedium = validCountMed / (double)LOH_MOVING_ITEM_COUNT;
        ValidRateLow = validCountLow / (double)LOH_MOVING_ITEM_COUNT;

        if(ValidRateHigh < VALID_BAF_RATE_HIGH && ValidRateMedium < VALID_BAF_RATE_MED && ValidRateLow <= VALID_BAF_RATE_LOW)
        {
            LohCandidate = true;
        }
    }

    public String toString()
    {
        return format("valid(h=%s m=%s l=%s) avgRate(h=%.2f m=%.3f l=%.2f) loh(candidate=%s recovered=%s)",
                ValidBafHigh, ValidBafMedium, ValidBafLow, ValidRateHigh, ValidRateMedium, ValidRateLow, LohCandidate, LohRecovered);
    }
}
