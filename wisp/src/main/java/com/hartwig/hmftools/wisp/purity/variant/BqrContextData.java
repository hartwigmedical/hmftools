package com.hartwig.hmftools.wisp.purity.variant;

import static java.lang.Math.round;
import static java.lang.String.format;

import com.hartwig.hmftools.common.qual.BaseQualAdjustment;

public class BqrContextData
{
    public final String TrinucleotideContext;
    public final String Alt;

    public int AltCount;
    public int TotalCount;

    public BqrContextData(final String trinucleotideContext, final String alt)
    {
        TrinucleotideContext = trinucleotideContext;
        Alt = alt;

        AltCount = 0;
        TotalCount = 0;
    }

    private static final int PER_MILLION = 1_000_000;

    public double errorRatePerMillion()
    {
        return TotalCount > 0 ? round(100.0 * AltCount / TotalCount * PER_MILLION) / 100.0 : 0;
    }
    public double errorRate() { return TotalCount > 0 ? AltCount / (double)TotalCount : 0;  }

    public double calculatedQual()
    {
        return TotalCount > 0 ? BaseQualAdjustment.probabilityToPhredQual(AltCount / (double) TotalCount) : 0;
    }

    public String toString()
    {
        return format("key(%s:%s) counts(total=%d alt=%d) errorPerM(%.4f) calcQual(%.4f)",
                TrinucleotideContext, Alt, TotalCount, AltCount, errorRatePerMillion(), calculatedQual());
    }
}
