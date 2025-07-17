package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import java.text.DecimalFormat;

import org.jetbrains.annotations.NotNull;

// Criteria for how to find acceptable probes and choose the best probes.
public record ProbeEvalCriteria(
        // Quality score must be >= this value.
        double qualityScoreMin,
        // Target GC content.
        double gcContentTarget,
        // How much +/- gcContentTarget to accept.
        double gcContentTolerance
)
{
    private static final DecimalFormat DECIMAL_FORMAT;

    static
    {
        DECIMAL_FORMAT = new DecimalFormat();
        DECIMAL_FORMAT.setMinimumFractionDigits(0);
    }

    public ProbeEvalCriteria
    {
        if(!(qualityScoreMin > 0 && qualityScoreMin <= 1))
        {
            // Note quality score is always required, quality=0 is never acceptable.
            throw new IllegalArgumentException("qualityScoreMin must be in range (0, 1]");
        }
        if(!(gcContentMax() >= 0 && gcContentMin() <= 1))
        {
            throw new IllegalArgumentException("GC content range must overlap range [0, 1]");
        }
    }

    public double gcContentMin()
    {
        return gcContentTarget - gcContentTolerance;
    }

    public double gcContentMax()
    {
        return gcContentTarget + gcContentTolerance;
    }

    @NotNull
    @Override
    public String toString()
    {
        String str = format("quality>=%s", DECIMAL_FORMAT.format(qualityScoreMin));
        if(gcContentMin() > 0 || gcContentMax() < 1)
        {
            // Only show GC criteria if it does anything.
            str += format(" gc=%s+-%s", DECIMAL_FORMAT.format(gcContentTarget), DECIMAL_FORMAT.format(gcContentTolerance));
        }
        return str;
    }
}
