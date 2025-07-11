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

    static {
        DECIMAL_FORMAT = new DecimalFormat();
        DECIMAL_FORMAT.setMinimumFractionDigits(0);
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
        return format("quality>=%s gc=%s+-%s",
                DECIMAL_FORMAT.format(qualityScoreMin),
                DECIMAL_FORMAT.format(gcContentTarget),
                DECIMAL_FORMAT.format(gcContentTolerance));
    }
}
