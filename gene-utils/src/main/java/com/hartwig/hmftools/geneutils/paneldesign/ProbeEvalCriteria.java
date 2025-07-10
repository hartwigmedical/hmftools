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
        DecimalFormat df = new DecimalFormat();
        df.setMinimumFractionDigits(0);
        return format("quality>=%s gc=%s+-%s", df.format(qualityScoreMin), df.format(gcContentTarget), df.format(gcContentTolerance));
    }
}
