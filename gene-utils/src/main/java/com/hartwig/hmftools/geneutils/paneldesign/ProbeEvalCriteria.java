package com.hartwig.hmftools.geneutils.paneldesign;

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
    public double gcContentMin() {
        return gcContentTarget - gcContentTolerance;
    }

    public double gcContentMax() {
        return gcContentTarget + gcContentTolerance;
    }
}
