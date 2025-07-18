package com.hartwig.hmftools.geneutils.paneldesign;

// TODO: put this inside ProbeSelector
// When there are multiple acceptable candidate probes and only 1 probe is desired, how to select the best?
public sealed interface ProbeSelectStrategy
        permits ProbeSelectStrategy.FirstAcceptable, ProbeSelectStrategy.MaxQuality, ProbeSelectStrategy.BestGc
{
    // Pick the first probe that is acceptable.
    record FirstAcceptable() implements ProbeSelectStrategy
    {
    }

    // Pick the acceptable probe with the highest quality score.
    record MaxQuality(
            // Consider the max quality to have been found if exceeding this value.
            // Useful to improve runtime performance.
            double optimalQuality
    ) implements ProbeSelectStrategy
    {
        public MaxQuality()
        {
            this(1);
        }
    }

    // Pick the acceptable probe with the GC content closest to gcContentTarget.
    record BestGc(
            // Consider the best GC to have been found if within this tolerance.
            // Useful to improve runtime performance.
            double gcToleranceOptimal
    ) implements ProbeSelectStrategy
    {
        public BestGc()
        {
            this(1);
        }
    }
}
