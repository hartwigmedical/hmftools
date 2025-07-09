package com.hartwig.hmftools.geneutils.paneldesign;

// When there are multiple acceptable candidate probes and only 1 probe is desired, how to select the best?
public enum ProbeSelectStrategy
{
    // Pick the first probe that is acceptable.
    FIRST_ACCEPTABLE,
    // Pick the acceptable probe with the highest quality score.
    MAX_QUALITY,
    // Pick the acceptable probe with the GC content closest to gcContentTarget.
    BEST_GC
}
