package com.hartwig.hmftools.panelbuilder;

// Per-gene probe generation statistics, for both DNA and RNA gene probes.
public record GeneStats(
        String geneName,
        int probeCount
)
{
}
