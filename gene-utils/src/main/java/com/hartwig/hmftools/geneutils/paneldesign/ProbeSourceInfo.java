package com.hartwig.hmftools.geneutils.paneldesign;

// Information about what caused a probe to be generated.
public record ProbeSourceInfo(
        ProbeSourceType type,
        String extra
)
{
}
