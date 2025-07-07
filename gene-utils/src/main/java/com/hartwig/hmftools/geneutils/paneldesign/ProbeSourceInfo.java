package com.hartwig.hmftools.geneutils.paneldesign;

// Information about what caused a probe to be generated.
public record ProbeSourceInfo(
        ProbeSource source,
        String extra
)
{
}
