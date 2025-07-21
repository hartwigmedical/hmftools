package com.hartwig.hmftools.geneutils.paneldesign;

public record ProbeSelectCriteria(
        ProbeEvaluator.Criteria eval,
        ProbeSelector.Strategy select
)
{
}
