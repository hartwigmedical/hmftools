package com.hartwig.hmftools.panelbuilder;

public record ProbeSelectCriteria(
        ProbeEvaluator.Criteria eval,
        ProbeSelector.Strategy select
)
{
}
