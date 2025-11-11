package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;

import java.util.List;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Defines how probes are to be generated, at a high level.
// Useful to have this as data to decouple the probe generation code.
public sealed interface ProbeGenerationSpec
        permits ProbeGenerationSpec.CoverRegion, ProbeGenerationSpec.CoverOneSubregion, ProbeGenerationSpec.CoverOnePosition,
        ProbeGenerationSpec.SingleProbe
{
    // Generate the best acceptable probes to cover an entire region (general purpose).
    record CoverRegion(
            ChrBaseRegion region,
            TargetMetadata metadata,
            ProbeEvaluator.Criteria evalCriteria,
            ProbeSelector.Strategy localSelectStrategy
    ) implements ProbeGenerationSpec
    {
    }

    // Generate the one best acceptable probe contained within the specified region.
    record CoverOneSubregion(
            ChrBaseRegion region,
            TargetMetadata metadata,
            ProbeEvaluator.Criteria evalCriteria,
            ProbeSelector.Strategy selectStrategy
    ) implements ProbeGenerationSpec
    {
    }

    // Generate the one best acceptable probe centered on one of the given positions.
    record CoverOnePosition(
            List<BasePosition> positions,
            TargetMetadata metadata,
            ProbeEvaluator.Criteria evalCriteria,
            ProbeSelector.Strategy selectStrategy
    ) implements ProbeGenerationSpec
    {
    }

    // A single probe with the specified data.
    record SingleProbe(
            SequenceDefinition sequenceDefinition,
            TargetedRange targetedRange,
            TargetMetadata metadata,
            ProbeEvaluator.Criteria evalCriteria
    ) implements ProbeGenerationSpec
    {
        public SingleProbe
        {
            if(sequenceDefinition.baseLength() != PROBE_LENGTH)
            {
                throw new IllegalArgumentException("region length must be equal to probe length");
            }
        }
    }
}
