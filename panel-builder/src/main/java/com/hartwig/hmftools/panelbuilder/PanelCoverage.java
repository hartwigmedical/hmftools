package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.VARIANT_NOVEL_SEQUENCE_BASES_MIN;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.sequenceIndelSize;

import java.util.stream.Stream;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Interface for querying regions already covered by the panel probes.
public interface PanelCoverage
{
    // Checks if a region is fully covered by probes in the panel.
    boolean isCovered(final ChrBaseRegion region);

    // Checks if a probe sequence is fully covered by probes in the panel.
    default boolean isCovered(final SequenceDefinition sequenceDefinition)
    {
        if(needsCoverageCheck(sequenceDefinition))
        {
            return sequenceDefinition.regions().stream().allMatch(this::isCovered);
        }
        else
        {
            // If the coverage check is not applicable (e.g. because it's a large variant) say the probe is not covered, so it gets included
            // in the panel.
            // Theoretically, the probe sequence could already be covered by the panel if multiple probes happen to have the same sequence,
            // but that is so unlikely that we'll assume it doesn't happen.
            return false;
        }
    }

    // Gets all regions covered by probes in the panel.
    Stream<ChrBaseRegion> coveredRegions();

    private static boolean needsCoverageCheck(final SequenceDefinition sequenceDefinition)
    {
        // Only do the coverage check for variants where the probe is similar to the ref genome.
        // If the probe is similar (e.g. SNV) then that region could be captured by probing the ref genome sequence,
        // so the variant probe is not needed.
        // If the probe is dissimilar (e.g. large INDEL or SV) then we need the variant probe to capture the variant.
        return sequenceIndelSize(sequenceDefinition).orElse(Integer.MAX_VALUE) < VARIANT_NOVEL_SEQUENCE_BASES_MIN;
    }
}
