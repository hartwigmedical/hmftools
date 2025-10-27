package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.VARIANT_NOVEL_SEQUENCE_BASES_MIN;

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

    static boolean needsCoverageCheck(final SequenceDefinition sequenceDefinition)
    {
        return sequenceDefinition.isExactRegion() || variantProbeNeedsCoverageCheck(sequenceDefinition);
    }

    private static boolean variantProbeNeedsCoverageCheck(final SequenceDefinition definition)
    {
        // Only do the coverage check for variants where the probe is similar to the ref genome.
        // If the probe is similar (e.g. SNV) then that region could be captured by probing the ref genome sequence,
        // so the variant probe is not needed.
        // If the probe is dissimilar (e.g. large INDEL or SV) then we need the variant probe to capture the variant.

        ChrBaseRegion start = definition.startRegion();
        ChrBaseRegion end = definition.endRegion();
        int insertLength = definition.insertSequence() == null ? 0 : definition.insertSequence().length();
        if(start == null && end == null)
        {
            // Unknown region, assume novel sequence.
            return true;
        }
        else if(start != null && end != null)
        {
            if(start.chromosome().equals(end.chromosome()))
            {
                // SNV, INDEL, or SV on same chromosome.
                if(start.start() > end.start())
                {
                    // Ensure start and end are ordered correctly to calculate the delete length.
                    start = definition.endRegion();
                    end = definition.startRegion();
                }
                // Clamp to >=0 because theoretically the regions could overlap in the case of an SV.
                int deleteLength = max(end.start() - start.end() - 1, 0);
                int difference = abs(insertLength - deleteLength);
                return difference >= VARIANT_NOVEL_SEQUENCE_BASES_MIN;
            }
            else
            {
                // SV across different chromosomes.
                return true;
            }
        }
        else
        {
            // Single ended SV.
            return true;
        }
    }
}
