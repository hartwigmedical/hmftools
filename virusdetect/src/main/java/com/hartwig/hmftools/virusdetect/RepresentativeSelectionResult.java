package com.hartwig.hmftools.virusdetect;

import java.util.List;
import java.util.Map;

// The outcome of selecting a representative contig for each oncology group.
public record RepresentativeSelectionResult(
        List<ContigClassification> classifications,
        Map<String, Double> oncologyGroupVoteTotals
)
{
}
