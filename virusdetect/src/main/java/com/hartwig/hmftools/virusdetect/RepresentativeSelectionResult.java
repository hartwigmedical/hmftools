package com.hartwig.hmftools.virusdetect;

import java.util.List;

// The outcome of selecting a representative contig for each oncology group.
public record RepresentativeSelectionResult(
        List<ContigClassification> classifications
)
{
}
