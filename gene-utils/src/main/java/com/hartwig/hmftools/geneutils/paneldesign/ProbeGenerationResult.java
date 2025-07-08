package com.hartwig.hmftools.geneutils.paneldesign;

import java.util.Collections;
import java.util.List;
import java.util.stream.Stream;

public record ProbeGenerationResult(
        List<ProbeCandidate> probes,
        List<RejectedRegion> rejectedRegions
)
{
    public ProbeGenerationResult()
    {
        this(Collections.emptyList(), Collections.emptyList());
    }

    public ProbeGenerationResult add(ProbeGenerationResult other)
    {
        // Naive merge, does not try to deduplicate data.
        return new ProbeGenerationResult(
                Stream.concat(probes.stream(), other.probes.stream()).toList(),
                Stream.concat(rejectedRegions.stream(), other.rejectedRegions.stream()).toList()
        );
    }
}
