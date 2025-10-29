package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyList;

import java.util.List;
import java.util.stream.Stream;

public record ProbeGenerationResult(
        List<Probe> probes,
        // Regions which were potentially targeted to be covered (and may or may not be covered). For informational purposes only.
        List<TargetRegion> candidateTargetRegions,
        // TODO: should remove this and store the target subrange in each probe. calculate total target regions later
        // Regions which were targeted to be covered and are covered by probe regions. Guaranteed to be a subset of the probe regions.
        List<TargetRegion> coveredTargetRegions,
        // Regions which could not be covered due for various reasons. For informational purposes only.
        List<RejectedRegion> rejectedRegions
        // There is no exact relationship between the types of regions stored here. Don't try to calculate anything based off
        // candidateTargetRegions or rejectedRegions.
)
{
    public ProbeGenerationResult
    {
        if(probes.stream().anyMatch(probe -> !probe.accepted()))
        {
            throw new IllegalArgumentException("Should only add accepted probes to result");
        }
    }

    public ProbeGenerationResult()
    {
        this(emptyList(), emptyList(), emptyList(), emptyList());
    }

    public ProbeGenerationResult add(ProbeGenerationResult other)
    {
        // Naive merge, does not try to deduplicate data.
        return new ProbeGenerationResult(
                Stream.concat(probes.stream(), other.probes.stream()).toList(),
                Stream.concat(candidateTargetRegions.stream(), other.candidateTargetRegions.stream()).toList(),
                Stream.concat(coveredTargetRegions.stream(), other.coveredTargetRegions.stream()).toList(),
                Stream.concat(rejectedRegions.stream(), other.rejectedRegions.stream()).toList()
        );
    }

    // Convenience method for creating a result from candidate targets which got no probes since they were already covered.
    public static ProbeGenerationResult alreadyCoveredTargets(final List<TargetRegion> targets)
    {
        return new ProbeGenerationResult(emptyList(), List.copyOf(targets), emptyList(), emptyList());
    }

    // Convenience method for creating a result from rejecting multiple entire target regions.
    public static ProbeGenerationResult rejectTargets(final List<TargetRegion> targets)
    {
        List<RejectedRegion> rejectedRegions = targets.stream()
                .map(target -> new RejectedRegion(target.region(), target.metadata()))
                .toList();
        return new ProbeGenerationResult(
                emptyList(),
                List.copyOf(targets),
                emptyList(),
                rejectedRegions);
    }
}
