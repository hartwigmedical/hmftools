package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.Utils.regionIntersection;

import java.util.List;
import java.util.stream.Stream;

public record ProbeGenerationResult(
        List<Probe> probes,
        // Regions which potentially targeted to be covered.
        List<TargetRegion> candidateTargetRegions,
        // Regions which targeted to be covered and are covered.
        List<TargetRegion> coveredTargetRegions,
        // Regions which could not be covered due for various reasons.
        List<RejectedRegion> rejectedRegions
)
{
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

    // Convenience method for creating a result from covering a target region with exactly 1 probe.
    public static ProbeGenerationResult coveredTarget(final TargetRegion candidateTarget, final Probe probe)
    {
        TargetRegion covered = new TargetRegion(
                regionIntersection(candidateTarget.region(), probe.region()).orElseThrow(),
                candidateTarget.metadata());
        return new ProbeGenerationResult(
                List.of(probe),
                List.of(candidateTarget),
                List.of(covered),
                emptyList()
        );
    }

    // Convenience method for creating a result from rejecting an entire target region.
    public static ProbeGenerationResult rejectTarget(final TargetRegion target, final String rejectionReason)
    {
        return new ProbeGenerationResult(
                emptyList(),
                List.of(target),
                emptyList(),
                List.of(RejectedRegion.rejectTarget(target, rejectionReason)));
    }
}
