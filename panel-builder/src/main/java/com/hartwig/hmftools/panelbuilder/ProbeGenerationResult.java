package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyList;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionIntersection;

import java.util.List;
import java.util.stream.Stream;

public record ProbeGenerationResult(
        List<Probe> probes,
        // Regions which where potentially targeted to be covered (and may or may not be covered).
        List<TargetRegion> candidateTargetRegions,
        // Regions which targeted to be covered and are covered.
        List<TargetRegion> coveredTargetRegions,
        // Regions which could not be covered due for various reasons.
        List<RejectedRegion> rejectedRegions
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

    // Convenience method for creating a result from covering a target region with exactly one probe.
    public static ProbeGenerationResult coveredTarget(final TargetRegion candidateTarget, final Probe probe)
    {
        TargetRegion covered = new TargetRegion(
                regionIntersection(candidateTarget.region(), requireNonNull(probe.region())).orElseThrow(),
                candidateTarget.metadata());
        return new ProbeGenerationResult(
                List.of(probe),
                List.of(candidateTarget),
                List.of(covered),
                emptyList()
        );
    }

    // Convenience method for creating a result from a candidate target which got no probes since it was already covered.
    public static ProbeGenerationResult alreadyCoveredTarget(final TargetRegion candidateTarget)
    {
        return new ProbeGenerationResult(emptyList(), List.of(candidateTarget), emptyList(), emptyList());
    }

    // Convenience method for creating a result from candidate targets which got no probes since they were already covered.
    public static ProbeGenerationResult alreadyCoveredTargets(final List<TargetRegion> targets)
    {
        return new ProbeGenerationResult(emptyList(), List.copyOf(targets), emptyList(), emptyList());
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

    // Convenience method for creating a result from rejecting multiple entire target regions.
    public static ProbeGenerationResult rejectTargets(final List<TargetRegion> targets, final String rejectionReason)
    {
        return new ProbeGenerationResult(
                emptyList(),
                List.copyOf(targets),
                emptyList(),
                targets.stream().map(target -> RejectedRegion.rejectTarget(target, rejectionReason)).toList());
    }
}
