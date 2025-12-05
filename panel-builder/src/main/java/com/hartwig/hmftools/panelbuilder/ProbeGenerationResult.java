package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.panelbuilder.ProbeUtils.probeTargetRegions;

import java.util.List;
import java.util.stream.Stream;

public record ProbeGenerationResult(
        List<Probe> probes,
        // Regions which were potentially targeted to be covered (and may or may not be covered). For informational purposes only.
        List<TargetRegion> candidateTargetRegions,
        // Features which could not be covered for various reasons. For informational purposes only.
        List<RejectedFeature> rejectedFeatures
        // There is no exact relationship between the types of regions stored here. Don't try to calculate anything based off them.
)
{
    public ProbeGenerationResult
    {
        if(!probes.stream().allMatch(Probe::accepted))
        {
            throw new IllegalArgumentException("Should only add accepted probes to result");
        }
    }

    public ProbeGenerationResult()
    {
        this(emptyList(), emptyList(), emptyList());
    }

    public ProbeGenerationResult add(ProbeGenerationResult other)
    {
        // Naive merge, does not try to deduplicate data.
        return new ProbeGenerationResult(
                Stream.concat(probes.stream(), other.probes.stream()).toList(),
                Stream.concat(candidateTargetRegions.stream(), other.candidateTargetRegions.stream()).toList(),
                Stream.concat(rejectedFeatures.stream(), other.rejectedFeatures.stream()).toList()
        );
    }

    // Convenience method for creating a result from candidate targets which got no probes since they were already covered by the panel.
    public static ProbeGenerationResult alreadyCoveredTargets(final List<TargetRegion> targets)
    {
        return new ProbeGenerationResult(emptyList(), List.copyOf(targets), emptyList());
    }

    // Convenience method for creating a result from rejecting multiple entire target regions.
    public static ProbeGenerationResult rejectTargets(final List<TargetRegion> targets)
    {
        List<RejectedFeature> rejectedFeatures = targets.stream()
                .map(target -> RejectedFeature.fromRegion(target.region(), target.metadata()))
                .toList();
        return new ProbeGenerationResult(emptyList(), List.copyOf(targets), rejectedFeatures);
    }

    // Convenience method for creating a result from a single probe whose regions were already covered by the panel.
    public static ProbeGenerationResult alreadyCoveredProbe(final Probe probe)
    {
        return ProbeGenerationResult.alreadyCoveredTargets(probeTargetRegions(probe));
    }

    // Convenience method for creating a result from a single accepted probe.
    public static ProbeGenerationResult acceptProbe(final Probe probe)
    {
        return new ProbeGenerationResult(List.of(probe), probeTargetRegions(probe), emptyList());
    }

    // Convenience method for creating a result from a single rejected probe.
    public static ProbeGenerationResult rejectProbe(final Probe probe)
    {
        // If the probe is entirely determined by a single reference genome region, add the region as rejected.
        // Otherwise, if it's a variant probe, add the full probe so the user has more information to inspect.
        if(probe.definition().isSingleRegion())
        {
            return ProbeGenerationResult.rejectTargets(probeTargetRegions(probe));
        }
        else
        {
            return new ProbeGenerationResult(emptyList(), probeTargetRegions(probe), List.of(RejectedFeature.fromProbe(probe)));
        }
    }
}
