package com.hartwig.hmftools.panelbuilder;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

// Region which was requested to be covered by probes but couldn't be due to probe selection criteria.
// Just for information purposes for the user.
// Either:
//   - region is provided, indicating a genome region which is rejected; or
//   - probe is provided, indicating a unique probe which is rejected (for probes which are not fully based on the reference genome).
public record RejectedFeature(
        // The exact region that could not be covered.
        @Nullable ChrBaseRegion region,
        // The probe which was rejected.
        @Nullable Probe probe,
        TargetMetadata metadata
)
{
    public RejectedFeature
    {
        if((region == null) == (probe == null))
        {
            throw new IllegalArgumentException("Exactly one of region or probe must be provided");
        }
        if(region != null && !region.hasValidPositions())
        {
            throw new IllegalArgumentException("Invalid region");
        }
    }

    public static RejectedFeature fromRegion(final ChrBaseRegion region, final TargetMetadata metadata)
    {
        return new RejectedFeature(region, null, metadata);
    }

    public static RejectedFeature fromProbe(final Probe probe)
    {
        return new RejectedFeature(null, probe, probe.metadata());
    }
}
