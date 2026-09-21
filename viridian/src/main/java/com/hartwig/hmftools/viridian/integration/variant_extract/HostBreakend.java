package com.hartwig.hmftools.viridian.integration.variant_extract;

import static java.util.Objects.requireNonNull;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

// Where a candidate viral integration joins the host genome, as ESVEE called it.
public record HostBreakend(
        // VCF ID field.
        String id,
        BasePosition basePosition,
        Orientation orientation,
        BreakendSupport support
)
{
    public HostBreakend
    {
        if(id.isEmpty())
        {
            throw new IllegalArgumentException("id cannot be empty");
        }
    }

    public static HostBreakend from(String id, StructuralVariantLeg leg)
    {
        BreakendSupport support = new BreakendSupport(
                // The tumor counts are always present, the tumor genotype having been resolved by name before any SV was built.
                requireNonNull(leg.tumorVariantFragmentCount()),
                requireNonNull(leg.tumorReferenceFragmentCount()),
                requireNonNull(leg.alleleFrequency()),
                leg.normalVariantFragmentCount(),
                leg.normalReferenceFragmentCount());
        return new HostBreakend(id, new BasePosition(leg.chromosome(), leg.position()), Orientation.fromByte(leg.orientation()), support);
    }
}
