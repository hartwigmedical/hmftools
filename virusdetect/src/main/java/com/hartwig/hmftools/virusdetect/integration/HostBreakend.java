package com.hartwig.hmftools.virusdetect.integration;

import static java.util.Objects.requireNonNull;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

// Where a candidate viral integration joins the host genome, as ESVEE called it.
public record HostBreakend(
        // TODO: use hmf-common BasePosition
        String chromosome,
        int position,
        Orientation orientation,
        BreakendSupport support
)
{
    public HostBreakend
    {
        if(chromosome.isEmpty())
        {
            throw new IllegalArgumentException("breakend chromosome is empty");
        }
        if(position < 1)
        {
            throw new IllegalArgumentException("invalid breakend position: " + position);
        }
    }

    public static HostBreakend from(StructuralVariantLeg leg)
    {
        BreakendSupport support = new BreakendSupport(
                // The tumor counts are always present, the tumor genotype having been resolved by name before any SV was built.
                requireNonNull(leg.tumorVariantFragmentCount()),
                requireNonNull(leg.tumorReferenceFragmentCount()),
                requireNonNull(leg.alleleFrequency()),
                leg.normalVariantFragmentCount(),
                leg.normalReferenceFragmentCount());

        return new HostBreakend(leg.chromosome(), leg.position(), Orientation.fromByte(leg.orientation()), support);
    }
}
