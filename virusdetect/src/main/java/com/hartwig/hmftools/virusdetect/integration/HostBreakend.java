package com.hartwig.hmftools.virusdetect.integration;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

// Where a candidate integration joins the host genome, as ESVEE called it.
public record HostBreakend(
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
                leg.tumorVariantFragmentCount(), leg.tumorReferenceFragmentCount(), leg.alleleFrequency(),
                leg.normalVariantFragmentCount(), leg.normalReferenceFragmentCount());

        return new HostBreakend(leg.chromosome(), leg.position(), Orientation.fromByte(leg.orientation()), support);
    }
}
