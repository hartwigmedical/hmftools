package com.hartwig.hmftools.purple.segment;

import java.util.Collections;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

final class SVSegmentFactory
{
    public static List<SVSegment> create(final List<StructuralVariant> variants)
    {
        final List<SVSegment> positions = Lists.newArrayList();
        for(StructuralVariant variant : variants)
        {
            if(variant.type() != StructuralVariantType.INS)
            {
                positions.add(create(variant.type(), variant.start()));
                Optional.ofNullable(variant.end()).map(x -> create(variant.type(), x)).ifPresent(positions::add);
            }
        }

        Collections.sort(positions);
        return positions;
    }

    private static SVSegment create(final StructuralVariantType type, final StructuralVariantLeg leg)
    {
        return new SVSegment(leg.chromosome(), leg.cnaPosition(), type);
    }
}
