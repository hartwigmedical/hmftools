package com.hartwig.hmftools.purple.segment;

import java.util.Collections;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;

final class SVSegmentFactory {

    @NotNull
    public static List<SVSegment> create(@NotNull final List<StructuralVariant> variants) {
        final List<SVSegment> positions = Lists.newArrayList();
        for (StructuralVariant variant : variants) {

            if (variant.type() != StructuralVariantType.INS) {
                positions.add(create(variant.type(),  variant.start()));
                Optional.ofNullable(variant.end()).map(x -> create(variant.type(), x)).ifPresent(positions::add);
            }
        }

        Collections.sort(positions);
        return positions;
    }

    @NotNull
    private static SVSegment create(@NotNull final StructuralVariantType type, @NotNull final StructuralVariantLeg leg) {
        return ImmutableSVSegment.builder().chromosome(leg.chromosome()).position(leg.cnaPosition()).type(type).build();
    }
}
