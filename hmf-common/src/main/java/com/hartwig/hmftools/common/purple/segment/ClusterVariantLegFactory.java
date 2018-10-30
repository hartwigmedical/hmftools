package com.hartwig.hmftools.common.purple.segment;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;

final class ClusterVariantLegFactory {

    @NotNull
    public static List<ClusterVariantLeg> create(@NotNull final List<StructuralVariant> variants) {
        final List<ClusterVariantLeg> positions = Lists.newArrayList();
        for (StructuralVariant variant : variants) {

            if (variant.type() != StructuralVariantType.INS) {
                final StructuralVariantLeg start = variant.start();
                positions.add(create(variant.type(), start));

                final StructuralVariantLeg end = variant.end();
                if (end != null) {
                    positions.add(create(variant.type(), end));
                }
            }
        }

        Collections.sort(positions);
        return positions;
    }

    @NotNull
    private static ClusterVariantLeg create(@NotNull final StructuralVariantType type, @NotNull final StructuralVariantLeg leg) {
        long position = leg.orientation() == 1 ? leg.position() - 1 : leg.position();

        position = leg.position();

        return ImmutableClusterVariantLeg.builder().chromosome(leg.chromosome()).position(position).type(type).build();
    }

}
