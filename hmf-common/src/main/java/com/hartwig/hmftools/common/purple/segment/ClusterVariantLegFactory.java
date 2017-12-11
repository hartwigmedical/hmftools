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
                positions.add(ImmutableClusterVariantLeg.builder().from(start).type(variant.type()).build());

                final StructuralVariantLeg end = variant.end();
                positions.add(ImmutableClusterVariantLeg.builder().from(end).type(variant.type()).build());
            }
        }

        Collections.sort(positions);
        return positions;
    }
}
