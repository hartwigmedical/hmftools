package com.hartwig.hmftools.common.purple.segment;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;

final class StructuralVariantPositionFactory {

    @NotNull
    public static List<StructuralVariantPosition> create(@NotNull final List<StructuralVariant> variants) {
        final List<StructuralVariantPosition> positions = Lists.newArrayList();
        for (StructuralVariant variant : variants) {

            if (variant.type() != StructuralVariantType.INS) {
                final StructuralVariantLeg start = variant.start();
                positions.add(ImmutableStructuralVariantPosition.builder()
                        .from(start)
                        .type(variant.type())
                        .orientation(start.orientation())
                        .alleleFrequency(start.alleleFrequency())
                        .build());

                final StructuralVariantLeg end = variant.end();
                positions.add(ImmutableStructuralVariantPosition.builder()
                        .from(end)
                        .type(variant.type())
                        .orientation(end.orientation())
                        .alleleFrequency(end.alleleFrequency())
                        .build());
            }
        }

        Collections.sort(positions);
        return positions;
    }

}
