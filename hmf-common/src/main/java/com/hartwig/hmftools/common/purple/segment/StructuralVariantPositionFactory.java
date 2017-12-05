package com.hartwig.hmftools.common.purple.segment;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;

final class StructuralVariantPositionFactory {

    private static final long MIN_BASES = 1000;

    @NotNull
    public static List<StructuralVariantPosition> create(@NotNull final List<StructuralVariant> variants) {
        final List<StructuralVariantPosition> positions = Lists.newArrayList();
        for (StructuralVariant variant : variants) {

            if (variant.type() != StructuralVariantType.INS) {
                positions.add(ImmutableStructuralVariantPosition.builder()
                        .id(variant.id())
                        .chromosome(variant.startChromosome())
                        .position(variant.startPosition())
                        .type(variant.type())
                        .orientation(variant.startOrientation())
                        .alleleFrequency(variant.startAF())
                        .build());

                positions.add(ImmutableStructuralVariantPosition.builder()
                        .id(variant.id())
                        .chromosome(variant.endChromosome())
                        .position(variant.endPosition())
                        .type(variant.type())
                        .orientation(variant.endOrientation())
                        .alleleFrequency(variant.endAF())
                        .build());
            }
        }

        Collections.sort(positions);
        return positions;
    }

    private static boolean include(@NotNull StructuralVariant variant) {
        return !variant.startChromosome().equals(variant.endChromosome()) || bases(variant) > MIN_BASES;
    }

    private static long bases(@NotNull StructuralVariant variant) {
        return 1 + variant.endPosition() - variant.startPosition();
    }
}
