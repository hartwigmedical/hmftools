package com.hartwig.hmftools.common.variant.structural;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

public class StructuralVariantPositionFactory {

    private static final long MIN_BASES = 2000;

    public static List<StructuralVariantPosition> create(List<StructuralVariant> variants) {
        final List<StructuralVariantPosition> positions = Lists.newArrayList();
        for (StructuralVariant variant : variants) {
            if (include(variant)) {

                positions.add(ImmutableStructuralVariantPosition.builder()
                        .chromosome(variant.startChromosome())
                        .position(variant.startPosition())
                        .type(variant.type())
                        .build());

                positions.add(ImmutableStructuralVariantPosition.builder()
                        .chromosome(variant.endChromosome())
                        .position(variant.endPosition())
                        .type(variant.type())
                        .build());
            }
        }

        Collections.sort(positions);
        return positions;
    }

    private static boolean include(StructuralVariant variant) {
        return !variant.startChromosome().equals(variant.endChromosome()) || bases(variant) > MIN_BASES;
    }

    private static long bases(StructuralVariant variant) {
        return 1 + variant.endPosition() - variant.startPosition();
    }

}
