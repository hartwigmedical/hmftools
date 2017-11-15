package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;

class StructuralVariantLegsFactory {

    @NotNull
    static List<StructuralVariantLegs> create(@NotNull List<StructuralVariant> variants) {
        final List<ModifiableStructuralVariantLegs> result = createLegs(variants);
        final Multimap<GenomePosition, ModifiableStructuralVariantLegs> duplicates = findDuplicates(result);
        for (GenomePosition duplicatePosition : duplicates.keySet()) {
            final StructuralVariantLeg approvedLeg = reduce(duplicatePosition, duplicates.get(duplicatePosition));
            boolean match = false;
            for (ModifiableStructuralVariantLegs legs : duplicates.get(duplicatePosition)) {
                if (legs.start().filter(x -> isDuplicate(duplicatePosition, x)).isPresent()) {
                    if (legs.start().get().equals(approvedLeg)) {
                        match = true;
                    } else {
                        legs.setStart(Optional.empty());
                    }
                }

                if (legs.end().filter(x -> isDuplicate(duplicatePosition, x)).isPresent()) {
                    if (legs.end().get().equals(approvedLeg)) {
                        match = true;
                    } else {
                        legs.setEnd(Optional.empty());
                    }
                }
            }

            if (!match) {
                result.add(ModifiableStructuralVariantLegs.create().setStart(approvedLeg).setEnd(Optional.empty()));
            }
        }

        return new ArrayList<>(result);
    }

    @NotNull
    private static StructuralVariantLeg reduce(@NotNull final GenomePosition position,
            @NotNull final Collection<ModifiableStructuralVariantLegs> legs) {
        final List<StructuralVariantLeg> leg = Lists.newArrayList();
        for (ModifiableStructuralVariantLegs modifiableStructuralVariantLegs : legs) {
            modifiableStructuralVariantLegs.start().filter(x -> isDuplicate(position, x)).ifPresent(leg::add);
            modifiableStructuralVariantLegs.end().filter(x -> isDuplicate(position, x)).ifPresent(leg::add);
        }

        return reduce(leg);
    }

    private static boolean isDuplicate(@NotNull final GenomePosition position, @NotNull final StructuralVariantLeg leg) {
        return position.equals(GenomePositions.create(leg));
    }

    @NotNull
    static StructuralVariantLeg reduce(@NotNull final List<StructuralVariantLeg> legs) {

        double maxPositive = 0;
        double maxNegative = 0;

        for (StructuralVariantLeg leg : legs) {
            if (leg.orientation() == 1) {
                maxPositive = Math.max(maxPositive, leg.vaf());
            } else {
                maxNegative = Math.max(maxNegative, leg.vaf());
            }
        }

        if (Doubles.isZero(maxNegative)) {
            for (StructuralVariantLeg leg : legs) {
                if (Doubles.equal(maxPositive, leg.vaf())) {
                    return leg;
                }
            }
        }

        if (Doubles.isZero(maxPositive)) {
            for (StructuralVariantLeg leg : legs) {
                if (Doubles.equal(maxNegative, leg.vaf())) {
                    return leg;
                }
            }
        }

        int orientation = Doubles.greaterThan(maxPositive, maxNegative) ? 1 : -1;
        double vaf = Math.abs(maxPositive - maxNegative);
        return ImmutableStructuralVariantLeg.builder().from(legs.get(0)).orientation(orientation).vaf(vaf).build();
    }

    @NotNull
    private static Multimap<GenomePosition, ModifiableStructuralVariantLegs> findDuplicates(
            @NotNull final Collection<ModifiableStructuralVariantLegs> legs) {

        final ListMultimap<GenomePosition, ModifiableStructuralVariantLegs> result = ArrayListMultimap.create();
        for (ModifiableStructuralVariantLegs leg : legs) {
            leg.start().ifPresent(x -> result.put(GenomePositions.create(x), leg));
            leg.end().ifPresent(x -> result.put(GenomePositions.create(x), leg));
        }

        result.keySet().removeIf(key -> result.get(key).size() <= 1);
        return result;
    }

    @NotNull
    private static List<ModifiableStructuralVariantLegs> createLegs(@NotNull final List<StructuralVariant> variants) {
        final List<ModifiableStructuralVariantLegs> result = Lists.newArrayList();

        for (StructuralVariant variant : variants) {

            if (variant.type() != StructuralVariantType.INS) {

                final Optional<StructuralVariantLeg> start = Optional.ofNullable(variant.startAF())
                        .map(startVaf -> create(variant.startChromosome(), variant.startPosition(), variant.startOrientation(), startVaf));

                final Optional<StructuralVariantLeg> end = Optional.ofNullable(variant.endAF())
                        .map(endVaf -> create(variant.endChromosome(), variant.endPosition(), variant.endOrientation(), endVaf));

                if (start.isPresent() || end.isPresent()) {
                    final ModifiableStructuralVariantLegs legs = ModifiableStructuralVariantLegs.create().setStart(start).setEnd(end);
                    result.add(legs);
                }
            }
        }
        return result;

    }

    private static StructuralVariantLeg create(@NotNull String chromosome, long position, int orientation, double vaf) {
        return ImmutableStructuralVariantLeg.builder().chromosome(chromosome).position(position).orientation(orientation).vaf(vaf).build();
    }
}
