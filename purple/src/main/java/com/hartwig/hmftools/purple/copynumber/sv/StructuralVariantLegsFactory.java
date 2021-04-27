package com.hartwig.hmftools.purple.copynumber.sv;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;

public final class StructuralVariantLegsFactory
{
    @NotNull
    public static List<StructuralVariantLegs> create(@NotNull final StructuralVariant variants)
    {
        final List<ModifiableStructuralVariantLegs> legs = createLegs(true, Collections.singletonList(variants));
        return legs.stream().map(x -> (StructuralVariantLegs) x).collect(Collectors.toList());
    }

    @NotNull
    public static List<StructuralVariantLegs> create(@NotNull final List<StructuralVariant> variants)
    {
        final List<ModifiableStructuralVariantLegs> result = createLegs(false, variants);
        final Set<GenomePosition> duplicatePositions = findDuplicatePositions(result);
        for(GenomePosition duplicatePosition : duplicatePositions)
        {
            final List<ModifiableStructuralVariantLegs> duplicates =
                    result.stream().filter(x -> matches(duplicatePosition, x)).collect(Collectors.toList());

            final StructuralVariantLeg approvedLeg = reduce(duplicatePosition, duplicates);
            boolean match = false;
            for(ModifiableStructuralVariantLegs legs : duplicates)
            {
                Optional<StructuralVariantLeg> start = legs.start();
                if(start.filter(x -> isDuplicate(duplicatePosition, x)).isPresent())
                {
                    if(!match && start.get().equals(approvedLeg))
                    {
                        match = true;
                    }
                    else
                    {
                        legs.setStart(Optional.empty());
                    }
                }

                Optional<StructuralVariantLeg> end = legs.end();
                if(end.filter(x -> isDuplicate(duplicatePosition, x)).isPresent())
                {
                    if(!match && end.get().equals(approvedLeg))
                    {
                        match = true;
                    }
                    else
                    {
                        legs.setEnd(Optional.empty());
                    }
                }
            }

            if(!match)
            {
                result.add(ModifiableStructuralVariantLegs.create().setStart(approvedLeg).setEnd(Optional.empty()));
            }
        }

        return result.stream().filter(x -> x.start().isPresent() || x.end().isPresent()).collect(Collectors.toList());
    }

    private static boolean matches(@NotNull final GenomePosition position, @NotNull final ModifiableStructuralVariantLegs legs)
    {
        return legs.start().filter(x -> cnaPosition(x).equals(position)).isPresent() || legs.end()
                .filter(x -> cnaPosition(x).equals(position))
                .isPresent();
    }

    @NotNull
    private static StructuralVariantLeg reduce(@NotNull final GenomePosition position,
            @NotNull final Collection<ModifiableStructuralVariantLegs> legs)
    {
        final List<StructuralVariantLeg> leg = Lists.newArrayList();
        for(ModifiableStructuralVariantLegs modifiableStructuralVariantLegs : legs)
        {
            modifiableStructuralVariantLegs.start().filter(x -> isDuplicate(position, x)).ifPresent(leg::add);
            modifiableStructuralVariantLegs.end().filter(x -> isDuplicate(position, x)).ifPresent(leg::add);
        }

        return reduce(position, leg);
    }

    private static boolean isDuplicate(@NotNull final GenomePosition position, @NotNull final StructuralVariantLeg leg)
    {
        return position.equals(cnaPosition(leg));
    }

    @NotNull
    public static StructuralVariantLeg reduce(@NotNull final GenomePosition cnaPosition, @NotNull final List<StructuralVariantLeg> legs)
    {
        double maxPositive = 0;
        double maxNegative = 0;

        for(StructuralVariantLeg leg : legs)
        {
            if(leg.orientation() == 1)
            {
                maxPositive = Math.max(maxPositive, leg.alleleFrequency());
            }
            else
            {
                maxNegative = Math.max(maxNegative, leg.alleleFrequency());
            }
        }

        if(Doubles.isZero(maxNegative))
        {
            for(StructuralVariantLeg leg : legs)
            {
                if(Doubles.equal(maxPositive, leg.alleleFrequency()))
                {
                    return leg;
                }
            }
        }

        if(Doubles.isZero(maxPositive))
        {
            for(StructuralVariantLeg leg : legs)
            {
                if(Doubles.equal(maxNegative, leg.alleleFrequency()))
                {
                    return leg;
                }
            }
        }

        int maxTumorReferenceFragmentCount = 0;
        int cumulativeTumorVariantFragmentCount = 0;
        for(StructuralVariantLeg leg : legs)
        {
            int orientation = leg.orientation();
            Integer tumorVariantFragmentCount = leg.tumorVariantFragmentCount();
            Integer tumorReferenceFragmentCount = leg.tumorReferenceFragmentCount();

            if(tumorVariantFragmentCount != null)
            {
                cumulativeTumorVariantFragmentCount += orientation * tumorVariantFragmentCount;
            }

            if(tumorReferenceFragmentCount != null)
            {
                maxTumorReferenceFragmentCount = Math.max(maxTumorReferenceFragmentCount, tumorReferenceFragmentCount);
            }

        }

        // double cumulativeVaf = 0;
        // cumulativeVaf += orientation * leg.alleleFrequency() / (1 - leg.alleleFrequency());
        // double vaf = cumulativeVaf / (cumulativeVaf + 1);

        byte orientation = (byte) (Doubles.greaterThan(maxPositive, maxNegative) ? 1 : -1);
        double vaf = Math.abs(maxPositive - maxNegative);
        return ImmutableStructuralVariantLegImpl.builder()
                .chromosome(cnaPosition.chromosome())
                .position(orientation == -1 ? cnaPosition.position() : cnaPosition.position() - 1)
                .orientation(orientation)
                .alleleFrequency(Math.abs(vaf))
                .tumorVariantFragmentCount(cumulativeTumorVariantFragmentCount == 0 ? null : Math.abs(cumulativeTumorVariantFragmentCount))
                .tumorReferenceFragmentCount(maxTumorReferenceFragmentCount == 0 ? null : Math.abs(maxTumorReferenceFragmentCount))
                .alleleFrequency(Math.abs(vaf))
                .homology("")
                .anchoringSupportDistance(0)
                .build();
    }

    @NotNull
    private static Set<GenomePosition> findDuplicatePositions(@NotNull final Collection<ModifiableStructuralVariantLegs> legs)
    {
        final ListMultimap<GenomePosition, ModifiableStructuralVariantLegs> result = ArrayListMultimap.create();
        for(ModifiableStructuralVariantLegs leg : legs)
        {
            leg.start().ifPresent(x -> result.put(cnaPosition(x), leg));
            leg.end().ifPresent(x -> result.put(cnaPosition(x), leg));
        }

        result.keySet().removeIf(key -> result.get(key).size() <= 1);
        return result.keySet();
    }

    @NotNull
    private static List<ModifiableStructuralVariantLegs> createLegs(boolean allowInserts, @NotNull final List<StructuralVariant> variants)
    {
        final List<ModifiableStructuralVariantLegs> result = Lists.newArrayList();

        for(StructuralVariant variant : variants)
        {

            if(allowInserts || variant.type() != StructuralVariantType.INS)
            {
                final Optional<StructuralVariantLeg> start = Optional.of(variant.start()).filter(x -> x.alleleFrequency() != null);
                final Optional<StructuralVariantLeg> end = Optional.ofNullable(variant.end()).filter(x -> x.alleleFrequency() != null);

                if(start.isPresent() || end.isPresent())
                {
                    final ModifiableStructuralVariantLegs legs = ModifiableStructuralVariantLegs.create().setStart(start).setEnd(end);
                    result.add(legs);
                }
            }
        }

        return result;
    }

    @NotNull
    private static GenomePosition cnaPosition(@NotNull final StructuralVariantLeg leg)
    {
        return GenomePositions.create(leg.chromosome(), leg.cnaPosition());
    }

}
