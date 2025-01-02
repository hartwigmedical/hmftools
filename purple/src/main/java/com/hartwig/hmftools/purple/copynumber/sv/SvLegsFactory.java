package com.hartwig.hmftools.purple.copynumber.sv;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

public final class SvLegsFactory
{
    public static List<StructuralVariantLegs> create(final StructuralVariant variants)
    {
        List<StructuralVariantLegs> legs = createLegs(true, Collections.singletonList(variants));
        return legs;
    }

    public static List<StructuralVariantLegs> create(final List<StructuralVariant> variants)
    {
        List<StructuralVariantLegs> result = createLegs(false, variants);
        Set<GenomePosition> duplicatePositions = findDuplicatePositions(result);
        
        for(GenomePosition duplicatePosition : duplicatePositions)
        {
            List<StructuralVariantLegs> duplicates = result.stream()
                    .filter(x -> matches(duplicatePosition, x)).collect(Collectors.toList());

            StructuralVariantLeg approvedLeg = reduce(duplicatePosition, duplicates);
            boolean match = false;

            for(StructuralVariantLegs legs : duplicates)
            {
                StructuralVariantLeg start = legs.start();

                if(start != null && isDuplicate(duplicatePosition, start))
                {
                    if(!match && start.equals(approvedLeg))
                    {
                        match = true;
                    }
                    else
                    {
                        legs.setStart(null);
                    }
                }

                StructuralVariantLeg end = legs.end();
                if(end != null && isDuplicate(duplicatePosition, end))
                {
                    if(!match && end.equals(approvedLeg))
                    {
                        match = true;
                    }
                    else
                    {
                        legs.setEnd(null);
                    }
                }
            }

            if(!match)
            {
                result.add(new StructuralVariantLegs(approvedLeg, null));
            }
        }

        return result.stream().filter(x -> x.hasEither()).collect(Collectors.toList());
    }

    private static boolean matches(final GenomePosition position, final StructuralVariantLegs legs)
    {
        return (legs.start() != null && cnaPosition(legs.start()).equals(position))
            || (legs.end() != null && cnaPosition(legs.end()).equals(position));
    }

    private static StructuralVariantLeg reduce(final GenomePosition position, final Collection<StructuralVariantLegs> legs)
    {
        List<StructuralVariantLeg> leg = Lists.newArrayList();
        
        for(StructuralVariantLegs svLegs : legs)
        {
            if(svLegs.start() != null && isDuplicate(position, svLegs.start()))
                leg.add(svLegs.start());

            if(svLegs.end() != null && isDuplicate(position, svLegs.end()))
                leg.add(svLegs.end());
        }

        return reduce(position, leg);
    }

    private static boolean isDuplicate(final GenomePosition position, final StructuralVariantLeg leg)
    {
        return position.equals(cnaPosition(leg));
    }

    public static StructuralVariantLeg reduce(final GenomePosition cnaPosition, final List<StructuralVariantLeg> legs)
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

    private static Set<GenomePosition> findDuplicatePositions(final Collection<StructuralVariantLegs> legs)
    {
        ListMultimap<GenomePosition, StructuralVariantLegs> result = ArrayListMultimap.create();

        for(StructuralVariantLegs svLegs : legs)
        {
            if(svLegs.start() != null)
                result.put(cnaPosition(svLegs.start()), svLegs);

            if(svLegs.end() != null)
                result.put(cnaPosition(svLegs.end()), svLegs);
        }

        result.keySet().removeIf(key -> result.get(key).size() <= 1);
        return result.keySet();
    }

    private static List<StructuralVariantLegs> createLegs(boolean allowInserts, final List<StructuralVariant> variants)
    {
        final List<StructuralVariantLegs> result = Lists.newArrayList();

        for(StructuralVariant variant : variants)
        {
            if(allowInserts || variant.type() != StructuralVariantType.INS)
            {
                StructuralVariantLeg start = variant.start() != null && variant.start().alleleFrequency() != null ? variant.start() : null;
                StructuralVariantLeg end = variant.end() != null && variant.end().alleleFrequency() != null ? variant.end() : null;

                if(start != null || end != null)
                {
                    StructuralVariantLegs legs = new StructuralVariantLegs(start, end);
                    result.add(legs);
                }
            }
        }

        return result;
    }

    private static GenomePosition cnaPosition(final StructuralVariantLeg leg)
    {
        return GenomePositions.create(leg.chromosome(), leg.cnaPosition());
    }
}
