package com.hartwig.hmftools.purple.copynumber.sv;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

import org.jetbrains.annotations.Nullable;

public class StructuralVariantLegCopyNumberChangeFactory
{
    private final Map<GenomePosition, CopyNumberChange> mCopyNumberChangeMap;

    public StructuralVariantLegCopyNumberChangeFactory(final PurityAdjuster purityAdjuster,
            final Multimap<Chromosome, PurpleCopyNumber> copyNumbers, final Collection<StructuralVariant> variants)
    {
        mCopyNumberChangeMap = copyNumberChangeMap(purityAdjuster, copyNumbers, variants);
    }

    public double copyNumberChange(final StructuralVariantLegPloidy leg)
    {
        GenomePosition cnaPosition = GenomePositions.create(leg.chromosome(), leg.cnaPosition());

        if(!mCopyNumberChangeMap.containsKey(cnaPosition))
        {
            return CopyNumberChange.copyNumberChangeSimple(leg);
        }

        return mCopyNumberChangeMap.get(cnaPosition).copyNumberChange(leg);
    }

    @Nullable
    public Double copyNumberChange(final StructuralVariantLegCopyNumber copyNumber)
    {
        if(!mCopyNumberChangeMap.containsKey(GenomePositions.create(copyNumber.chromosome(), copyNumber.cnaPosition())))
        {
            return CopyNumberChange.copyNumberChangeSimple(copyNumber);
        }

        return null;
    }

    private static Map<GenomePosition, CopyNumberChange> copyNumberChangeMap(
            final PurityAdjuster purityAdjuster, final Multimap<Chromosome, PurpleCopyNumber> copyNumbers,
            final Collection<StructuralVariant> variants)
    {

        final StructuralVariantLegPloidyFactory<PurpleCopyNumber> ploidyFactory =
                new StructuralVariantLegPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);

        final ListMultimap<GenomePosition, StructuralVariantLegPloidy> ploidyMap = ArrayListMultimap.create();
        for(final StructuralVariant variant : variantsAtDuplicateLocations(variants))
        {
            List<StructuralVariantLegPloidy> ploidies = ploidyFactory.create(variant, copyNumbers);
            for(StructuralVariantLegPloidy ploidy : ploidies)
            {
                ploidyMap.put(GenomePositions.create(ploidy.chromosome(), ploidy.cnaPosition()), ploidy);
            }
        }

        Map<GenomePosition, CopyNumberChange> copyNumberChangeMap = Maps.newHashMap();
        for(GenomePosition genomePosition : ploidyMap.keySet())
        {
            final List<StructuralVariantLegPloidy> legs = ploidyMap.get(genomePosition);
            if(legs.size() > 1)
            {
                copyNumberChangeMap.put(genomePosition, new CopyNumberChange(legs));
            }
        }

        return copyNumberChangeMap;
    }

    private static Set<StructuralVariant> variantsAtDuplicateLocations(final Collection<StructuralVariant> variants)
    {
        Set<StructuralVariant> result = Sets.newHashSet();
        Multimap<GenomePosition, StructuralVariant> cnaMap = cnaMap(variants);
        for(GenomePosition genomePosition : cnaMap.keySet())
        {
            final Collection<StructuralVariant> variantsAtLocation = cnaMap.get(genomePosition);
            if(variantsAtLocation.size() > 1)
            {
                result.addAll(variantsAtLocation);
            }

        }
        return result;
    }

    private static ListMultimap<GenomePosition,StructuralVariant> cnaMap(final Collection<StructuralVariant> variants)
    {
        final ListMultimap<GenomePosition, StructuralVariant> result = ArrayListMultimap.create();
        for(final StructuralVariant variant : variants)
        {
            final StructuralVariantLeg start = variant.start();
            result.put(GenomePositions.create(start.chromosome(), start.cnaPosition()), variant);

            @Nullable
            final StructuralVariantLeg end = variant.end();
            if(end != null)
            {
                result.put(GenomePositions.create(end.chromosome(), end.cnaPosition()), variant);
            }
        }

        return result;

    }
}
