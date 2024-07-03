package com.hartwig.hmftools.purple.copynumber;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegCopyNumber;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegCopyNumberChangeFactory;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegCopyNumberFactory;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidyFactory;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.Nullable;

public final class CopyNumberEnrichedStructuralVariantFactory
{
    private final PurityAdjuster mPurityAdjuster;
    private final Multimap<Chromosome,PurpleCopyNumber> mCopyNumbers;

    public CopyNumberEnrichedStructuralVariantFactory(
            final PurityAdjuster purityAdjuster, final Multimap<Chromosome, PurpleCopyNumber> copyNumbers)
    {
        mPurityAdjuster = purityAdjuster;
        mCopyNumbers = copyNumbers;
    }

    public List<EnrichedStructuralVariant> enrich(final List<StructuralVariant> variants)
    {
        final StructuralVariantLegCopyNumberChangeFactory changeFactory = new StructuralVariantLegCopyNumberChangeFactory(
                mPurityAdjuster, mCopyNumbers, variants);

        final StructuralVariantLegPloidyFactory<PurpleCopyNumber> ploidyFactory = new StructuralVariantLegPloidyFactory<>(
                mPurityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);

        final StructuralVariantLegCopyNumberFactory<PurpleCopyNumber> copyNumberFactory = new StructuralVariantLegCopyNumberFactory<>(
                PurpleCopyNumber::averageTumorCopyNumber);

        final List<EnrichedStructuralVariant> result = Lists.newArrayList();
        for(final StructuralVariant variant : variants)
        {
            ImmutableEnrichedStructuralVariant.Builder builder = ImmutableEnrichedStructuralVariant.builder().from(variant);
            ImmutableEnrichedStructuralVariantLeg.Builder startBuilder = createBuilder(variant.start());

            @Nullable
            final StructuralVariantLeg endLeg = variant.end();
            ImmutableEnrichedStructuralVariantLeg.Builder endBuilder = createBuilder(endLeg);

            List<StructuralVariantLegPloidy> ploidies = ploidyFactory.create(variant, mCopyNumbers);
            if(!ploidies.isEmpty())
            {
                // The implied ploidy should be equal between start and end, so doesn't matter what we pick.
                builder.junctionCopyNumber((ploidies.get(0).averageImpliedPloidy()));

                StructuralVariantLegPloidy startPloidy = ploidies.get(0);
                StructuralVariantLegPloidy endPloidy = ploidies.size() <= 1 ? null : ploidies.get(1);

                startBuilder.adjustedAlleleFrequency((startPloidy.adjustedVaf()));
                startBuilder.adjustedCopyNumber((startPloidy.adjustedCopyNumber()));
                startBuilder.adjustedCopyNumberChange((changeFactory.copyNumberChange(startPloidy)));

                if(endPloidy != null)
                {
                    endBuilder.adjustedAlleleFrequency((endPloidy.adjustedVaf()));
                    endBuilder.adjustedCopyNumber((endPloidy.adjustedCopyNumber()));
                    endBuilder.adjustedCopyNumberChange((changeFactory.copyNumberChange(endPloidy)));
                }
            }
            else
            {
                // Can't always get ploidies (if no VAF for example) but we can still get copy number info
                final StructuralVariantLegCopyNumber startCopyNumber = copyNumberFactory.create(variant.start(), mCopyNumbers);
                startBuilder.adjustedCopyNumber((startCopyNumber.adjustedCopyNumber()));
                startBuilder.adjustedCopyNumberChange((changeFactory.copyNumberChange(startCopyNumber)));

                // Lacking anything else, inferred variants can use copy number change as ploidy
                if(variant.type() == StructuralVariantType.INF)
                {
                    builder.junctionCopyNumber((changeFactory.copyNumberChange(startCopyNumber)));
                }

                if(endLeg != null)
                {
                    final StructuralVariantLegCopyNumber endCopyNumber = copyNumberFactory.create(endLeg, mCopyNumbers);
                    endBuilder.adjustedCopyNumber((endCopyNumber.adjustedCopyNumber()));
                    endBuilder.adjustedCopyNumberChange((changeFactory.copyNumberChange(endCopyNumber)));
                }
            }

            result.add(builder.start(startBuilder.build()).end(endBuilder == null ? null : endBuilder.build()).build());
        }

        return result;
    }

    @Nullable
    private ImmutableEnrichedStructuralVariantLeg.Builder createBuilder(@Nullable StructuralVariantLeg leg)
    {
        if(leg == null)
            return null;

        return ImmutableEnrichedStructuralVariantLeg.builder().from(leg);
    }
}
