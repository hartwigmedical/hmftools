package com.hartwig.hmftools.common.variant.structural;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidyFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public final class EnrichedStructuralVariantFactory {

    private static final int DISTANCE = 10;

    private final PurityAdjuster purityAdjuster;
    private final IndexedFastaSequenceFile reference;
    private final Multimap<Chromosome, PurpleCopyNumber> copyNumbers;

    public EnrichedStructuralVariantFactory(@NotNull final IndexedFastaSequenceFile reference, @NotNull final PurityAdjuster purityAdjuster,
            @NotNull final Multimap<Chromosome, PurpleCopyNumber> copyNumbers) {
        this.reference = reference;
        this.purityAdjuster = purityAdjuster;
        this.copyNumbers = copyNumbers;
    }

    @NotNull
    public List<EnrichedStructuralVariant> enrich(@NotNull final List<StructuralVariant> variants) {
        final StructuralVariantLegPloidyFactory<PurpleCopyNumber> ploidyFactory =
                new StructuralVariantLegPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);

        final List<EnrichedStructuralVariant> result = Lists.newArrayList();
        for (final StructuralVariant variant : variants) {
            ImmutableEnrichedStructuralVariant.Builder builder = ImmutableEnrichedStructuralVariant.builder().from(variant);
            ImmutableEnrichedStructuralVariantLeg.Builder startBuilder = createBuilder(variant.start());
            ImmutableEnrichedStructuralVariantLeg.Builder endBuilder = null;

            final List<StructuralVariantLegPloidy> ploidies = ploidyFactory.create(variant, copyNumbers);
            if (!ploidies.isEmpty()) {
                Double roundedPloidy = round(ploidies.get(0).averageImpliedPloidy());
                builder.ploidy(roundedPloidy);
            }

            if (ploidies.size() > 0) {
                final StructuralVariantLegPloidy start = ploidies.get(0);
                final StructuralVariantLegPloidy end = ploidies.size() <= 1 ? null : ploidies.get(1);

                startBuilder.adjustedAlleleFrequency(round(adjustedVAF(purityAdjuster, start)));
                startBuilder.adjustedCopyNumber(round(adjustedCopyNumber(start)));
                startBuilder.adjustedCopyNumberChange(round(adjustedCopyNumberChange(start)));

                if (end != null) {
                    endBuilder = createBuilder(variant.end());
                    endBuilder.adjustedAlleleFrequency(round(adjustedVAF(purityAdjuster, end)));
                    endBuilder.adjustedCopyNumber(round(adjustedCopyNumber(end)));
                    endBuilder.adjustedCopyNumberChange(round(adjustedCopyNumberChange(end)));
                }
            }

            result.add(builder.start(startBuilder.build()).end(endBuilder == null ? null : endBuilder.build()).build());
        }

        return result;
    }

    @NotNull
    private ImmutableEnrichedStructuralVariantLeg.Builder createBuilder(@Nullable StructuralVariantLeg leg) {
        assert leg != null;

        final ImmutableEnrichedStructuralVariantLeg.Builder builder = ImmutableEnrichedStructuralVariantLeg.builder().from(leg);
        builder.refGenomeContext(context(leg.chromosome(), leg.position()));
        return builder;
    }

    @Nullable
    private static Double round(@Nullable Double value) {
        return value == null ? null : Math.round(value * 1000d) / 1000d;
    }

    private static double adjustedVAF(@NotNull final PurityAdjuster purityAdjuster, @NotNull final StructuralVariantLegPloidy ploidy) {
        final Double adjustedCopyNumber = adjustedCopyNumber(ploidy);
        return purityAdjuster.purityAdjustedVAF(ploidy.chromosome(), Math.max(0.001, adjustedCopyNumber), ploidy.vaf());
    }

    @VisibleForTesting
    static double adjustedCopyNumber(@NotNull final StructuralVariantLegPloidy ploidy) {
        if (ploidy.orientation() == 1) {
            return ploidy.leftCopyNumber().orElse(0D);
        } else {
            return ploidy.rightCopyNumber().orElse(0D);
        }
    }

    @VisibleForTesting
    static double adjustedCopyNumberChange(@NotNull final StructuralVariantLegPloidy ploidy) {
        double leftCopyNumber = ploidy.leftCopyNumber().orElse(0D);
        double rightCopyNumber = ploidy.rightCopyNumber().orElse(0D);

        return ploidy.orientation() == 1 ? leftCopyNumber - rightCopyNumber : rightCopyNumber - leftCopyNumber;
    }

    @NotNull
    private String context(@NotNull String chromosome, long position) {
        final int chromosomeLength = reference.getSequenceDictionary().getSequence(chromosome).getSequenceLength();
        final ReferenceSequence sequence =
                reference.getSubsequenceAt(chromosome, Math.max(1, position - DISTANCE), Math.min(position + DISTANCE, chromosomeLength));
        return sequence.getBaseString();
    }
}
