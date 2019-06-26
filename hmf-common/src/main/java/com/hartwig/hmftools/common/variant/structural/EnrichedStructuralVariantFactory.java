package com.hartwig.hmftools.common.variant.structural;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;

public final class EnrichedStructuralVariantFactory {

    private static final int DISTANCE = 10;
    private final static String PURPLE_AF = "PURPLE_AF";
    private final static String PURPLE_CN = "PURPLE_CN";
    private final static String PURPLE_CN_CHANGE = "PURPLE_CN_CHANGE";
    private final static String PURPLE_PLOIDY = "PURPLE_PLOIDY";

    @Nullable
    private final IndexedFastaSequenceFile reference;

    public EnrichedStructuralVariantFactory(@Nullable final IndexedFastaSequenceFile reference) {
        this.reference = reference;
    }

    @NotNull
    public List<EnrichedStructuralVariant> enrich(@NotNull final List<StructuralVariant> variants) {
        final List<EnrichedStructuralVariant> result = Lists.newArrayList();

        for (StructuralVariant variant : variants) {
            final VariantContext context = variant.startContext();
            if (context != null) {

                final Double purplePloidy = context.hasAttribute(PURPLE_PLOIDY) ? context.getAttributeAsDouble(PURPLE_PLOIDY, 0) : null;

                final ImmutableEnrichedStructuralVariant.Builder builder = ImmutableEnrichedStructuralVariant.builder()
                        .from(variant)
                        .ploidy(purplePloidy)
                        .start(createBuilder(context, 0, variant.start()));

                @Nullable
                final StructuralVariantLeg endLeg = variant.end();
                if (endLeg != null) {
                    builder.end(createBuilder(context, 1, endLeg));
                }

                result.add(builder.build());
            }
        }

        return result;
    }

    @NotNull
    private ImmutableEnrichedStructuralVariantLeg createBuilder(@NotNull final VariantContext context, int index,
            @NotNull final StructuralVariantLeg leg) {

        final List<Double> purpleAF =
                context.hasAttribute(PURPLE_AF) ? context.getAttributeAsDoubleList(PURPLE_AF, 0.0) : Collections.emptyList();

        final List<Double> purpleCN =
                context.hasAttribute(PURPLE_CN) ? context.getAttributeAsDoubleList(PURPLE_CN, 0.0) : Collections.emptyList();

        final List<Double> purpleCNChange =
                context.hasAttribute(PURPLE_CN_CHANGE) ? context.getAttributeAsDoubleList(PURPLE_CN_CHANGE, 0.0) : Collections.emptyList();

        final ImmutableEnrichedStructuralVariantLeg.Builder builder = ImmutableEnrichedStructuralVariantLeg.builder()
                .from(leg)
                .refGenomeContext(context(leg.chromosome(), leg.position()))
                .adjustedAlleleFrequency(purpleAF.size() > index ? purpleAF.get(index) : null)
                .adjustedCopyNumber(purpleCN.size() > index ? purpleCN.get(index) : null)
                .adjustedCopyNumberChange(purpleCNChange.size() > index ? purpleCNChange.get(index) : null);
        return builder.build();
    }

    @NotNull
    private String context(@NotNull String chromosome, long position) {
        if (reference == null) {
            return Strings.EMPTY;
        }

        final int chromosomeLength = reference.getSequenceDictionary().getSequence(chromosome).getSequenceLength();
        final ReferenceSequence sequence =
                reference.getSubsequenceAt(chromosome, Math.max(1, position - DISTANCE), Math.min(position + DISTANCE, chromosomeLength));
        return sequence.getBaseString();
    }
}
