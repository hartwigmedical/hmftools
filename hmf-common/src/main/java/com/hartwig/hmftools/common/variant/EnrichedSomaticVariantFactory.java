package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.Builder;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class EnrichedSomaticVariantFactory extends RefGenomeEnrichedSomaticVariantFactory {

    public EnrichedSomaticVariantFactory(@NotNull final IndexedFastaSequenceFile reference) {
        super(reference);
    }

    @NotNull
    @Override
    protected Builder enrich(@NotNull final SomaticVariant variant) {
        return super.enrich(variant);
    }

}
