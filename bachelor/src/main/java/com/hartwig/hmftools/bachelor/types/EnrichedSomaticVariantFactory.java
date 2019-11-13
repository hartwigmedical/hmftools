package com.hartwig.hmftools.bachelor.types;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class EnrichedSomaticVariantFactory extends RefGenomeEnrichedSomaticVariantFactory {

    public EnrichedSomaticVariantFactory(@NotNull final IndexedFastaSequenceFile reference) {
        super(reference);
    }

    @NotNull
    @Override
    protected ImmutableEnrichedSomaticVariant.Builder enrich(@NotNull final SomaticVariant variant) {
        return super.enrich(variant);
    }
}
