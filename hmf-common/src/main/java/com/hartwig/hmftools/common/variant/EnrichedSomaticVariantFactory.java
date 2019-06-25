package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant.Builder;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class EnrichedSomaticVariantFactory extends RefGenomeEnrichedSomaticVariantFactory {

    @NotNull
    private final GenomeRegionSelector<GenomeRegion> highConfidenceSelector;
    @NotNull
    private final ClonalityFactory clonalityFactory;

    public EnrichedSomaticVariantFactory(@NotNull final Multimap<String, GenomeRegion> highConfidenceRegions,
            @NotNull final IndexedFastaSequenceFile reference, @NotNull final ClonalityFactory clonalityFactory) {
        super(reference);
        this.highConfidenceSelector = GenomeRegionSelectorFactory.create(highConfidenceRegions);
        this.clonalityFactory = clonalityFactory;
    }

    @NotNull
    @Override
    protected Builder enrich(@NotNull final SomaticVariant variant) {
        final Builder builder = super.enrich(variant);

        highConfidenceSelector.select(variant).ifPresent(x -> inHighConfidenceRegion(builder));
        builder.clonality(clonalityFactory.determineClonalityForVariant(variant));

        return builder;
    }

    @NotNull
    private static Builder inHighConfidenceRegion(@NotNull final Builder builder) {
        return builder.highConfidenceRegion(true);
    }
}
