package com.hartwig.hmftools.cobalt.ratio;

import java.util.Optional;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;

import org.jetbrains.annotations.NotNull;

class GCRatioSupplier {

    private final GCMedianReadCount tumorGCMedianReadCount;
    private final GCMedianReadCount referenceGCMedianReadCount;
    private final ListMultimap<Chromosome, ReadRatio> tumorRatios;
    private final ListMultimap<Chromosome, ReadRatio> referenceRatios;

    GCRatioSupplier(@NotNull final Multimap<Chromosome, GCProfile> gcProfiles, @NotNull final Multimap<Chromosome, CobaltCount> counts) {
        final GenomeRegionSelector<GCProfile> gcProfileSelector = GenomeRegionSelectorFactory.createImproved(gcProfiles);

        final GCRatioNormalization tumorRatiosBuilder = new GCRatioNormalization();
        final GCRatioNormalization referenceRatiosBuilder = new GCRatioNormalization();

        for (Chromosome chromosome : counts.keySet()) {
            for (CobaltCount cobaltPosition : counts.get(chromosome)) {
                final Optional<GCProfile> optionalGCProfile = gcProfileSelector.select(cobaltPosition);
                if (optionalGCProfile.isPresent()) {
                    final GCProfile gcProfile = optionalGCProfile.get();
                    referenceRatiosBuilder.addPosition(chromosome, gcProfile, cobaltPosition.referenceReadCount());
                    tumorRatiosBuilder.addPosition(chromosome, gcProfile, cobaltPosition.tumorReadCount());
                }
            }
        }

        referenceGCMedianReadCount = referenceRatiosBuilder.gcMedianReadCount();
        referenceRatios = referenceRatiosBuilder.build(referenceGCMedianReadCount);

        tumorGCMedianReadCount = tumorRatiosBuilder.gcMedianReadCount();
        tumorRatios = tumorRatiosBuilder.build(tumorGCMedianReadCount);
    }

    @NotNull
    ListMultimap<Chromosome, ReadRatio> referenceRatios() {
        return referenceRatios;
    }

    @NotNull
    GCMedianReadCount referenceGCMedianReadCount() {
        return referenceGCMedianReadCount;
    }

    @NotNull
    ListMultimap<Chromosome, ReadRatio> tumorRatios() {
        return tumorRatios;
    }

    @NotNull
    GCMedianReadCount tumorGCMedianReadCount() {
        return tumorGCMedianReadCount;
    }
}
