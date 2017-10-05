package com.hartwig.hmftools.cobalt.ratio;

import java.util.Collection;
import java.util.Optional;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.cobalt.ImmutableReadCount;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.gc.GCProfile;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class GCRatioSupplier {

    private static final Logger LOGGER = LogManager.getLogger(GCRatioSupplier.class);
    private final GCMedianReadCount tumorGCMedianReadCount;
    private final ListMultimap<String, ReadRatio> tumorRatios;
    private final GCMedianReadCount referenceGCMedianReadCount;
    private final ListMultimap<String, ReadRatio> referenceRatios;

    GCRatioSupplier(final Multimap<String, GCProfile> gcProfiles, final Multimap<String, ReadCount> reference,
            final Multimap<String, ReadCount> tumor) {

        final GenomeRegionSelector<GCProfile> gcProfileSelector = GenomeRegionSelectorFactory.create(gcProfiles);
        final GenomePositionSelector<ReadCount> tumorReadCountSelector = GenomePositionSelectorFactory.create(tumor);

        final GCRatioNormalization tumorRatiosBuilder = new GCRatioNormalization();
        final GCRatioNormalization referenceRatiosBuilder = new GCRatioNormalization();

        for (String chromosomeName : reference.keySet()) {
            if (HumanChromosome.contains(chromosomeName)) {
                final Chromosome chromosome = HumanChromosome.fromString(chromosomeName);
                final Collection<ReadCount> referenceReadCount = reference.get(chromosomeName);
                for (final ReadCount referenceCount : referenceReadCount) {

                    final Optional<GCProfile> optionalGCProfile = gcProfileSelector.select(referenceCount);
                    if (optionalGCProfile.isPresent()) {
                        final GCProfile gcProfile = optionalGCProfile.get();
                        referenceRatiosBuilder.addPosition(chromosome, gcProfile, referenceCount);

                        final ReadCount tumorCount =
                                tumorReadCountSelector.select(referenceCount).orElseGet(() -> empty(referenceCount));
                        tumorRatiosBuilder.addPosition(chromosome, gcProfile, tumorCount);
                    }
                }

            } else {
                LOGGER.info("Excluding unsupported {} chromosome from read ratios", chromosomeName);
            }
        }

        referenceGCMedianReadCount = referenceRatiosBuilder.gcMedianReadCount();
        referenceRatios = referenceRatiosBuilder.build(referenceGCMedianReadCount);

        tumorGCMedianReadCount = tumorRatiosBuilder.gcMedianReadCount();
        tumorRatios = tumorRatiosBuilder.build(tumorGCMedianReadCount);
    }

    ListMultimap<String, ReadRatio> referenceRatios() {
        return referenceRatios;
    }

    GCMedianReadCount referenceGCMedianReadCount() {
        return referenceGCMedianReadCount;
    }

    ListMultimap<String, ReadRatio> tumorRatios() {
        return tumorRatios;
    }

    GCMedianReadCount tumorGCMedianReadCount() {
        return tumorGCMedianReadCount;
    }

    private ReadCount empty(@NotNull final GenomePosition position) {
        return ImmutableReadCount.builder().from(position).readCount(0).build();
    }
}
