package com.hartwig.hmftools.cobalt.ratio;

import java.util.Collection;
import java.util.Optional;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.gc.GCProfile;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

class GCRatioSupplier {


    private static final Logger LOGGER = LogManager.getLogger(GCRatioSupplier.class);
    private final ListMultimap<String, ReadRatio> ratios;
    private final GCMedianReadCount gcMedianReadCount;

    GCRatioSupplier(final Multimap<String, GCProfile> gcProfiles, final Multimap<String, ReadCount> readCounts) {

        final GenomeRegionSelector<GCProfile> gcProfileSelector = GenomeRegionSelectorFactory.create(gcProfiles);
        final GCRatioNormalization ratiosBuilder = new GCRatioNormalization();

        for (String chromosomeName : readCounts.keySet()) {
            if (HumanChromosome.contains(chromosomeName)) {
                final Chromosome chromosome = HumanChromosome.fromString(chromosomeName);
                final Collection<ReadCount> chromosomeReadCounts = readCounts.get(chromosomeName);
                for (final ReadCount readCount : chromosomeReadCounts) {
                    final Optional<GCProfile> optionalGCProfile = gcProfileSelector.select(readCount);
                    if (optionalGCProfile.isPresent()) {
                        final GCProfile gcProfile = optionalGCProfile.get();
                        ratiosBuilder.addPosition(chromosome, gcProfile, readCount);
                    }
                }
            } else {
                LOGGER.info("Excluding unsupported {} chromosome from read ratios", chromosomeName);
            }
        }

        gcMedianReadCount = ratiosBuilder.gcMedianReadCount();
        ratios = ratiosBuilder.build(gcMedianReadCount);
    }

    ListMultimap<String, ReadRatio> ratios() {
        return ratios;
    }

    GCMedianReadCount gcMedianReadCount() {
        return gcMedianReadCount;
    }
}
