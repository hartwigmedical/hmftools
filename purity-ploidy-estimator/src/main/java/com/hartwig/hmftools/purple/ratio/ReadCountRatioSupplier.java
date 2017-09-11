package com.hartwig.hmftools.purple.ratio;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Optional;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.cobalt.ImmutableReadCount;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gc.GCMedianReadCountFile;
import com.hartwig.hmftools.common.gc.GCProfile;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.purple.ratio.NormalizedRatios;
import com.hartwig.hmftools.common.purple.ratio.NormalizedRatiosBuilder;
import com.hartwig.hmftools.common.purple.ratio.ReadRatio;
import com.hartwig.hmftools.common.purple.ratio.ReadRatioFile;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ReadCountRatioSupplier implements RatioSupplier {

    private static final Logger LOGGER = LogManager.getLogger(ReadCountRatioSupplier.class);

    private final Multimap<String, ReadRatio> referenceRatios;
    private final Multimap<String, ReadRatio> tumorRatios;

    public ReadCountRatioSupplier(final CommonConfig config, final Multimap<String, GCProfile> gcContent)
            throws IOException, HartwigException {

        final String tumorRatioFile = ReadRatioFile.generateFilename(config.outputDirectory(), config.tumorSample());
        final String referenceRatioFile = ReadRatioFile.generateFilename(config.outputDirectory(), config.refSample());

        if (!config.forceSegmentation() && new File(tumorRatioFile).exists() && new File(referenceRatioFile).exists()) {
            LOGGER.info("Loading reference ratios from {}", referenceRatioFile);
            referenceRatios = ReadRatioFile.read(referenceRatioFile);

            LOGGER.info("Loading tumor ratios from {}", tumorRatioFile);
            tumorRatios = ReadRatioFile.read(tumorRatioFile);
        } else {

            final ReadCountSupplier readCountSupplier = new ReadCountSupplier(config);
            final Multimap<String, ReadCount> tumorReadCount = readCountSupplier.tumorReadCount();
            final Multimap<String, ReadCount> normalReadCount = readCountSupplier.referenceReadCount();

            final GenomeRegionSelector<GCProfile> gcProfileSelector = GenomeRegionSelectorFactory.create(gcContent);
            final GenomePositionSelector<ReadCount> tumorReadCountSelector = GenomePositionSelectorFactory.create(tumorReadCount);

            LOGGER.info("Generating gc normalized read ratios");
            final NormalizedRatiosBuilder normalRatiosBuilder = new NormalizedRatiosBuilder();
            final NormalizedRatiosBuilder tumorRatiosBuilder = new NormalizedRatiosBuilder();
            for (String chromosomeName : normalReadCount.keySet()) {
                if (HumanChromosome.contains(chromosomeName)) {
                    final Chromosome chromosome = HumanChromosome.fromString(chromosomeName);
                    final Collection<ReadCount> referenceReadCount = normalReadCount.get(chromosomeName);
                    for (final ReadCount referenceCount : referenceReadCount) {

                        final Optional<GCProfile> optionalGCProfile = gcProfileSelector.select(referenceCount);
                        if (optionalGCProfile.isPresent()) {
                            final GCProfile gcProfile = optionalGCProfile.get();
                            normalRatiosBuilder.addPosition(chromosome, gcProfile, referenceCount);

                            final ReadCount tumorCount =
                                    tumorReadCountSelector.select(referenceCount).orElseGet(() -> empty(referenceCount));
                            tumorRatiosBuilder.addPosition(chromosome, gcProfile, tumorCount);
                        }
                    }

                } else {
                    LOGGER.info("Excluding unsupported {} chromosome from read ratios", chromosomeName);
                }
            }

            NormalizedRatios normalizedReferenceRatios = normalRatiosBuilder.build();
            NormalizedRatios normalizedTumorRatios = tumorRatiosBuilder.build();

            referenceRatios = normalizedReferenceRatios.normalisedRatios();
            tumorRatios = normalizedTumorRatios.normalisedRatios();

            final String tumorGCMedianFileName = GCMedianReadCountFile.generateFilename(config.outputDirectory(), config.tumorSample());
            LOGGER.info("Persisting read count medians to {}", tumorGCMedianFileName);
            GCMedianReadCountFile.write(tumorGCMedianFileName, normalizedTumorRatios.medianReadCount());

            final String referenceGCMedianFileName = GCMedianReadCountFile.generateFilename(config.outputDirectory(), config.refSample());
            LOGGER.info("Persisting read count medians to {}", referenceGCMedianFileName);
            GCMedianReadCountFile.write(referenceGCMedianFileName, normalizedReferenceRatios.medianReadCount());

            LOGGER.info("Persisting gc normalized read ratios to file");
            ReadRatioFile.write(config.outputDirectory(), config.refSample(), referenceRatios);
            ReadRatioFile.write(config.outputDirectory(), config.tumorSample(), tumorRatios);
        }
    }

    @Override
    @NotNull
    public Multimap<String, ReadRatio> tumorRatios() {
        return tumorRatios;
    }

    @Override
    @NotNull
    public Multimap<String, ReadRatio> referenceRatios() {
        return referenceRatios;
    }

    private ReadCount empty(@NotNull final GenomePosition position) {
        return ImmutableReadCount.builder().from(position).readCount(0).build();
    }

}
