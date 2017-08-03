package com.hartwig.hmftools.purple.ratio;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.position.GenomePositionSelector;
import com.hartwig.hmftools.common.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.purple.ratio.GCContent;
import com.hartwig.hmftools.common.purple.ratio.GCMedianFile;
import com.hartwig.hmftools.common.purple.ratio.NormalizedRatios;
import com.hartwig.hmftools.common.purple.ratio.NormalizedRatiosBuilder;
import com.hartwig.hmftools.common.purple.ratio.ReadRatio;
import com.hartwig.hmftools.common.purple.ratio.ReadRatioFile;
import com.hartwig.hmftools.purple.config.CommonConfig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ReadCountRatioSupplier implements RatioSupplier {

    private static final Logger LOGGER = LogManager.getLogger(ReadCountRatioSupplier.class);

    private final Multimap<String, ReadRatio> referenceRatios;
    private final Multimap<String, ReadRatio> tumorRatios;

    public ReadCountRatioSupplier(final CommonConfig config, final Multimap<String, GCContent> gcContent)
            throws IOException, HartwigException {

        final String tumorRatioFile = ReadRatioFile.generateFilename(config.outputDirectory(), config.tumorSample());
        final String referenceRatioFile = ReadRatioFile.generateFilename(config.outputDirectory(), config.refSample());

        if (!config.forceRecalculate() && new File(tumorRatioFile).exists() && new File(referenceRatioFile).exists()) {
            LOGGER.info("Loading reference ratios from {}", referenceRatioFile);
            referenceRatios = ReadRatioFile.read(referenceRatioFile);

            LOGGER.info("Loading tumor ratios from {}", tumorRatioFile);
            tumorRatios = ReadRatioFile.read(tumorRatioFile);
        } else {

            final ReadCountSupplier readCountSupplier = new ReadCountSupplier(config);
            final Multimap<String, ReadCount> tumorReadCount = readCountSupplier.tumorReadCount();
            final Multimap<String, ReadCount> normalReadCount = readCountSupplier.referenceReadCount();

            final GenomePositionSelector<ReadCount> referenceReadCountSelector = GenomePositionSelectorFactory.create(normalReadCount);
            final GenomePositionSelector<ReadCount> tumorReadCountSelector = GenomePositionSelectorFactory.create(tumorReadCount);

            LOGGER.info("Generating gc normalized read ratios");
            final NormalizedRatiosBuilder normalRatiosBuilder = new NormalizedRatiosBuilder();
            final NormalizedRatiosBuilder tumorRatiosBuilder = new NormalizedRatiosBuilder();
            for (String chromosomeName : normalReadCount.keySet()) {
                final Chromosome chromosome = HumanChromosome.fromString(chromosomeName);

                List<GCContent> chromosomeGCContent = (List<GCContent>) gcContent.get(chromosomeName);
                for (GCContent windowGCContent : chromosomeGCContent) {

                    referenceReadCountSelector.select(windowGCContent)
                            .ifPresent(x -> normalRatiosBuilder.addPosition(chromosome, windowGCContent, x));
                    tumorReadCountSelector.select(windowGCContent)
                            .ifPresent(x -> tumorRatiosBuilder.addPosition(chromosome, windowGCContent, x));
                }
            }

            NormalizedRatios normalizedReferenceRatios = normalRatiosBuilder.build();
            NormalizedRatios normalizedTumorRatios = tumorRatiosBuilder.build();

            referenceRatios = normalizedReferenceRatios.normalisedRatios();
            tumorRatios = normalizedTumorRatios.normalisedRatios();

            final String tumorGCMedianFileName = GCMedianFile.generateFilename(config.outputDirectory(), config.tumorSample());
            LOGGER.info("Persisting read count medians to {}", tumorGCMedianFileName);
            GCMedianFile.write(tumorGCMedianFileName, normalizedTumorRatios.medianReadCount());

            final String referenceGCMedianFileName = GCMedianFile.generateFilename(config.outputDirectory(), config.refSample());
            LOGGER.info("Persisting read count medians to {}", referenceGCMedianFileName);
            GCMedianFile.write(referenceGCMedianFileName, normalizedReferenceRatios.medianReadCount());

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

}
