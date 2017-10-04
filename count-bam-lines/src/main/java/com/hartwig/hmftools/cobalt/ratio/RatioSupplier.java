package com.hartwig.hmftools.cobalt.ratio;

import java.io.IOException;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatioFile;
import com.hartwig.hmftools.common.gc.GCMedianReadCountFile;
import com.hartwig.hmftools.common.gc.GCProfile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RatioSupplier  {

    private static final Logger LOGGER = LogManager.getLogger(GCRatioSupplier.class);

    private final String sample;
    private final String outputDirectory;
    private final boolean reference;

    public RatioSupplier(@NotNull final String sample, @NotNull final String outputDirectory, final boolean reference) {
        this.sample = sample;
        this.outputDirectory = outputDirectory;
        this.reference = reference;
    }

    public void generateRatios(@NotNull final Multimap<String, GCProfile> gcProfiles, @NotNull final Multimap<String, ReadCount> readCounts)
            throws IOException {

        LOGGER.info("Generating gc normalized read ratios");
        final GCRatioSupplier gcRatioSupplier = new GCRatioSupplier(gcProfiles, readCounts);

        final Multimap<String, ReadRatio> ratios;
        if (reference) {
            LOGGER.info("Applying diploid normalization");
            ListMultimap<String, ReadRatio> intermediateRatios = gcRatioSupplier.ratios();
            ratios = new DiploidRatioSupplier(intermediateRatios).result();

            final String ratioFileName = ReadRatioFile.generateFilename(outputDirectory, sample) + ".gc";
            LOGGER.info("Persisting gc corrected read ratios to {}", ratioFileName);
            ReadRatioFile.write(ratioFileName, intermediateRatios);

        } else {
            ratios = gcRatioSupplier.ratios();
        }

        final String ratioFileName = ReadRatioFile.generateFilename(outputDirectory, sample);
        LOGGER.info("Persisting final read ratios to {}", ratioFileName);
        ReadRatioFile.write(ratioFileName, ratios);

        final String gcMedianFileName = GCMedianReadCountFile.generateFilename(outputDirectory, sample);
        LOGGER.info("Persisting read count medians to {}", gcMedianFileName);
        GCMedianReadCountFile.write(gcMedianFileName, gcRatioSupplier.gcMedianReadCount());
    }
}

