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

public class RatioSupplier {

    private static final Logger LOGGER = LogManager.getLogger(GCRatioSupplier.class);

    private final String tumor;
    private final String reference;
    private final String outputDirectory;

    public RatioSupplier(final String reference, final String tumor, final String outputDirectory) {
        this.tumor = tumor;
        this.reference = reference;
        this.outputDirectory = outputDirectory;
    }

    public void generateRatios(
            @NotNull final Multimap<String, GCProfile> gcProfiles,
            @NotNull final Multimap<String, ReadCount> referenceCounts,
            @NotNull final Multimap<String, ReadCount> tumorCounts)
            throws IOException {

        LOGGER.info("Applying ratio gc normalization");
        final GCRatioSupplier gcRatioSupplier = new GCRatioSupplier(gcProfiles, referenceCounts, tumorCounts);
        final ListMultimap<String, ReadRatio> gcTumorRatios = gcRatioSupplier.tumorRatios();
        final ListMultimap<String, ReadRatio> gcReferenceRatios = gcRatioSupplier.referenceRatios();

        LOGGER.info("Applying ratio diploid normalization");
        final ListMultimap<String, ReadRatio> referenceRatios = new DiploidRatioSupplier(gcReferenceRatios).result();

        LOGGER.info("Persisting ratios {}", outputDirectory);
        ReadRatioFile.write(ReadRatioFile.generateFilename(outputDirectory, tumor), gcTumorRatios);
        ReadRatioFile.write(ReadRatioFile.generateFilename(outputDirectory, reference), referenceRatios);
        ReadRatioFile.write(ReadRatioFile.generateFilename(outputDirectory, reference) + ".gc", gcReferenceRatios);

        LOGGER.info("Persisting gc read count medians to {}", outputDirectory);
        final String tumorGCMedianFilename = GCMedianReadCountFile.generateFilename(outputDirectory, tumor);
        final String referenceGCMedianFilename = GCMedianReadCountFile.generateFilename(outputDirectory, reference);
        GCMedianReadCountFile.write(tumorGCMedianFilename, gcRatioSupplier.tumorGCMedianReadCount());
        GCMedianReadCountFile.write(referenceGCMedianFilename, gcRatioSupplier.referenceGCMedianReadCount());
    }
}

