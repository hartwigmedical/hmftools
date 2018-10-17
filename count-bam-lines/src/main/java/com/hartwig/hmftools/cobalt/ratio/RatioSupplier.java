package com.hartwig.hmftools.cobalt.ratio;

import java.io.IOException;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.cobalt.CobaltCount;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFactory;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
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

    @NotNull
    public Multimap<Chromosome, CobaltRatio> generateRatios(@NotNull final Multimap<Chromosome, GCProfile> gcProfiles,
            @NotNull final Multimap<Chromosome, CobaltCount> readCounts) throws IOException {
        LOGGER.info("Applying ratio gc normalization");
        final GCRatioSupplier gcRatioSupplier = new GCRatioSupplier(gcProfiles, readCounts);
        final ListMultimap<Chromosome, ReadRatio> tumorGCRatio = gcRatioSupplier.tumorRatios();
        final ListMultimap<Chromosome, ReadRatio> referenceGCRatio = gcRatioSupplier.referenceRatios();

        LOGGER.info("Applying ratio diploid normalization");
        final ListMultimap<Chromosome, ReadRatio> referenceGCDiploidRatio = new DiploidRatioSupplier(referenceGCRatio).result();

        LOGGER.info("Persisting gc read count medians to {}", outputDirectory);
        final String tumorGCMedianFilename = GCMedianReadCountFile.generateFilename(outputDirectory, tumor);
        final String referenceGCMedianFilename = GCMedianReadCountFile.generateFilename(outputDirectory, reference);
        GCMedianReadCountFile.write(tumorGCMedianFilename, gcRatioSupplier.tumorGCMedianReadCount());
        GCMedianReadCountFile.write(referenceGCMedianFilename, gcRatioSupplier.referenceGCMedianReadCount());

        return CobaltRatioFactory.merge(readCounts, referenceGCRatio, tumorGCRatio, referenceGCDiploidRatio);
    }
}

