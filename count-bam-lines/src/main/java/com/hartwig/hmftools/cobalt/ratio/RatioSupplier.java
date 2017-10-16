package com.hartwig.hmftools.cobalt.ratio;

import java.io.IOException;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltPositionFactory;
import com.hartwig.hmftools.common.cobalt.CobaltPositionFile;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatioFile;
import com.hartwig.hmftools.common.gc.GCMedianReadCountFile;
import com.hartwig.hmftools.common.gc.GCProfile;
import com.hartwig.hmftools.common.purple.gender.Gender;

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

    public void generateRatios(@NotNull final Multimap<String, GCProfile> gcProfiles,
            @NotNull final Multimap<String, ReadCount> referenceCounts, @NotNull final Multimap<String, ReadCount> tumorCounts)
            throws IOException {

        final CobaltPositionFactory cobaltPositionFactory = new CobaltPositionFactory();
        cobaltPositionFactory.addReadCounts(referenceCounts, tumorCounts);

        LOGGER.info("Applying ratio gc normalization");
        final GCRatioSupplier gcRatioSupplier = new GCRatioSupplier(gcProfiles, cobaltPositionFactory.build());
        final ListMultimap<String, ReadRatio> tumorGCRatio = gcRatioSupplier.tumorRatios();
        final ListMultimap<String, ReadRatio> referenceGCRatio = gcRatioSupplier.referenceRatios();

        LOGGER.info("Applying ratio diploid normalization");
        final Gender gender = Gender.fromReferenceReadRatios(referenceGCRatio);
        final ListMultimap<String, ReadRatio> referenceGCDiploidRatio = new DiploidRatioSupplier(gender, referenceGCRatio).result();

        cobaltPositionFactory.addRatios(referenceGCRatio, tumorGCRatio, referenceGCDiploidRatio);
        final String outputFile = CobaltPositionFile.generateFilename(outputDirectory, tumor);
        LOGGER.info("Persisting ratios to {}", outputFile);
        CobaltPositionFile.write(outputFile, cobaltPositionFactory.build());

        LOGGER.info("Persisting ratios {}", outputDirectory);
        ReadRatioFile.write(ReadRatioFile.generateFilename(outputDirectory, tumor), tumorGCRatio);
        ReadRatioFile.write(ReadRatioFile.generateFilename(outputDirectory, reference), referenceGCDiploidRatio);
        ReadRatioFile.write(ReadRatioFile.generateFilename(outputDirectory, reference) + ".gc", referenceGCRatio);


        LOGGER.info("Persisting gc read count medians to {}", outputDirectory);
        final String tumorGCMedianFilename = GCMedianReadCountFile.generateFilename(outputDirectory, tumor);
        final String referenceGCMedianFilename = GCMedianReadCountFile.generateFilename(outputDirectory, reference);
        GCMedianReadCountFile.write(tumorGCMedianFilename, gcRatioSupplier.tumorGCMedianReadCount());
        GCMedianReadCountFile.write(referenceGCMedianFilename, gcRatioSupplier.referenceGCMedianReadCount());
    }
}

