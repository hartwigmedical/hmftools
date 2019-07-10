package com.hartwig.hmftools.purple.config;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.pcf.PCFFile;
import com.hartwig.hmftools.common.pcf.PCFPosition;
import com.hartwig.hmftools.common.pcf.PCFSource;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface CobaltData {

    Logger LOGGER = LogManager.getLogger(CobaltData.class);

    @NotNull
    Gender gender();

    @NotNull
    ListMultimap<Chromosome, CobaltRatio> ratios();

    @NotNull
    Multimap<Chromosome, PCFPosition> tumorSegments();

    @NotNull
    Multimap<Chromosome, PCFPosition> referenceSegments();

    @NotNull
    static CobaltData createCobaltData(@NotNull final CommonConfig commonConfig) throws ParseException, IOException {
        final String cobaltDirectory = commonConfig.cobaltDirectory();
        final String cobaltFilename = CobaltRatioFile.generateFilenameForReading(cobaltDirectory, commonConfig.tumorSample());
        if (!new File(cobaltFilename).exists()) {
            throw new ParseException("Unable to open cobalt ratio file: " + cobaltFilename);
        }

        final String referenceSegmentFile = PCFFile.generateRatioFilename(cobaltDirectory, commonConfig.refSample());
        if (!new File(referenceSegmentFile).exists()) {
            throw new ParseException("Unable to open cobalt reference pcf file: " + referenceSegmentFile);
        }

        final String tumorSegmentFile = PCFFile.generateRatioFilename(cobaltDirectory, commonConfig.tumorSample());
        if (!new File(tumorSegmentFile).exists()) {
            throw new ParseException("Unable to open cobalt tumor pcf file: " + tumorSegmentFile);
        }

        LOGGER.info("Reading cobalt ratios from {}", cobaltFilename);
        final ListMultimap<Chromosome, CobaltRatio> ratios = CobaltRatioFile.read(cobaltFilename);
        final Gender gender = Gender.fromCobalt(ratios);

        LOGGER.info("Reading cobalt reference segments from {}", referenceSegmentFile);
        final Multimap<Chromosome, PCFPosition> referenceSegments =
                PCFFile.readPositions(commonConfig.windowSize(), PCFSource.REFERENCE_RATIO, referenceSegmentFile);

        LOGGER.info("Reading cobalt tumor segments from {}", tumorSegmentFile);
        final Multimap<Chromosome, PCFPosition> tumorSegments =
                PCFFile.readPositions(commonConfig.windowSize(), PCFSource.TUMOR_RATIO, tumorSegmentFile);

        return ImmutableCobaltData.builder()
                .ratios(ratios)
                .gender(gender)
                .tumorSegments(tumorSegments)
                .referenceSegments(referenceSegments)
                .build();
    }

}
