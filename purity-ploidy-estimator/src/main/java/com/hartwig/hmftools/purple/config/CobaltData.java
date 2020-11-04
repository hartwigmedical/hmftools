package com.hartwig.hmftools.purple.config;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatioFactory;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.utils.pcf.PCFFile;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;

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
    CobaltChromosomes cobaltChromosomes();

    @NotNull
    ListMultimap<Chromosome, CobaltRatio> ratios();

    @NotNull
    Multimap<Chromosome, PCFPosition> tumorSegments();

    @NotNull
    Multimap<Chromosome, PCFPosition> referenceSegments();

    @NotNull
    default Gender gender() {
        return cobaltChromosomes().gender();
    }

    @NotNull
    static CobaltData createCobaltData(@NotNull final CommonConfig commonConfig, @NotNull final Gender amberGender)
            throws ParseException, IOException {
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
        final ListMultimap<Chromosome, CobaltRatio> ratios = commonConfig.tumorOnly()
                ? CobaltRatioFile.readTumorOnly(cobaltFilename, amberGender)
                : CobaltRatioFile.read(cobaltFilename);

        LOGGER.info("Reading cobalt reference segments from {}", referenceSegmentFile);
        final Multimap<Chromosome, PCFPosition> referenceSegments =
                PCFFile.readPositions(commonConfig.windowSize(), PCFSource.REFERENCE_RATIO, referenceSegmentFile);

        LOGGER.info("Reading cobalt tumor segments from {}", tumorSegmentFile);
        final Multimap<Chromosome, PCFPosition> tumorSegments =
                PCFFile.readPositions(commonConfig.windowSize(), PCFSource.TUMOR_RATIO, tumorSegmentFile);

        final List<MedianRatio> medianRatios = MedianRatioFactory.create(ratios);
        final CobaltChromosomes cobaltChromosomes = new CobaltChromosomes(medianRatios);

        return ImmutableCobaltData.builder()
                .ratios(ratios)
                .cobaltChromosomes(cobaltChromosomes)
                .tumorSegments(tumorSegments)
                .referenceSegments(referenceSegments)
                .build();
    }
}
