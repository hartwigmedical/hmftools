package com.hartwig.hmftools.purple.config;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.amber.qc.AmberQCFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.purple.segment.ExpectedBAF;
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
public interface AmberData {

    Logger LOGGER = LogManager.getLogger(AmberData.class);
    int DEFAULT_READ_DEPTH = 100;

    @NotNull
    Gender gender();

    @NotNull
    Multimap<Chromosome, AmberBAF> bafs();

    @NotNull
    Multimap<Chromosome, PCFPosition> tumorSegments();

    int averageTumorDepth();

    double contamination();

    @NotNull
    static AmberData createAmberData(@NotNull final CommonConfig commonConfig)
            throws ParseException, IOException {
        final String amberDirectory = commonConfig.amberDirectory();
        final String amberFilename = AmberBAFFile.generateAmberFilenameForReading(amberDirectory, commonConfig.tumorSample());
        if (!new File(amberFilename).exists()) {
            throw new ParseException("Unable to open amber baf file: " + amberFilename);
        }

        final String pcfFilename = PCFFile.generateBAFFilename(amberDirectory, commonConfig.tumorSample());
        if (!new File(pcfFilename).exists()) {
            throw new ParseException("Unable to open amber pcf file: " + pcfFilename);
        }

        final String qcFile = AmberQCFile.generateFilename(amberDirectory, commonConfig.tumorSample());
        if (!new File(qcFile).exists()) {
            throw new ParseException("Unable to open amber qc file: " + qcFile);
        }

        LOGGER.info("Reading amber QC from {}", qcFile);
        final double contamination = AmberQCFile.read(qcFile).contamination();

        LOGGER.info("Reading amber bafs from {}", amberFilename);
        final Multimap<Chromosome, AmberBAF> bafs = AmberBAFFile.read(amberFilename);

        LOGGER.info("Reading amber pcfs from {}", pcfFilename);
        final Multimap<Chromosome, PCFPosition> tumorSegments =
                PCFFile.readPositions(commonConfig.windowSize(), PCFSource.TUMOR_BAF, pcfFilename);

        int averageTumorDepth = (int) Math.round(bafs.values()
                .stream()
                .mapToInt(AmberBAF::tumorDepth)
                .filter(x -> x > 0)
                .average()
                .orElse(DEFAULT_READ_DEPTH));
        LOGGER.info("Average amber tumor depth is {} reads implying an ambiguous BAF of {}",
                averageTumorDepth,
                new DecimalFormat("0.000").format(ExpectedBAF.expectedBAF(averageTumorDepth)));

        final Gender gender = Gender.fromAmber(bafs);

        return ImmutableAmberData.builder()
                .contamination(contamination)
                .averageTumorDepth(averageTumorDepth)
                .bafs(bafs)
                .tumorSegments(tumorSegments)
                .gender(gender)
                .build();
    }
}
