package com.hartwig.hmftools.patientreporter.actionability;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class ActionabilityAnalyzer {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityAnalyzer.class);
    static final String DELIMITER = "\t";

    @NotNull
    public static ActionabilityAnalyzer loadFromFile(String file) throws IOException {
        final String line =  Files.readAllLines(new File(file).toPath()).get(1);
        fromLine(line);
        LOGGER.info(fromLine(line));
        return null;
    }

    @NotNull
    static ActionableVariant fromLine(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableVariant.builder().gene(values[0]).chromosome(values[1]).build();
    }

    public boolean isActionable(@NotNull SomaticVariant variant) {
        String variantGene = variant.gene();
        String variantChromosome = variant.chromosome();
        String variantPosition = variant.chromosomePosition();
        String variantAlt= variant.alt();
        String variantRef = variant.ref();
        LOGGER.info(variantGene + variantChromosome + variantPosition + variantAlt + variantRef);
        // TODO (LIST): Implement
        return false;
    }
}
