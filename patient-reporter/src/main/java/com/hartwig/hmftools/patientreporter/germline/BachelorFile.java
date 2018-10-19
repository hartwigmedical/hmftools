package com.hartwig.hmftools.patientreporter.germline;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class BachelorFile {

    private static final Logger LOGGER = LogManager.getLogger(BachelorFile.class);
    private static final String DELIMITER = ",";

    private BachelorFile() {
    }

    @NotNull
    public static List<GermlineVariant> loadBachelorFile(@NotNull String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    @NotNull
    private static List<GermlineVariant> fromLines(@NotNull List<String> lines) {
        List<GermlineVariant> germlineVariants = Lists.newArrayList();
        // KODU: Skip header line
        for (String line : lines.subList(1, lines.size())) {
            germlineVariants.add(fromString(line));
        }
        return germlineVariants;
    }

    @NotNull
    private static GermlineVariant fromString(@NotNull final String germlineVariant) {
        String[] values = germlineVariant.split(DELIMITER);

        String program = values[1];
        if (!program.equalsIgnoreCase("lynparza")) {
            LOGGER.warn("Unexpected bachelor program found: " + program);
        }

        return ImmutableGermlineVariant.builder()
                .gene(values[8].trim())
                .hgvsCodingImpact(values[24].trim())
                .hgvsProteinImpact(values[23].trim())
                .totalReadCount(Integer.valueOf(values[15]))
                .alleleReadCount(Integer.valueOf(values[14]))
                .germlineStatus(Strings.EMPTY)
                .adjustedVAF(Double.valueOf(values[17]))
                .adjustedCopyNumber(Double.valueOf(values[16]))
                .minorAllelePloidy(Double.valueOf(values[29]))
                .biallelic(Boolean.valueOf(values[25].trim()))
                .build();
    }
}
