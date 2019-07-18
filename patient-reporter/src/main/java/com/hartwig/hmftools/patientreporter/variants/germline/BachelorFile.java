package com.hartwig.hmftools.patientreporter.variants.germline;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class BachelorFile {

    private static final Logger LOGGER = LogManager.getLogger(BachelorFile.class);

    private static final String DELIMITER = ",";

    private BachelorFile() {
    }

    @NotNull
    public static List<GermlineVariant> loadBachelorCsv(@NotNull String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    @NotNull
    private static List<GermlineVariant> fromLines(@NotNull List<String> lines) {
        List<GermlineVariant> germlineVariants = Lists.newArrayList();
        // Skip header line
        for (String line : lines.subList(1, lines.size())) {
            germlineVariants.add(fromString(line));
        }
        return germlineVariants;
    }

    @NotNull
    private static GermlineVariant fromString(@NotNull String line) {
        String[] values = line.split(DELIMITER);

        String program = values[1];
        if (!program.equalsIgnoreCase("hmf")) {
            LOGGER.warn("Unexpected bachelor program found: " + program);
        }

        String filter = values[32].trim();
        int altReadCount = Integer.valueOf(values[16]);
        int totalReadCount = Integer.valueOf(values[17]);

        return ImmutableGermlineVariant.builder()
                .passFilter(filter.equalsIgnoreCase("pass"))
                .gene(values[7].trim())
                .chromosome(values[2].trim())
                .position(Integer.valueOf(values[3]))
                .hgvsCodingImpact(values[26].trim())
                .hgvsProteinImpact(values[25].trim())
                .totalReadCount(totalReadCount)
                .alleleReadCount(altReadCount)
                .adjustedVAF(Double.valueOf(values[19]))
                .adjustedCopyNumber(Double.valueOf(values[18]))
                .biallelic(Boolean.valueOf(values[27].trim()))
                .build();
    }
}
