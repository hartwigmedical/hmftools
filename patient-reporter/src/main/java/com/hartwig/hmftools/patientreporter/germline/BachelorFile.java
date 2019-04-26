package com.hartwig.hmftools.patientreporter.germline;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class BachelorFile {

    private static final Logger LOGGER = LogManager.getLogger(BachelorFile.class);

    private static final String DELIMITER = ",";
    private static final String BACHELOR_FILE_EXTENSION = "_germline_variants.csv";
    private static final String BACHELOR_HAS_RUN_FILE = "bachelor_found_no_variants";

    private BachelorFile() {
    }

    public static boolean hasBachelorRun(@NotNull String bachelorDirectory, @NotNull String sample) {
        if (findBachelorFilePath(bachelorDirectory, sample) == null) {
            String bachelorHasRunPath = bachelorDirectory + File.separator + BACHELOR_HAS_RUN_FILE;
            return Files.exists(new File(bachelorHasRunPath).toPath());
        }

        return true;
    }

    @Nullable
    public static String findBachelorFilePath(@NotNull String bachelorDirectory, @NotNull String sample) {
        String path = bachelorDirectory + File.separator + sample + BACHELOR_FILE_EXTENSION;
        return Files.exists(new File(path).toPath()) ? path : null;
    }

    @NotNull
    public static List<GermlineVariant> loadBachelorFile(@NotNull String filePath) throws IOException {
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
        // TODO Lynparza can be removed once all samples are on bachelor v1.5+
        if (!program.equalsIgnoreCase("lynparza") && !program.equalsIgnoreCase("hmf")) {
            LOGGER.warn("Unexpected bachelor program found: " + program);
        }

        String filter = values[32].trim();
        int altReadCount = Integer.valueOf(values[16]);
        int totalReadCount = Integer.valueOf(values[17]);

        return ImmutableGermlineVariant.builder()
                .passFilter(filter.equalsIgnoreCase("pass"))
                .gene(values[7].trim())
                .hgvsCodingImpact(values[26].trim())
                .hgvsProteinImpact(values[25].trim())
                .totalReadCount(totalReadCount)
                .alleleReadCount(altReadCount)
                .germlineStatus(values[30].trim())
                .adjustedVAF(Double.valueOf(values[19]))
                .adjustedCopyNumber(Double.valueOf(values[18]))
                .minorAllelePloidy(Double.valueOf(values[31]))
                .biallelic(Boolean.valueOf(values[27].trim()))
                .build();
    }
}
