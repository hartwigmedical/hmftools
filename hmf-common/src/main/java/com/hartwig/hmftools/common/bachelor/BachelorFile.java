package com.hartwig.hmftools.common.bachelor;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class BachelorFile {

    private static final Logger LOGGER = LogManager.getLogger(BachelorFile.class);

    private static final String DELIMITER = "\t";

    private BachelorFile() {
    }

    @NotNull
    public static List<GermlineVariant> loadBachelorTsv(@NotNull String filePath) throws IOException {
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
    public static GermlineVariant fromString(@NotNull String line) {
        String[] values = line.split(DELIMITER);

        String program = values[24];
        if (!program.equalsIgnoreCase("HMF")) {
            LOGGER.warn("Unexpected bachelor program found: {}", program);
        }

        String filter = values[2].trim();

        return ImmutableGermlineVariant.builder()
                .passFilter(filter.equalsIgnoreCase("PASS"))
                .gene(values[6].trim())
                .ref(values[4].trim())
                .alt(values[5].trim())
                .codingEffect(CodingEffect.valueOf(values[11].trim()))
                .chromosome(values[0].trim())
                .position(Integer.parseInt(values[1].trim()))
                .hgvsCodingImpact(values[21].trim())
                .hgvsProteinImpact(values[20].trim())
                .totalReadCount(Integer.parseInt(values[16]))
                .alleleReadCount(Integer.parseInt(values[15]))
                .adjustedVAF(Double.parseDouble(values[17].trim()))
                .adjustedCopyNumber(Double.parseDouble(values[18].trim()))
                .biallelic(Boolean.parseBoolean(values[22].trim()))
                .build();
    }
}
