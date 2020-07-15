package com.hartwig.hmftools.patientreporter.variants.germline;

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

    private static final String PROGRAM_COLUMN = "program";
    private static final String GENE_COLUMN = "gene";
    private static final String CHROMOSOME_COLUMN = "chromosome";
    private static final String POSITION_COLUMN = "position";
    private static final String REF_COLUMN = "ref";
    private static final String ALT_COLUMN = "alts";
    private static final String FILTER_COLUMN = "filter";
    private static final String CODING_EFFECT_COLUMN = "codingEffect";
    private static final String HGVS_CODING_COLUMN = "hgvsCoding";
    private static final String HGVS_PROTEIN_COLUMN = "hgvsProtein";
    private static final String ALLELE_READ_COLUMN = "alleleReadCount";
    private static final String TOTAL_READ_COLUMN = "totalReadCount";
    private static final String ADJUSTED_VAF_COLUMN = "adjustedVaf";
    private static final String ADJUSTED_COPY_NUMBER_COLUMN = "adjustedCopyNumber";
    private static final String BIALLELIC_COLUMN = "biallelic";

    private BachelorFile() {
    }

    @NotNull
    public static List<GermlineVariant> loadBachelorTsv(@NotNull String filePath) throws IOException {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    @NotNull
    private static List<GermlineVariant> fromLines(@NotNull List<String> lines) {
        List<GermlineVariant> germlineVariants = Lists.newArrayList();
        String[] headers = lines.get(0).split(DELIMITER);
        // Skip header line
        for (String line : lines.subList(1, lines.size())) {
            germlineVariants.add(fromString(headers, line));
        }
        return germlineVariants;
    }

    @NotNull
    private static GermlineVariant fromString(@NotNull String[] headers, @NotNull String line) {
        String[] values = line.split(DELIMITER);

        String program = values[findIndexForValue(headers, PROGRAM_COLUMN)];
        if (!program.equalsIgnoreCase("HMF")) {
            LOGGER.warn("Unexpected bachelor program found: {}", program);
        }

        String filter = values[findIndexForValue(headers, FILTER_COLUMN)].trim();

        return ImmutableGermlineVariant.builder()
                .gene(values[findIndexForValue(headers, GENE_COLUMN)].trim())
                .chromosome(values[findIndexForValue(headers, CHROMOSOME_COLUMN)].trim())
                .position(Integer.parseInt(values[findIndexForValue(headers, POSITION_COLUMN)].trim()))
                .ref(values[findIndexForValue(headers, REF_COLUMN)].trim())
                .alt(values[findIndexForValue(headers, ALT_COLUMN)].trim())
                .passFilter(filter.equalsIgnoreCase("PASS"))
                .codingEffect(CodingEffect.valueOf(values[findIndexForValue(headers, CODING_EFFECT_COLUMN)].trim()))
                .hgvsCodingImpact(values[findIndexForValue(headers, HGVS_CODING_COLUMN)].trim())
                .hgvsProteinImpact(values[findIndexForValue(headers, HGVS_PROTEIN_COLUMN)].trim())
                .alleleReadCount(Integer.parseInt(values[findIndexForValue(headers, ALLELE_READ_COLUMN)]))
                .totalReadCount(Integer.parseInt(values[findIndexForValue(headers, TOTAL_READ_COLUMN)]))
                .adjustedVAF(Double.parseDouble(values[findIndexForValue(headers, ADJUSTED_VAF_COLUMN)].trim()))
                .adjustedCopyNumber(Double.parseDouble(values[findIndexForValue(headers, ADJUSTED_COPY_NUMBER_COLUMN)].trim()))
                .biallelic(Boolean.parseBoolean(values[findIndexForValue(headers, BIALLELIC_COLUMN)].trim()))
                .build();
    }

    private static int findIndexForValue(@NotNull String[] headers, @NotNull String value) {
        for (int i = 0; i < headers.length; i++) {
            if (headers[i].equals(value)) {
                return i;
            }
        }

        throw new IllegalStateException(String.format("Value %s not found in headers of BachelorFile!", value));
    }
}
