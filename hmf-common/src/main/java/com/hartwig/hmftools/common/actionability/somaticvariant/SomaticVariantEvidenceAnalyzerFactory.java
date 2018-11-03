package com.hartwig.hmftools.common.actionability.somaticvariant;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.actionability.util.MultiDrugCurator;

import org.jetbrains.annotations.NotNull;

public final class SomaticVariantEvidenceAnalyzerFactory {

    private static final String DELIMITER = "\t";

    private SomaticVariantEvidenceAnalyzerFactory() {
    }

    @NotNull
    public static SomaticVariantEvidenceAnalyzer loadFromFileVariantsAndFileRanges(@NotNull String fileVariants, @NotNull String fileRanges)
            throws IOException {
        final List<ActionableSomaticVariant> variants = new ArrayList<>();
        final List<ActionableRange> ranges = new ArrayList<>();
        final List<String> lineVariants = Files.readAllLines(new File(fileVariants).toPath());
        final List<String> lineRanges = Files.readAllLines(new File(fileRanges).toPath());

        // KODU: Skip header line
        for (String lineVariant : lineVariants.subList(1, lineVariants.size())) {
            variants.add(fromLineVariants(lineVariant));
        }

        // KODU: Skip header line
        for (String lineRange : lineRanges.subList(1, lineRanges.size())) {
            ranges.add(fromLineRanges(lineRange));
        }

        return new SomaticVariantEvidenceAnalyzer(variants, ranges);
    }

    @NotNull
    private static ActionableSomaticVariant fromLineVariants(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableSomaticVariant.builder()
                .gene(values[0])
                .chromosome(values[1])
                .position(Long.valueOf(values[2]))
                .ref(values[3])
                .alt(values[4])
                .source(values[5])
                .reference(values[6])
                .drug(MultiDrugCurator.reformat(values[7]))
                .drugsType(values[8])
                .cancerType(values[9])
                .level(values[11])
                .response(values[14])
                .build();
    }

    @NotNull
    private static ActionableRange fromLineRanges(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableRange.builder()
                .gene(values[0])
                .chromosome(values[2])
                .start(Long.valueOf(values[3]))
                .end(Long.valueOf(values[4]))
                .source(values[6])
                .reference(values[7])
                .drug(MultiDrugCurator.reformat(values[8]))
                .drugsType(values[9])
                .cancerType(values[10])
                .level(values[12])
                .response(values[15])
                .build();
    }
}
