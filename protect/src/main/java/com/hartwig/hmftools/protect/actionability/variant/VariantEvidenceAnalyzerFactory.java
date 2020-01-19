package com.hartwig.hmftools.protect.actionability.variant;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.protect.actionability.util.MultiDrugCurator;

import org.jetbrains.annotations.NotNull;

public final class VariantEvidenceAnalyzerFactory {

    private static final String DELIMITER = "\t";

    private VariantEvidenceAnalyzerFactory() {
    }

    @NotNull
    public static VariantEvidenceAnalyzer loadFromFileVariantsAndFileRanges(@NotNull String actionableVariantTsv,
            @NotNull String actionableRangeTsv) throws IOException {
        final List<ActionableVariant> variants = Lists.newArrayList();
        final List<ActionableRange> ranges = Lists.newArrayList();
        final List<String> lineVariants = Files.readAllLines(new File(actionableVariantTsv).toPath());
        final List<String> lineRanges = Files.readAllLines(new File(actionableRangeTsv).toPath());

        // Skip header line for variants
        for (String lineVariant : lineVariants.subList(1, lineVariants.size())) {
            variants.add(fromLineVariants(lineVariant));
        }

        // Skip header line for ranges
        for (String lineRange : lineRanges.subList(1, lineRanges.size())) {
            ranges.add(fromLineRanges(lineRange));
        }

        return new VariantEvidenceAnalyzer(variants, ranges);
    }

    @NotNull
    private static ActionableVariant fromLineVariants(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableVariant.builder()
                .gene(values[0])
                .chromosome(values[1])
                .position(Long.parseLong(values[2]))
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
                .start(Long.parseLong(values[3]))
                .end(Long.parseLong(values[4]))
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
