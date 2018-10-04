package com.hartwig.hmftools.actionability.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.actionability.cancerTypeMapping.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ActionabilityVariantsAnalyzer {
    private static final Logger LOGGER = LogManager.getLogger(ActionabilityVariantsAnalyzer.class);
    private static final String DELIMITER = "\t";

    @NotNull
    private final List<ActionabilityVariant> variants;
    private final List<ActionabilityRange> variantsRanges;

    ActionabilityVariantsAnalyzer(@NotNull final List<ActionabilityVariant> variants, @NotNull List<ActionabilityRange> variantsRanges) {
        this.variants = variants;
        this.variantsRanges = variantsRanges;
    }

    @NotNull
    public Set<String> actionableGenes() {
        Set<String> genes = Sets.newHashSet();
        for (ActionabilityVariant variant : variants) {
            genes.add(variant.gene());
        }
        for (ActionabilityRange range : variantsRanges) {
            genes.add(range.gene());
        }
        return genes;
    }

    public boolean actionableVariants(@NotNull SomaticVariant variant, @NotNull CancerTypeAnalyzer cancerTypeAnalyzer,
            @Nullable String doidsPrimaryTumorLocation) {
        boolean booleanValue = false;
        for (ActionabilityVariant actionabilityVariant : variants) {
            if (variant.gene().equals(actionabilityVariant.gene()) && variant.chromosome().equals(actionabilityVariant.chromosome())
                    && variant.position() == actionabilityVariant.position() && variant.ref().equals(actionabilityVariant.ref())
                    && variant.alt().equals(actionabilityVariant.alt())) {
                LOGGER.info(cancerTypeAnalyzer.foundTumorLocation(actionabilityVariant.cancerType(), doidsPrimaryTumorLocation));

                if (cancerTypeAnalyzer.foundTumorLocation(actionabilityVariant.cancerType(), doidsPrimaryTumorLocation)) {
                    printVariantRow(actionabilityVariant, "yes");
                } else {
                    printVariantRow(actionabilityVariant, "no");
                }
                booleanValue = true;
            } else {
                booleanValue = false;
            }
        }
        return booleanValue;
    }

    private static void printVariantRow(@NotNull ActionabilityVariant variant, @NotNull String isActionable) {
        LOGGER.info(variant.gene() + "\t" + variant.chromosome() + "\t" + variant.position() + "\t" + variant.ref() + "\t" + variant.alt()
                + "\t" + variant.drug() + "\t" + variant.drugsType() + "\t" + variant.cancerType() + "\t" + variant.level() + "\t"
                + variant.response() + "\t" + isActionable);
    }

    public boolean actionableRange(@NotNull SomaticVariant variant, @NotNull CancerTypeAnalyzer cancerTypeAnalyzer,
            @Nullable String doidsPrimaryTumorLocation) {
        boolean booleanValue = false;
        for (ActionabilityRange range : variantsRanges) {
            if (variant.gene().equals(range.gene()) && variant.chromosome()
                    .equals(range.chromosome()) && variant.position() >= range.start() && variant.position() <= range.end()) {

                if (cancerTypeAnalyzer.foundTumorLocation(range.cancerType(), doidsPrimaryTumorLocation)) {
                    printVariantRangeRow(range, "yes");
                } else {
                    printVariantRangeRow(range, "no");
                }
                booleanValue = true;
            } else {
                booleanValue = false;
            }
        }
        return booleanValue;
    }

    private static void printVariantRangeRow(@NotNull ActionabilityRange range, @NotNull String isActionable) {
        LOGGER.info(range.gene() + "\t" + range.chromosome() + "\t" + range.start() + "\t" + range.end() + "\t" + range.drug() + "\t"
                + range.drugsType() + "\t" + range.cancerType() + "\t" + range.level() + "\t" + range.response() + "\t" + isActionable);
    }

    @NotNull
    public static ActionabilityVariantsAnalyzer loadFromFileVariantsAndFileRanges(String fileVariants, String fileRanges)
            throws IOException {
        final List<ActionabilityVariant> variants = new ArrayList<>();
        final List<ActionabilityRange> ranges = new ArrayList<>();
        final List<String> lineVariants = Files.readAllLines(new File(fileVariants).toPath());
        final List<String> lineRanges = Files.readAllLines(new File(fileRanges).toPath());

        for (String line : lineVariants) {
            if (!line.contains("event") || !line.contains("actionability")) {
                variants.add(fromLineVariants(line));
            }
        }

        for (String line : lineRanges) {
            if (!line.contains("event") || !line.contains("actionability")) {
                ranges.add(fromLineRanges(line));
            }
        }
        return new ActionabilityVariantsAnalyzer(variants, ranges);
    }

    @NotNull
    private static ActionabilityVariant fromLineVariants(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityVariant.builder()
                .gene(values[0])
                .chromosome(values[1])
                .position(Long.valueOf(values[2]))
                .ref(values[3])
                .alt(values[4])
                .source(values[5])
                .reference(values[6])
                .drug(values[7])
                .drugsType(values[8])
                .cancerType(values[9])
                .level(values[10])
                .response(values[11])
                .build();
    }

    @NotNull
    private static ActionabilityRange fromLineRanges(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityRange.builder()
                .gene(values[0])
                .chromosome(values[1])
                .start(Long.valueOf(values[2]))
                .end(Long.valueOf(values[3]))
                .source(values[4])
                .reference(values[5])
                .drug(values[6])
                .drugsType(values[7])
                .cancerType(values[8])
                .level(values[9])
                .response(values[10])
                .build();
    }
}
