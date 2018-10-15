package com.hartwig.hmftools.common.actionability.somaticvariant;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ActionabilityVariantsAnalyzer {
    private static final Logger LOGGER = LogManager.getLogger(ActionabilityVariantsAnalyzer.class);
    private static final String DELIMITER = "\t";

    @NotNull
    private final List<EvidenceItem> variants;
    private final List<ActionabilityRange> variantsRanges;

    ActionabilityVariantsAnalyzer(@NotNull final List<EvidenceItem> variants, @NotNull List<ActionabilityRange> variantsRanges) {
        this.variants = variants;
        this.variantsRanges = variantsRanges;
    }

    @NotNull
    public Set<String> actionableGenes() {
        Set<String> genes = Sets.newHashSet();
        for (EvidenceItem variant : variants) {
            genes.add(variant.gene());
        }
        for (ActionabilityRange range : variantsRanges) {
            genes.add(range.gene());
        }
        return genes;
    }

    public VariantEvidenceItems actionableVariants(@NotNull SomaticVariant variant, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) {
        List<EvidenceItem> onLabel = Lists.newArrayList();
        List<EvidenceItem> offLabel = Lists.newArrayList();
        for (EvidenceItem evidenceItem : variants) {
            if (variant.gene().equals(evidenceItem.gene()) && variant.chromosome().equals(evidenceItem.chromosome())
                    && variant.position() == evidenceItem.position() && variant.ref().equals(evidenceItem.ref()) && variant.alt()
                    .equals(evidenceItem.alt())) {
                if (actionabilityAnalyzerData.cancerTypeAnalyzer().foundTumorLocation(evidenceItem.cancerType(), doidsPrimaryTumorLocation)) {
                    onLabel.add(evidenceItem);
                    printVariantRow(evidenceItem, "yes");

                } else {
                    offLabel.add(evidenceItem);
                    printVariantRow(evidenceItem, "no");

                }
            }
        }
        return ImmutableVariantEvidenceItems.of(onLabel, offLabel);
    }

    private static void printVariantRow(@NotNull EvidenceItem evidence, @NotNull String isActionable) {
        LOGGER.info(
                evidence.gene() + "\t" + evidence.chromosome() + "\t" + evidence.position() + "\t" + evidence.ref() + "\t" + evidence.alt()
                        + "\t" + evidence.drug() + "\t" + evidence.drugsType() + "\t" + evidence.cancerType() + "\t" + evidence.level()
                        + "\t" + evidence.response() + "\t" + evidence.source() + "\t" + isActionable);
    }

    public ActionabilityRangeEvidenceItem actionableRange(@NotNull SomaticVariant variant, @Nullable String doidsPrimaryTumorLocation,
            @NotNull ActionabilityAnalyzer actionabilityAnalyzerData) {
        List<ActionabilityRange> onLabel = Lists.newArrayList();
        List<ActionabilityRange> offLabel = Lists.newArrayList();
        for (ActionabilityRange range : variantsRanges) {
            if (variant.gene().equals(range.gene()) && variant.chromosome().equals(range.chromosome())
                    && variant.position() >= range.start() && variant.position() <= range.end()) {
                if (actionabilityAnalyzerData.cancerTypeAnalyzer().foundTumorLocation(range.cancerType(), doidsPrimaryTumorLocation)) {
                    printVariantRangeRow(range, "yes");
                    onLabel.add(range);
                } else {
                    offLabel.add(range);
                    printVariantRangeRow(range, "no");
                }
            }
        }
        return ImmutableActionabilityRangeEvidenceItem.of(onLabel, offLabel);
    }

    private static void printVariantRangeRow(@NotNull ActionabilityRange range, @NotNull String isActionable) {
        LOGGER.info(range.gene() + "\t" + range.chromosome() + "\t" + range.start() + "\t" + range.end() + "\t" + range.source() + "\t"
                + range.drug() + "\t" + range.drugsType() + "\t" + range.cancerType() + "\t" + range.level() + "\t" + range.response()
                + "\t" + isActionable);
    }

    @NotNull
    public static ActionabilityVariantsAnalyzer loadFromFileVariantsAndFileRanges(String fileVariants, String fileRanges)
            throws IOException {
        final List<EvidenceItem> variants = new ArrayList<>();
        final List<ActionabilityRange> ranges = new ArrayList<>();
        final List<String> lineVariants = Files.readAllLines(new File(fileVariants).toPath());
        final List<String> lineRanges = Files.readAllLines(new File(fileRanges).toPath());

        for (String lineVariant : lineVariants) {
            if (!lineVariant.contains("event") || !lineVariant.contains("actionability")) {
                variants.add(fromLineVariants(lineVariant));
            }
        }

        for (String lineRange : lineRanges) {
            if (!lineRange.contains("event") || !lineRange.contains("actionability")) {
                ranges.add(fromLineRanges(lineRange));
            }
        }
        return new ActionabilityVariantsAnalyzer(variants, ranges);
    }

    @NotNull
    private static EvidenceItem fromLineVariants(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableEvidenceItem.builder()
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
                .level(values[11])
                .response(values[14])
                .build();
    }

    @NotNull
    private static ActionabilityRange fromLineRanges(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityRange.builder()
                .gene(values[0])
                .chromosome(values[2])
                .start(Long.valueOf(values[3]))
                .end(Long.valueOf(values[4]))
                .source(values[6])
                .reference(values[7])
                .drug(values[8])
                .drugsType(values[9])
                .cancerType(values[10])
                .level(values[12])
                .response(values[15])
                .build();
    }
}