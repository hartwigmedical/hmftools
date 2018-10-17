package com.hartwig.hmftools.common.actionability.somaticvariant;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SomaticVariantEvidenceAnalyzer {

    private static final String DELIMITER = "\t";

    @NotNull
    private final List<ActionableSomaticVariant> actionableVariants;
    private final List<ActionableRange> actionableRanges;

    SomaticVariantEvidenceAnalyzer(@NotNull final List<ActionableSomaticVariant> actionableVariants,
            @NotNull List<ActionableRange> actionableRanges) {
        this.actionableVariants = actionableVariants;
        this.actionableRanges = actionableRanges;
    }

    @NotNull
    public Set<String> actionableGenes() {
        Set<String> genes = Sets.newHashSet();
        for (ActionableSomaticVariant variant : actionableVariants) {
            genes.add(variant.gene());
        }
        for (ActionableRange range : actionableRanges) {
            genes.add(range.gene());
        }
        return genes;
    }

    @NotNull
    public List<EvidenceItem> evidenceForSomaticVariant(@NotNull SomaticVariant variant, @Nullable String doidsPrimaryTumorLocation,
            @NotNull CancerTypeAnalyzer cancerTypeAnalyzer) {
        List<EvidenceItem> evidenceItems = Lists.newArrayList();
        for (ActionableSomaticVariant actionableVariant : actionableVariants) {
            if (variant.gene().equals(actionableVariant.gene()) && variant.chromosome().equals(actionableVariant.chromosome())
                    && variant.position() == actionableVariant.position() && variant.ref().equals(actionableVariant.ref()) && variant.alt()
                    .equals(actionableVariant.alt())) {
                ImmutableEvidenceItem.Builder evidenceBuilder = fromActionableVariant(actionableVariant);

                evidenceBuilder.event(variant.canonicalHgvsCodingImpact());
                evidenceBuilder.isOnLabel(cancerTypeAnalyzer.foundTumorLocation(actionableVariant.cancerType(), doidsPrimaryTumorLocation));

                evidenceItems.add(evidenceBuilder.build());
            }
        }

        for (ActionableRange actionableRange : actionableRanges) {
            if (variant.gene().equals(actionableRange.gene()) && variant.chromosome().equals(actionableRange.chromosome())
                    && variant.position() >= actionableRange.start() && variant.position() <= actionableRange.end()) {
                ImmutableEvidenceItem.Builder evidenceBuilder = fromActionableRange(actionableRange);

                evidenceBuilder.event(variant.canonicalHgvsCodingImpact());
                evidenceBuilder.isOnLabel(cancerTypeAnalyzer.foundTumorLocation(actionableRange.cancerType(), doidsPrimaryTumorLocation));

                evidenceItems.add(evidenceBuilder.build());
            }
        }
        return evidenceItems;
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder fromActionableVariant(@NotNull ActionableSomaticVariant actionableVariant) {
        return ImmutableEvidenceItem.builder()
                .gene(actionableVariant.gene())
                .reference(actionableVariant.reference())
                .source(actionableVariant.source())
                .drug(actionableVariant.drug())
                .drugsType(actionableVariant.drugsType())
                .level(actionableVariant.level())
                .response(actionableVariant.response());
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder fromActionableRange(@NotNull ActionableRange actionableRange) {
        return ImmutableEvidenceItem.builder()
                .gene(actionableRange.gene())
                .reference(actionableRange.reference())
                .source(actionableRange.source())
                .drug(actionableRange.drug())
                .drugsType(actionableRange.drugsType())
                .level(actionableRange.level())
                .response(actionableRange.response());
    }

    @NotNull
    public static SomaticVariantEvidenceAnalyzer loadFromFileVariantsAndFileRanges(String fileVariants, String fileRanges)
            throws IOException {
        final List<ActionableSomaticVariant> variants = new ArrayList<>();
        final List<ActionableRange> ranges = new ArrayList<>();
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
                .drug(values[7])
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
                .drug(values[8])
                .drugsType(values[9])
                .cancerType(values[10])
                .level(values[12])
                .response(values[15])
                .build();
    }
}