package com.hartwig.hmftools.protect.actionability_v2;

import java.io.File;

import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.protect.actionability_v2.range.ImmutableActionableRange;
import com.hartwig.hmftools.protect.actionability_v2.signature.ImmutableActionableSignature;
import com.hartwig.hmftools.protect.actionability_v2.variant.ActionableVariant;
import com.hartwig.hmftools.protect.actionability_v2.fusion.ActionableFusion;
import com.hartwig.hmftools.protect.actionability_v2.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.protect.actionability_v2.gene.ActionableGene;
import com.hartwig.hmftools.protect.actionability_v2.gene.ImmutableActionableGene;
import com.hartwig.hmftools.protect.actionability_v2.range.ActionableRange;
import com.hartwig.hmftools.protect.actionability_v2.signature.ActionableSignature;
import com.hartwig.hmftools.protect.actionability_v2.variant.ImmutableActionableVariant;

import org.jetbrains.annotations.NotNull;

public final class ReadActionabilityFiles {

    private static final String DELIMITER = "\t";

    private ReadActionabilityFiles() {

    }

    @NotNull
    public static List<ActionableFusion> loadFromFileFusion(@NotNull String actionableFusionTsv) throws IOException {
        final List<ActionableFusion> fusions = Lists.newArrayList();
        final List<String> lineFusion = Files.readAllLines(new File(actionableFusionTsv).toPath());

        // Skip header line for fusions
        for (String lineFusions : lineFusion.subList(1, lineFusion.size())) {
            fusions.add(fromLineFusion(lineFusions));
        }
        return fusions;
    }

    @NotNull
    public static List<ActionableGene> loadFromFileGene(@NotNull String actionableGeneTsv) throws IOException {
        final List<ActionableGene> genes = Lists.newArrayList();
        final List<String> lineGene = Files.readAllLines(new File(actionableGeneTsv).toPath());

        // Skip header line for gene
        for (String lineGenes : lineGene.subList(1, lineGene.size())) {
            genes.add(fromLineGene(lineGenes));
        }
        return genes;
    }

    @NotNull
    public static List<ActionableRange> loadFromFileRange(@NotNull String actionableRangeTsv) throws IOException {
        final List<ActionableRange> ranges = Lists.newArrayList();
        final List<String> lineRange = Files.readAllLines(new File(actionableRangeTsv).toPath());

        // Skip header line for range
        for (String lineRanges : lineRange.subList(1, lineRange.size())) {
            ranges.add(fromLineRange(lineRanges));
        }
        return ranges;
    }

    @NotNull
    public static List<ActionableSignature> loadFromFileSignature(@NotNull String actionableSignatureTsv) throws IOException {
        final List<ActionableSignature> signatures = Lists.newArrayList();
        final List<String> lineSignature = Files.readAllLines(new File(actionableSignatureTsv).toPath());

        // Skip header line for signatures
        for (String lineSignatures : lineSignature.subList(1, lineSignature.size())) {
            signatures.add(fromLineSignatures(lineSignatures));
        }
        return signatures;
    }

    @NotNull
    public static List<ActionableVariant> loadFromFileVariant(@NotNull String actionableVariantTsv) throws IOException {
        final List<ActionableVariant> variants = Lists.newArrayList();
        final List<String> lineVariant = Files.readAllLines(new File(actionableVariantTsv).toPath());

        // Skip header line for variant
        for (String lineVariants : lineVariant.subList(1, lineVariant.size())) {
            variants.add(fromLineVariant(lineVariants));
        }
        return variants;
    }

    @NotNull
    private static ActionableFusion fromLineFusion(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableFusion.builder()
                .fusion(values[0])
                .source(values[1])
                .drug(values[2])
                .drugType(values[3])
                .cancerType(values[4])
                .level(values[5])
                .direction(values[6])
                .link(values[7])
                .build();
    }

    @NotNull
    private static ActionableGene fromLineGene(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableGene.builder()
                .gene(values[0])
                .type(values[1])
                .source(values[2])
                .drug(values[3])
                .drugType(values[4])
                .cancerType(values[5])
                .level(values[6])
                .direction(values[7])
                .link(values[8])
                .build();
    }

    @NotNull
    private static ActionableRange fromLineRange(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableRange.builder()
                .gene(values[0])
                .chromosome(values[1])
                .start(values[2])
                .end(values[3])
                .mutationType(values[4])
                .source(values[5])
                .drug(values[6])
                .drugType(values[7])
                .cancerType(values[8])
                .level(values[9])
                .direction(values[10])
                .link(values[11])
                .build();
    }

    @NotNull
    private static ActionableSignature fromLineSignatures(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableSignature.builder()
                .signature(values[0])
                .source(values[1])
                .drug(values[2])
                .drugType(values[3])
                .cancerType(values[4])
                .level(values[5])
                .direction(values[6])
                .link(values[7])
                .build();
    }

    @NotNull
    private static ActionableVariant fromLineVariant(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionableVariant.builder()
                .gene(values[0])
                .chromosome(values[1])
                .position(values[2])
                .ref(values[3])
                .alt(values[4])
                .source(values[5])
                .drug(values[6])
                .drugType(values[7])
                .cancerType(values[8])
                .level(values[9])
                .direction(values[10])
                .link(values[11])
                .build();
    }
}


