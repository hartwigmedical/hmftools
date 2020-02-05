package com.hartwig.hmftools.protect.actionability_v2;

import java.io.File;

import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.protect.actionability_v2.fusion.ActionableFusion;
import com.hartwig.hmftools.protect.actionability_v2.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.protect.actionability_v2.gene.ActionableGene;
import com.hartwig.hmftools.protect.actionability_v2.gene.ImmutableActionableGene;

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
        for (String lineVariant : lineFusion.subList(1, lineFusion.size())) {
            fusions.add(fromLineFusion(lineVariant));
        }
        return fusions;
    }

    @NotNull
    public static List<ActionableGene> loadFromFileGene(@NotNull String actionableGeneTsv) throws IOException {
        final List<ActionableGene> genes = Lists.newArrayList();
        final List<String> lineGene = Files.readAllLines(new File(actionableGeneTsv).toPath());

        // Skip header line for gene
        for (String lineVariant : lineGene.subList(1, lineGene.size())) {
            genes.add(fromLineGene(lineVariant));
        }
        return genes;
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

}


