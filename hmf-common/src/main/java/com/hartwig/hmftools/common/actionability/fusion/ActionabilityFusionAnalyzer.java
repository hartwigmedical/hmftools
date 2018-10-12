package com.hartwig.hmftools.common.actionability.fusion;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ActionabilityFusionAnalyzer {
    private static final Logger LOGGER = LogManager.getLogger(ActionabilityFusionAnalyzer.class);
    private static final String DELIMITER = "\t";


    private ActionabilityFusionAnalyzer(@NotNull final List<ActionabilityFusionPairs> fusionPairs,
            final List<ActionabilityPromiscuosThree> promiscuousThree, final List<ActionabilityPromiscuousFive> promiscuousFive) {
        this.fusionPairs = fusionPairs;
        this.promiscuousThree = promiscuousThree;
        this.promiscuousFive = promiscuousFive;
    }
    private final List<ActionabilityFusionPairs> fusionPairs;
    private final List<ActionabilityPromiscuosThree> promiscuousThree;
    private final List<ActionabilityPromiscuousFive> promiscuousFive;

    @NotNull
    public Set<String> actionableGenes() {
        Set<String> genes = Sets.newHashSet();
        for (ActionabilityFusionPairs fusionPairsSet : fusionPairs) {
            genes.add(fusionPairsSet.fiveGene());
        }
        for (ActionabilityPromiscuosThree promiscuousThreeSet: promiscuousThree) {
            genes.add(promiscuousThreeSet.gene());
        }
        for (ActionabilityPromiscuousFive promiscuousFiveSet : promiscuousFive) {
            genes.add(promiscuousFiveSet.gene());
        }
        return genes;
    }

    public void actionableFusions(@NotNull CancerTypeAnalyzer cancerTypeAnalyzer,
            @Nullable String doidsPrimaryTumorLocation) {
        List<FusionEvidenceItems> onLabel = Lists.newArrayList();
        List<FusionEvidenceItems> offLabel = Lists.newArrayList();

      //  return ImmutableFusionEvidenceItems.of(onLabel, offLabel);
    }

    @NotNull
    public static ActionabilityFusionAnalyzer loadFromFileFusions(String fileFusionPairs, String filePromiscuousFive,
            String filePromiscuousThree) throws IOException {
        final List<ActionabilityFusionPairs> fusionPairs = new ArrayList<>();
        final List<ActionabilityPromiscuosThree> promiscuousThree = new ArrayList<>();
        final List<ActionabilityPromiscuousFive> promiscuousFive = new ArrayList<>();

        final List<String> lineFusionPairs = Files.readAllLines(new File(fileFusionPairs).toPath());
        final List<String> linePromiscuousFives = Files.readAllLines(new File(filePromiscuousFive).toPath());
        final List<String> linePromiscuousThrees = Files.readAllLines(new File(filePromiscuousThree).toPath());


        for (String lineFusionPair : lineFusionPairs) {
            if (!lineFusionPair.contains("event") || !lineFusionPair.contains("actionability")) {
                fusionPairs.add(fromLineFusionPairs(lineFusionPair));
            }
        }

        for (String linePromiscuousFive : linePromiscuousFives) {
            if (!linePromiscuousFive.contains("event") || !linePromiscuousFive.contains("actionability")) {
                promiscuousFive.add(fromLinePromiscuousFive(linePromiscuousFive));
            }
        }

        for (String linePromiscuousThree : linePromiscuousThrees) {
            if (!linePromiscuousThree.contains("event") || !linePromiscuousThree.contains("actionability")) {
                promiscuousThree.add(fromLinePromiscuousThree(linePromiscuousThree));
            }
        }

        return new ActionabilityFusionAnalyzer(fusionPairs, promiscuousThree, promiscuousFive);
    }

    @NotNull
    private static ActionabilityPromiscuosThree fromLinePromiscuousThree(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityPromiscuosThree.builder()
                .gene(values[0])
                .source(values[1])
                .reference(values[2])
                .drugsName(values[3])
                .drugsType(values[4])
                .cancerType(values[5])
                .hmfLevel(values[7])
                .hmfResponse(values[10])
                .build();
    }

    @NotNull
    private static ActionabilityPromiscuousFive fromLinePromiscuousFive(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityPromiscuousFive.builder()
                .gene(values[0])
                .source(values[1])
                .reference(values[2])
                .drugsName(values[3])
                .drugsType(values[4])
                .cancerType(values[5])
                .hmfLevel(values[7])
                .hmfResponse(values[10])
                .build();
    }

    @NotNull
    private static ActionabilityFusionPairs fromLineFusionPairs(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityFusionPairs.builder()
                .fiveGene(values[0])
                .threeGene(values[1])
                .source(values[2])
                .reference(values[3])
                .drugsName(values[4])
                .drugsType(values[5])
                .cancerType(values[6])
                .hmfLevel(values[8])
                .hmfResponse(values[11])
                .build();
    }

}
