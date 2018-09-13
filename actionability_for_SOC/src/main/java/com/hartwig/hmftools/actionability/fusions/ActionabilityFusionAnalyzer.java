package com.hartwig.hmftools.actionability.fusions;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.actionability.variants.ActionabilityRanges;
import com.hartwig.hmftools.actionability.variants.ActionabilityVariantsAnalyzer;
import com.hartwig.hmftools.actionability.variants.ActionabilityVariantsSOC;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class ActionabilityFusionAnalyzer {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityFusionAnalyzer.class);
    static final String DELIMITER = "\t";

    @NotNull
    private final List<ActionabilityFusionPairs> fusionPairs;

    public ActionabilityFusionAnalyzer(@NotNull final List<ActionabilityFusionPairs> fusionPairs,
            final List<ActionabilityPromiscuosThree> promiscuousThree, final List<ActionabilityPromiscuousFive> promiscuousFive) {
        this.fusionPairs = fusionPairs;
        this.promiscuousThree = promiscuousThree;
        this.promiscuousFive = promiscuousFive;
    }

    private final List<ActionabilityPromiscuosThree> promiscuousThree;
    private final List<ActionabilityPromiscuousFive> promiscuousFive;

    public boolean actionableFusions(@NotNull SomaticVariant variant, @NotNull String primaryTumorLocation, @NotNull int sizeVariants) {
        return true;
    }

    @NotNull
    public static ActionabilityFusionAnalyzer loadFromFileFusions(String fileFusionPairs, String filePromiscuousFive, String filePromiscuousThree) throws
            IOException {
        final List<ActionabilityFusionPairs> fusionPairs = new ArrayList<>();
        final List<ActionabilityPromiscuosThree> promiscuousThree = new ArrayList<>();
        final List<ActionabilityPromiscuousFive> promiscuousFive = new ArrayList<>();

        final List<String> lineFusionPairs =  Files.readAllLines(new File(fileFusionPairs).toPath());
        final List<String> linePromiscuousFive =  Files.readAllLines(new File(filePromiscuousFive).toPath());
        final List<String> linePromiscuousThree =  Files.readAllLines(new File(filePromiscuousThree).toPath());

        for (int i = 1; i< lineFusionPairs.size(); i++) {
            fromLineFusionPairs(lineFusionPairs.get(i));
            fusionPairs.add(fromLineFusionPairs(lineFusionPairs.get(i)));
        }

        for (int i = 1; i< linePromiscuousFive.size(); i++) {
            fromLinePromiscuousFive(linePromiscuousFive.get(i));
            promiscuousFive.add(fromLinePromiscuousFive(linePromiscuousFive.get(i)));
        }

        for (int i = 1; i< linePromiscuousThree.size(); i++) {
            fromLinePromiscuousThree(linePromiscuousThree.get(i));
            promiscuousThree.add(fromLinePromiscuousThree(linePromiscuousThree.get(i)));
        }

        return new ActionabilityFusionAnalyzer(fusionPairs, promiscuousThree, promiscuousFive);
    }

    @NotNull
    static ActionabilityPromiscuosThree fromLinePromiscuousThree(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityPromiscuosThree.builder()
                .gene(values[0])
                .source(values[1])
                .reference(values[2])
                .drugsName(values[3])
                .drugsType(values[4])
                .cancerType(values[5])
                .level(values[6])
                .hmfLevel(values[7])
                .evidenceType(values[8])
                .significanceSource(values[9])
                .hmfResponse(values[10])
                .build();
    }

    @NotNull
    static ActionabilityPromiscuousFive fromLinePromiscuousFive(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityPromiscuousFive.builder()
                .gene(values[0])
                .source(values[1])
                .reference(values[2])
                .drugsName(values[3])
                .drugsType(values[4])
                .cancerType(values[5])
                .levelSource(values[6])
                .hmfLevel(values[7])
                .evidenceType(values[8])
                .significanceSource(values[9])
                .hmfResponse(values[10])
                .build();
    }

    @NotNull
    static ActionabilityFusionPairs fromLineFusionPairs(@NotNull String line){
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityFusionPairs.builder()
                .fiveGene(values[0])
                .threeGene(values[1])
                .source(values[2])
                .reference(values[3])
                .drugsName(values[4])
                .drugsName(values[5])
                .cancerType(values[6])
                .levelSource(values[7])
                .hmfLevel(values[8])
                .evidenceType(values[9])
                .significanceSource(values[10])
                .hmfResponse(values[11])
                .build();
    }

}
