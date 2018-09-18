package com.hartwig.hmftools.actionability.CNVs;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class ActionabilityCNVsAnalyzer {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityCNVsAnalyzer.class);
    static final String DELIMITER = "\t";

    @NotNull
    private final List<ActionabilityCNVs> CNVs;

    public ActionabilityCNVsAnalyzer(@NotNull final List<ActionabilityCNVs> CNVs) {
        this.CNVs = CNVs;
    }

    public boolean actionableCNVs(@NotNull GeneCopyNumber geneCopyNumber, @NotNull String primaryTumorLocation){
        Boolean booleanValue = true;
        for (int i=0; i< CNVs.size();i++) {
            if (CNVs.get(i).cancerType().contains(primaryTumorLocation) &&
                    CNVs.get(i).gene().equals(geneCopyNumber.gene()) &&
                    CNVs.get(i).cnvType().equals(geneCopyNumber.minRegionStartSupport())) {
                booleanValue = true;
                LOGGER.info(CNVs.get(i));
            } else {
                booleanValue = false;
            }
        }
        return booleanValue;
    }

    @NotNull
    public static ActionabilityCNVsAnalyzer loadFromFileCNVs(String fileCNVs) throws
            IOException {
        final List<ActionabilityCNVs> CNVs = new ArrayList<>();
        final List<String> lineCNVs = Files.readAllLines(new File(fileCNVs).toPath());

        for (int i = 1; i < lineCNVs.size(); i++) {
            fromLineCNVs(lineCNVs.get(i));
            CNVs.add(fromLineCNVs(lineCNVs.get(i)));
        }
        return new ActionabilityCNVsAnalyzer(CNVs);
    }

    @NotNull
    private static ActionabilityCNVs fromLineCNVs(@NotNull String line){
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityCNVs.builder()
                .gene(values[0])
                .cnvType(values[1])
                .source(values[2])
                .reference(values[3])
                .drugsName(values[4])
                .drugsType(values[5])
                .cancerType(values[6])
                .levelSource(values[7])
                .hmfLevel(values[8])
                .evidenceType(values[9])
                .significanceSource(values[10])
                .hmfResponse(values[11])
                .build();

    }
}
