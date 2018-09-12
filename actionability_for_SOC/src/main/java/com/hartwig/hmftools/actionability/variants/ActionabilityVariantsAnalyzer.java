package com.hartwig.hmftools.actionability.variants;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class ActionabilityVariantsAnalyzer {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(ActionabilityVariantsAnalyzer.class);


    @NotNull
    private final List<ActionabilityVariantsSOC> variants;

    public ActionabilityVariantsAnalyzer(@NotNull final List<ActionabilityVariantsSOC> variants) {
        this.variants = variants;
    }

    public boolean actionable(@NotNull SomaticVariant variant) {
//        variants.get(1).gene() == variant.gene();
        return true;
    }



//    @NotNull
//    static ActionabilityRanges fromLineRanges(@NotNull String line) {
//        final String [] values = line.split(DELIMITER);
//        return ImmutableActionabilityRanges.builder()
//                .mutationTranscript(values[0])
//                .chromosome(values[1])
//                .start(values[2])
//                .stop(values[3])
//                .geneTranscript(values[4])
//                .source(values[5])
//                .reference(values[6])
//                .drugsName(values[7])
//                .drugsType(values[8])
//                .cancerType(values[9])
//                .levelSource(values[10])
//                .hmfLevel(values[11])
//                .evidenceType(values[12])
//                .significanceSource(values[13])
//                .hmfResponse(values[14])
//                .build();
//    }
}
