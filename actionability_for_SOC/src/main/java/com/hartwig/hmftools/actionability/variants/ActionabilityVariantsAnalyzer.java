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
    private static final String DELIMITER = "\t";


    @NotNull
    private final List<ActionabilityVariantsSOC> variants;
    private final List<ActionabilityRanges> variantsRanges;

    public ActionabilityVariantsAnalyzer(@NotNull final List<ActionabilityVariantsSOC> variants, @NotNull List<ActionabilityRanges> variantsRanges) {
        this.variants = variants;
        this.variantsRanges = variantsRanges;
    }


    public boolean actionableVariants(@NotNull SomaticVariant variant, @NotNull String primaryTumorLocation, @NotNull int sizeVariants){
        Boolean booleanValue = true;
        LOGGER.info(sizeVariants);
        for (int i = 0; i < sizeVariants; i ++){
            LOGGER.info(primaryTumorLocation.contains(variants.get(i).cancerType()));
            LOGGER.info(variant.gene().equals(variants.get(i).gene()));
            LOGGER.info(variant.chromosome().equals(variants.get(i).chromosome()));
            LOGGER.info(Long.toString(variant.position()).equals(variants.get(i).position()));
            LOGGER.info(variant.ref().equals(variants.get(i).ref()));
            LOGGER.info(variant.alt().equals(variants.get(i).alt()));
            if (primaryTumorLocation.contains(variants.get(i).cancerType()) &&
                    variant.gene().equals(variants.get(i).gene()) &&
                    variant.chromosome().equals(variants.get(i).chromosome()) &&
                    Long.toString(variant.position()).equals(variants.get(i).position()) &&
                    variant.ref().equals(variants.get(i).ref()) &&
                    variant.alt().equals(variants.get(i).alt())) {
                booleanValue =  true;
                LOGGER.info(variants);
            } else {
                booleanValue =  false;
            }
        }
        return booleanValue;
    }

    public boolean actionableRange(@NotNull SomaticVariant variant, @NotNull String primaryTumorLocation, @NotNull int sizeVariants) {
        return true;
    }

    @NotNull
    public static ActionabilityVariantsAnalyzer loadFromFileVariantsAndFileRanges(String fileVariants, String fileRanges) throws IOException {
        final List<ActionabilityVariantsSOC> variants = new ArrayList<>();
        final List<ActionabilityRanges> ranges = new ArrayList<>();
        final List<String> lineVariants =  Files.readAllLines(new File(fileVariants).toPath());
        final List<String> lineRanges =  Files.readAllLines(new File(fileRanges).toPath());


        for (int i = 1; i< lineVariants.size(); i++) {
            fromLineVariants(lineVariants.get(i));
            variants.add(fromLineVariants(lineVariants.get(i)));
            LOGGER.info(fromLineVariants(lineVariants.get(i)));
        }

        for (int i = 1; i< lineRanges.size(); i++) {
            fromLineVariants(lineRanges.get(i));
            ranges.add(fromLineRanges(lineRanges.get(i)));
            LOGGER.info(fromLineVariants(lineRanges.get(i)));
        }
        return new ActionabilityVariantsAnalyzer(variants, ranges);
    }

    @NotNull
    private static ActionabilityVariantsSOC fromLineVariants(@NotNull String line) {
        final String[] values = line.split(DELIMITER);
        return ImmutableActionabilityVariantsSOC.builder()
                .gene(values[0])
                .chromosome(values[1])
                .position(values[2])
                .ref(values[3])
                .alt(values[4])
                .source(values[5])
                .reference(values[6])
                .drug(values[7])
                .drugsType(values[8])
                .cancerType(values[9])
                .levelSource(values[10])
                .levelHmf(values[11])
                .evidenceType(values[12])
                .significanceSource(values[13])
                .hmfResponse(values[14])
                .build();
    }

    @NotNull
    static ActionabilityRanges fromLineRanges(@NotNull String line) {
        final String [] values = line.split(DELIMITER);
        return ImmutableActionabilityRanges.builder()
                .mutationTranscript(values[0])
                .chromosome(values[1])
                .start(values[2])
                .stop(values[3])
                .geneTranscript(values[4])
                .source(values[5])
                .reference(values[6])
                .drugsName(values[7])
                .drugsType(values[8])
                .cancerType(values[9])
                .levelSource(values[10])
                .hmfLevel(values[11])
                .evidenceType(values[12])
                .significanceSource(values[13])
                .hmfResponse(values[14])
                .build();
    }
}
