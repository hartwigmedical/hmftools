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


    public boolean actionableVariants(@NotNull SomaticVariant variant, @NotNull String primaryTumorLocation, @NotNull int number){
        Boolean booleanValue;
        LOGGER.info(number);

        LOGGER.info("values");
        LOGGER.info(primaryTumorLocation);
        LOGGER.info(variants.get(number).cancerType());
        LOGGER.info(variant.gene());
        LOGGER.info(variants.get(number).gene());
        LOGGER.info(variant.chromosome());
        LOGGER.info(variants.get(number).chromosome());
        LOGGER.info(Long.toString(variant.position()));
        LOGGER.info(variants.get(number).position());
        LOGGER.info(variant.ref());
        LOGGER.info(variants.get(number).ref());
        LOGGER.info(variant.alt());
        LOGGER.info(variants.get(number).alt());
        LOGGER.info("Check if");
        LOGGER.info(variants.get(number).cancerType().contains(primaryTumorLocation));
        LOGGER.info(variant.gene().equals(variants.get(number).gene()));
        LOGGER.info(variant.chromosome().equals(variants.get(number).chromosome()));
        LOGGER.info(Long.toString(variant.position()).equals(variants.get(number).position()));
        LOGGER.info(variant.ref().equals(variants.get(number).ref()));
        LOGGER.info(variant.alt().equals(variants.get(number).alt()));
        if (variants.get(number).cancerType().contains(primaryTumorLocation) &&
                variant.gene().equals(variants.get(number).gene()) &&
                variant.chromosome().equals(variants.get(number).chromosome()) &&
                Long.toString(variant.position()).equals(variants.get(number).position()) &&
                variant.ref().equals(variants.get(number).ref()) &&
                variant.alt().equals(variants.get(number).alt())) {
            booleanValue =  true;
            LOGGER.info(variants);
        } else {
            booleanValue =  false;
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
        }

        for (int i = 1; i< lineRanges.size(); i++) {
            fromLineVariants(lineRanges.get(i));
            ranges.add(fromLineRanges(lineRanges.get(i)));
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
    private static ActionabilityRanges fromLineRanges(@NotNull String line) {
        final String [] values = line.split(DELIMITER);
        return ImmutableActionabilityRanges.builder()
                .gene(values[0])
                .mutationTranscript(values[1])
                .chromosome(values[2])
                .start(values[3])
                .stop(values[4])
                .geneTranscript(values[5])
                .source(values[6])
                .reference(values[7])
                .drugsName(values[8])
                .drugsType(values[9])
                .cancerType(values[10])
                .levelSource(values[11])
                .hmfLevel(values[12])
                .evidenceType(values[13])
                .significanceSource(values[14])
                .hmfResponse(values[15])
                .build();
    }
}
