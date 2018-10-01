package com.hartwig.hmftools.actionability.variants;

import java.awt.peer.LightweightPeer;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.actionability.cancerTypeMapping.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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

    public boolean actionableVariants(@NotNull SomaticVariant variant, @NotNull CancerTypeAnalyzer cancerTypeAnalyzer,
            @Nullable String doids, @NotNull String primaryTumorLocation) {
        boolean booleanValue = false;
        for (int i = 0; i < variants.size(); i++) {
            if (variants.get(i).cancerType().contains(primaryTumorLocation) && variant.gene().equals(variants.get(i).gene()) && variant.chromosome().equals(variants.get(i).chromosome()) && Long.toString(variant.position())
                    .equals(variants.get(i).position()) && variant.ref().equals(variants.get(i).ref()) && variant.alt()
                    .equals(variants.get(i).alt())) {
                if(variants.get(i).cancerType() != "") {
                    if (cancerTypeAnalyzer.foundTumorLocation(variants.get(i).cancerType(), doids)){
                        printTable(i, "yes");
                    } else {
                        printTable(i, "no");
                    }
                }
                booleanValue = true;
            } else {
                booleanValue = false;
            }
        }
        return booleanValue;
    }

    private void printTable(@NotNull int digit, @NotNull String isActionable) {
        LOGGER.info(variants.get(digit).gene() + "\t" + variants.get(digit).chromosome() + "\t" + variants.get(digit).position() + "\t"
                + variants.get(digit).ref() + "\t" + variants.get(digit).alt() + "\t" + variants.get(digit).drug()
                + "\t" + variants.get(digit).drugsType() + "\t" + variants.get(digit).cancerType()
                + "\t" + variants.get(digit).levelHmf() + "\t" + variants.get(digit).evidenceType() +
                "\t" + variants.get(digit).significanceSource() + "\t" + variants.get(digit).hmfResponse() + "\t" + isActionable);
    }

    public boolean actionableRange(@NotNull SomaticVariant variant, @NotNull CancerTypeAnalyzer cancerTypeAnalyzer,
            @Nullable String doids, @NotNull String primaryTumorLocation) {
        boolean booleanValue = false;
        for (int i=0; i< variantsRanges.size();i++) {
            if (variantsRanges.get(i).cancerType().contains(primaryTumorLocation) &&
                    variant.gene().equals(variantsRanges.get(i).gene()) &&
                    variant.chromosome().equals(variantsRanges.get(i).chromosome()) &&
                    Integer.parseInt(Long.toString(variant.position())) >= Integer.parseInt(variantsRanges.get(i).start()) &&
                    Integer.parseInt(Long.toString(variant.position())) <= Integer.parseInt(variantsRanges.get(i).stop())) {
                if (variantsRanges.get(i).cancerType() != "") {
                    if (cancerTypeAnalyzer.foundTumorLocation(variantsRanges.get(i).cancerType(), doids)) {
                        printTable(i, "yes");
                    } else {
                        printTable(i, "no");
                    }
                }
                booleanValue = true;
            } else {
                booleanValue = false;
            }
        }
        return booleanValue;
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
