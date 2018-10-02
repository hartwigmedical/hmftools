package com.hartwig.hmftools.actionability.CNVs;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.actionability.cancerTypeMapping.CancerTypeAnalyzer;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ActionabilityCNVsAnalyzer {
    private static final Logger LOGGER = LogManager.getLogger(ActionabilityCNVsAnalyzer.class);
    private static final String DELIMITER = "\t";
    private static final double REL_GAIN = 3.0;
    private static final double ABS_LOSS = 0.5;

    @NotNull
    private final List<ActionabilityCNVs> CNVs;

    public ActionabilityCNVsAnalyzer(@NotNull final List<ActionabilityCNVs> CNVs) {
        this.CNVs = CNVs;
    }

    public boolean actionableCNVs(@NotNull GeneCopyNumber geneCopyNumber, @NotNull CancerTypeAnalyzer cancerTypeAnalyzer,
            @Nullable String doids, @NotNull String primaryTumorLocation) {
        Double minCopyValue = (double) Math.max(0, Math.round(geneCopyNumber.minCopyNumber()));
        boolean booleanValue = false;
        for (int i = 0; i < CNVs.size(); i++) {
            if (CNVs.get(i).cancerType().contains(primaryTumorLocation) && checkCNVType(minCopyValue).equals(CNVs.get(i).cnvType())) {
                if (!CNVs.get(i).cancerType().isEmpty()) {
                    if (cancerTypeAnalyzer.foundTumorLocation(CNVs.get(i).cancerType(), doids)) {
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

    private void printTable(int digit, @NotNull String isActionable) {
        LOGGER.info(CNVs.get(digit).gene() + "\t" + CNVs.get(digit).cnvType() + "\t" + CNVs.get(digit).reference() + "\t" + CNVs.get(digit)
                .drugsName() + "\t" + CNVs.get(digit).drugsType() + "\t" + CNVs.get(digit).cancerType() + "\t" + CNVs.get(digit).hmfLevel()
                + "\t" + CNVs.get(digit).evidenceType() + "\t" + CNVs.get(digit).significanceSource() + "\t" + CNVs.get(digit).hmfResponse()
                + "\t" + isActionable);
    }

    private String checkCNVType(final double copyNumber) {
        Double relativeCopyNumber = copyNumber / 8.0;
        String CNVType = "";
        if (Doubles.lessOrEqual(copyNumber, ABS_LOSS)) {
            CNVType = "Deletion";
        } else if (Doubles.greaterOrEqual(relativeCopyNumber, REL_GAIN)) {
            CNVType = "Amplification";
        }
        return CNVType;
    }

    @NotNull
    public static ActionabilityCNVsAnalyzer loadFromFileCNVs(String fileCNVs) throws IOException {
        final List<ActionabilityCNVs> CNVs = new ArrayList<>();
        final List<String> lineCNVs = Files.readAllLines(new File(fileCNVs).toPath());

        for (int i = 1; i < lineCNVs.size(); i++) {
            fromLineCNVs(lineCNVs.get(i));
            CNVs.add(fromLineCNVs(lineCNVs.get(i)));
        }
        return new ActionabilityCNVsAnalyzer(CNVs);
    }

    @NotNull
    private static ActionabilityCNVs fromLineCNVs(@NotNull String line) {
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
