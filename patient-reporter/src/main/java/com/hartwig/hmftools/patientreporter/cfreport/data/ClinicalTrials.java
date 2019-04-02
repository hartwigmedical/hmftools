package com.hartwig.hmftools.patientreporter.cfreport.data;

import com.hartwig.hmftools.common.actionability.ClinicalTrial;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.stream.Collectors;

public class ClinicalTrials {

    /**
     * Get all trials from the list that are a trial source
     */
    @NotNull
    public static List<ClinicalTrial> filter(@NotNull List<ClinicalTrial> trials) {
        return trials.stream()
                .filter(trial -> trial.source().isTrialSource())
                .collect(Collectors.toList());
    }

    @NotNull
    public static List<ClinicalTrial> sort(@NotNull List<ClinicalTrial> trials) {
        return trials.stream().sorted((item1, item2) -> {
            if (item1.event().equals(item2.event())) {
                return item1.acronym().compareTo(item2.acronym());
            } else {
                return item1.event().compareTo(item2.event());
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static String CCMOId(@NotNull String reference) {
        // Expected format "EXT1 (CCMO)"
        String referenceWithoutParenthesis = reference.replace(")", "");
        String[] splitExtAndCCMO = referenceWithoutParenthesis.split("\\(");
        return splitExtAndCCMO[1];
    }

    @NotNull
    public static String sourceUrl(@NotNull final ClinicalTrial trial) {

        String source = trial.source().sourceName();
        String reference = trial.reference();
        String ext = EXTId(reference);

        switch (source.toLowerCase()) {
            case "iclusion":
                return "https://iclusion.org/hmf/" + ext;
            default:
                return Strings.EMPTY;
        }

    }

    @NotNull
    private static String EXTId(@NotNull String reference) {
        // Expected format "EXT1 (CCMO)"
        String[] splitExtAndCCMO = reference.split("\\(");
        String ext = splitExtAndCCMO[0];
        return ext.substring(3).trim();
    }

}
