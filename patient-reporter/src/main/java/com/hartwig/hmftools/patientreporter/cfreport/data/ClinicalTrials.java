package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ProtectEvidence;

import org.jetbrains.annotations.NotNull;

public final class ClinicalTrials {

    private ClinicalTrials() {
    }

    @NotNull
    public static List<ProtectEvidence> sort(@NotNull List<ProtectEvidence> trials) {
        return trials.stream().sorted((item1, item2) -> {
            if (item1.genomicEvent().equals(item2.genomicEvent())) {
                return item1.treatment().compareTo(item2.treatment());
            } else {
                return item1.genomicEvent().compareTo(item2.genomicEvent());
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

//    @NotNull
//    public static String sourceUrl(@NotNull ClinicalTrial trial) {
//        if (trial.source() == ActionabilitySource.ICLUSION) {
//            String ext = EXTId(trial.reference());
//            return "https://iclusion.org/hmf/" + ext;
//        }
//
//        return Strings.EMPTY;
//    }

    @NotNull
    private static String EXTId(@NotNull String reference) {
        // Expected format "EXT1 (CCMO)"
        String[] splitExtAndCCMO = reference.split("\\(");
        String ext = splitExtAndCCMO[0];
        return ext.substring(3).trim();
    }

    public static int uniqueEventCount(@NotNull List<ProtectEvidence> trials) {
        Set<String> events = Sets.newHashSet();
        for (ProtectEvidence trial : trials) {
            events.add(trial.genomicEvent());
        }
        return events.size();
    }

    public static int uniqueTrialCount(@NotNull List<ProtectEvidence> trials) {
        Set<String> acronyms = Sets.newHashSet();
        for (ProtectEvidence trial : trials) {
            acronyms.add(trial.treatment());
        }
        return acronyms.size();
    }
}
