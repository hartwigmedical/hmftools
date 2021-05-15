package com.hartwig.hmftools.vicc.annotation;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventPreprocessor;

import org.jetbrains.annotations.NotNull;

public class ViccProteinAnnotationExtractor implements EventPreprocessor {

    private static final Set<String> CAPITALIZED_STRINGS_TO_UNCAPITALIZE = Sets.newHashSet("DELINS", "DEL", "INS", "DUP", "FS");

    public ViccProteinAnnotationExtractor() {
    }

    @NotNull
    @Override
    public String apply(@NotNull String event) {
        String trimmedEvent = event.trim();
        // Many KBs include the gene in the feature name in some form (eg "EGFR E709K" or "EGFR:E709K").
        // Other KBs put the coding info behind the protein annotation ("V130L (c.388G>C)" rather than the gene in front of it)
        String proteinAnnotation;
        if (trimmedEvent.contains(" ")) {
            String[] trimmedParts = trimmedEvent.split(" ");
            if (trimmedParts[1].contains("(c.")) {
                proteinAnnotation = trimmedParts[0];
            } else {
                proteinAnnotation = trimmedParts[1];
            }
        } else if (trimmedEvent.contains(":")) {
            proteinAnnotation = trimmedEvent.split(":")[1];
        } else {
            proteinAnnotation = trimmedEvent;
        }

        // Some KBs include "p." in front of the protein annotation
        proteinAnnotation = proteinAnnotation.startsWith("p.") ? proteinAnnotation.substring(2) : proteinAnnotation;

        // Some KBs use DEL/INS/FS rather than del/ins/fs
        for (String stringToLookFor : CAPITALIZED_STRINGS_TO_UNCAPITALIZE) {
            int position = proteinAnnotation.indexOf(stringToLookFor);
            if (position > 0 && Character.isDigit(proteinAnnotation.charAt(position - 1))) {
                proteinAnnotation = proteinAnnotation.substring(0, position) + stringToLookFor.toLowerCase() + proteinAnnotation.substring(
                        position + stringToLookFor.length());
            }
        }

        // Cut out the trailing stop gained in case a stop gained is following on a frameshift
        int trailingStopGained = proteinAnnotation.indexOf("fs*");
        if (trailingStopGained > 0) {
            proteinAnnotation = proteinAnnotation.substring(0, trailingStopGained + 2);
        }

        return proteinAnnotation;
    }
}
