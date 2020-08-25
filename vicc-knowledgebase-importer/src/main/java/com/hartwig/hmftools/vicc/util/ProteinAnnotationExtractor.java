package com.hartwig.hmftools.vicc.util;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class ProteinAnnotationExtractor {

    private static final Set<String> CAPITALIZED_STRINGS_TO_UNCAPITALIZE = Sets.newHashSet("DELINS", "DEL", "INS", "DUP", "FS");

    private ProteinAnnotationExtractor() {
    }

    @NotNull
    public static String toProteinAnnotation(@NotNull String featureName) {
        String trimmedName = featureName.trim();
        // Many KBs include the gene in the feature name in some form (eg "EGFR E709K" or "EGFR:E709K").
        // Other KBs put the coding info behind the protein annotation ("V130L (c.388G>C)" rather than the gene in front of it)
        String proteinAnnotation;
        if (trimmedName.contains(" ")) {
            String[] trimmedParts = trimmedName.split(" ");
            proteinAnnotation = trimmedParts[1].contains("(c.") ? trimmedParts[0] : trimmedParts[1];
        } else if (trimmedName.contains(":")) {
            proteinAnnotation = trimmedName.split(":")[1];
        } else {
            proteinAnnotation = trimmedName;
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
        if (proteinAnnotation.endsWith("fs*")) {
            proteinAnnotation = proteinAnnotation.substring(0, proteinAnnotation.length() - 1);
        }

        return proteinAnnotation;
    }
}
