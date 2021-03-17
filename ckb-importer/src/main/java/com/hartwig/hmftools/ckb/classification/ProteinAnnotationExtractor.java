package com.hartwig.hmftools.ckb.classification;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventPreprocessor;

import org.jetbrains.annotations.NotNull;

public class ProteinAnnotationExtractor implements EventPreprocessor {

    private static final Set<String> CAPITALIZED_STRINGS_TO_UNCAPITALIZE = Sets.newHashSet("DELINS", "DEL", "INS", "DUP", "FS");

    public ProteinAnnotationExtractor() {
    }

    @NotNull
    @Override
    public String apply(@NotNull String event) {
        String proteinAnnotation = event;
        //TODO how to solve in CKB
        // iClusion tends to use DEL/INS/FS rather than del/ins/fs
        for (String stringToLookFor : CAPITALIZED_STRINGS_TO_UNCAPITALIZE) {
            int position = proteinAnnotation.indexOf(stringToLookFor);
            if (position > 0 && Character.isDigit(proteinAnnotation.charAt(position - 1))) {
                proteinAnnotation = proteinAnnotation.substring(0, position) + stringToLookFor.toLowerCase() + proteinAnnotation.substring(
                        position + stringToLookFor.length());
            }
        }

        return proteinAnnotation;
    }
}
