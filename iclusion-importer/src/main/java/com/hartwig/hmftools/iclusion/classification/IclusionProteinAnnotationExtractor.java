package com.hartwig.hmftools.iclusion.classification;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventPreprocessor;

import org.jetbrains.annotations.NotNull;

public class IclusionProteinAnnotationExtractor implements EventPreprocessor {

    private static final Set<String> CAPITALIZED_STRINGS_TO_UNCAPITALIZE = Sets.newHashSet("DELINS", "DEL", "INS", "DUP", "FS");

    public IclusionProteinAnnotationExtractor() {
    }

    @NotNull
    @Override
    public String apply(@NotNull String event) {
        String proteinAnnotation = event;

        if (proteinAnnotation.contains("ins")) {
            proteinAnnotation = proteinAnnotation.replace("-", "_");
        }

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
