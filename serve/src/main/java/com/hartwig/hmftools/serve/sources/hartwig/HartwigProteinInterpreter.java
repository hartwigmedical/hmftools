package com.hartwig.hmftools.serve.sources.hartwig;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.codon.AminoAcidFunctions;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class HartwigProteinInterpreter {

    @VisibleForTesting
    static final String IGNORE_PROTEIN_ANNOTATION = "-";

    private HartwigProteinInterpreter() {
    }

    @NotNull
    public static String interpretProtein(@NotNull String proteinAnnotation) {
        if (proteinAnnotation.equals(IGNORE_PROTEIN_ANNOTATION)) {
            return Strings.EMPTY;
        }

        String interpretedProteinAnnotation = proteinAnnotation;
        if (interpretedProteinAnnotation.startsWith("p.")) {
            interpretedProteinAnnotation = interpretedProteinAnnotation.substring(2);
        }

        return AminoAcidFunctions.forceSingleLetterProteinAnnotation(interpretedProteinAnnotation);
    }
}
