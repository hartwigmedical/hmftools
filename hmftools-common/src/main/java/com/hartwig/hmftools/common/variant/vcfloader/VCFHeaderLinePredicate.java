package com.hartwig.hmftools.common.variant.vcfloader;

import java.util.function.Predicate;

import org.jetbrains.annotations.NotNull;

class VCFHeaderLinePredicate implements Predicate<String> {

    private static final String HEADER_LINE_START_IDENTIFIER = "#CHROM";

    @Override
    public boolean test(@NotNull final String line) {
        return line.startsWith(HEADER_LINE_START_IDENTIFIER);
    }
}
