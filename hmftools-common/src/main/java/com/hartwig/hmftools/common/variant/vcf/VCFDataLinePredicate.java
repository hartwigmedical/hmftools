package com.hartwig.hmftools.common.variant.vcf;

import java.util.function.Predicate;

import org.jetbrains.annotations.NotNull;

class VCFDataLinePredicate implements Predicate<String> {

    private static final String NON_DATA_LINE_START_IDENTIFIER = "#";

    @Override
    public boolean test(@NotNull final String line) {
        return !line.startsWith(NON_DATA_LINE_START_IDENTIFIER);
    }
}
