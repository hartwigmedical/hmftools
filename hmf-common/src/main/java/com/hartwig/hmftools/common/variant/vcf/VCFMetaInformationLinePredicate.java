package com.hartwig.hmftools.common.variant.vcf;

import java.util.function.Predicate;

import org.jetbrains.annotations.NotNull;

class VCFMetaInformationLinePredicate implements Predicate<String> {

    private static final String META_INFORMATION_LINE_START = "##";

    @Override
    public boolean test(@NotNull final String line) {
        return line.startsWith(META_INFORMATION_LINE_START);
    }
}
