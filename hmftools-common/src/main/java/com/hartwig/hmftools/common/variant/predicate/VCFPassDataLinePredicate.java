package com.hartwig.hmftools.common.variant.predicate;

import org.jetbrains.annotations.NotNull;

public class VCFPassDataLinePredicate extends VCFDataLinePredicate {

    private static final String SEPARATOR_REGEX = "\t";
    private static final int FILTER_COLUMN = 6;

    private static final String PASS_IDENTIFIER_1 = "PASS";
    private static final String PASS_IDENTIFIER_2 = ".";

    @Override
    public boolean test(@NotNull final String line) {
        boolean isData = false;
        if (super.test(line)) {
            final String[] values = line.split(SEPARATOR_REGEX);
            final String filterValue = values[FILTER_COLUMN];
            if (filterValue.equals(PASS_IDENTIFIER_1) || filterValue.equals(PASS_IDENTIFIER_2)) {
                isData = true;
            }
        }
        return isData;
    }
}
