package com.hartwig.healthchecker.common.predicate;

import org.jetbrains.annotations.NotNull;

public class VCFPassDataLinePredicate extends VCFDataLinePredicate {

    private static final String SEPARATOR_REGEX = "\t";
    private static final int FILTER_INDEX = 6;
    private static final String DOT = ".";
    private static final String PASS = "PASS";

    @Override
    public boolean test(@NotNull final String line) {
        boolean isData = false;
        if (super.test(line)) {
            final String[] values = line.split(SEPARATOR_REGEX);
            final String filterValue = values[FILTER_INDEX];
            if (filterValue.equals(PASS) || filterValue.equals(DOT)) {
                isData = true;
            }
        }
        return isData;
    }
}
