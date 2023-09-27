package com.hartwig.hmftools.common.rna;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public enum RnaQcFilter
{
    PASS,
    FAIL_LOW_COVERAGE,
    WARN_DUPLICATE_RATE,
    WARN_SPLICED_GENE_COVERAGE;

    public static final String QC_FILTER_DELIM = ";";

    public static String qcFiltersToString(final List<RnaQcFilter> statusValues)
    {
        return statusValues.stream().map(x -> x.toString()).collect(Collectors.joining(QC_FILTER_DELIM));
    }

    public static List<RnaQcFilter> qcFiltersFromString(final String qcFiltersStr)
    {
        return Arrays.stream(qcFiltersStr.split(QC_FILTER_DELIM, -1)).map(x -> RnaQcFilter.valueOf(x)).collect(Collectors.toList());
    }
}
