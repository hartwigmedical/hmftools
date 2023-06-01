package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileDelimiters.ITEM_DELIM;

import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class SampleData
{
    public final String PatientId;
    public final String TumorId;
    public final List<String> CtDnaSamples;
    public final String VcfTag;

    public SampleData(final String patientId, final String tumorId, final List<String> ctDnaSamples, final String vcfTag)
    {
        PatientId = patientId;
        TumorId = tumorId;
        CtDnaSamples = ctDnaSamples;
        VcfTag = vcfTag;
    }

    public String toString()
    {
        StringJoiner sj = new StringJoiner(", ");
        CtDnaSamples.forEach(x -> sj.add(x));
        return format("patient(%s) tumor(%s) ctDnaSamples(%s)", PatientId, TumorId, sj);
    }

    public static List<String> ctDnaSamplesFromStr(final String ctDnaSamples)
    {
        return Arrays.stream(ctDnaSamples.split(ITEM_DELIM, -1)).collect(Collectors.toList());
    }
}
