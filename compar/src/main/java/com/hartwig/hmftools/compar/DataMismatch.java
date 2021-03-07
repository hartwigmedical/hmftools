package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.ComparConstants.DATA_DELIM;

import java.util.StringJoiner;

public class DataMismatch
{
    public final String SampleId;
    public final Category Category;
    public final MismatchType MismatchType;
    public final String RefName;
    public final String DiffName;
    public final String RefValue;
    public final String DiffValue;

    public DataMismatch(
            final String sampleId, final Category category, final MismatchType mismatchType, final String refName, final String diffName,
            final String refValue, final String diffValue)
    {
        SampleId = sampleId;
        Category = category;
        MismatchType = mismatchType;
        RefName = refName;
        DiffName = diffName;
        RefValue = refValue;
        DiffValue = diffValue;
    }

    public static String header()
    {
        return "SampleId,Category,MismatchType,RefName,DiffName,RefValue,DiffValue";
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(SampleId);
        sj.add(Category.toString());
        sj.add(MismatchType.toString());
        sj.add(RefName);
        sj.add(DiffName);
        sj.add(RefValue);
        sj.add(DiffValue);
        return sj.toString();
    }
}
