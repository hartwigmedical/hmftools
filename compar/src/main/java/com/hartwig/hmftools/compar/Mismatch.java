package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.CommonUtils.DATA_DELIM;

import java.util.StringJoiner;

public class Mismatch
{
    public final String SampleId;
    public final Category Category;
    public final MismatchType MismatchType;
    public final String RefSource;
    public final String OtherSource;
    public final String ItemName;
    public final String DiffValue;

    public Mismatch(
            final String sampleId, final Category category, final MismatchType mismatchType, final String refSource, final String otherSource,
            final String itemName, final String diffValue)
    {
        SampleId = sampleId;
        Category = category;
        MismatchType = mismatchType;
        RefSource = refSource;
        OtherSource = otherSource;
        ItemName = itemName;
        DiffValue = diffValue;
    }

    public static String header()
    {
        return "SampleId,Category,MismatchType,RefSource,OtherSource,Item,Diff";
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);
        sj.add(SampleId);
        sj.add(Category.toString());
        sj.add(MismatchType.toString());
        sj.add(RefSource);
        sj.add(OtherSource);
        sj.add(ItemName);
        sj.add(DiffValue);
        return sj.toString();
    }
}
