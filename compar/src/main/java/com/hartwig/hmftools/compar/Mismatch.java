package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.CommonUtils.DATA_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.diffsStr;

import java.util.List;
import java.util.StringJoiner;

public class Mismatch
{
    public final ComparableItem RefItem;
    public final ComparableItem OtherItem;
    public final MismatchType MismatchType;
    public final List<String> DiffValues; // list of the form: field(refValue/otherValue)

    public Mismatch(
            final ComparableItem refItem, final ComparableItem otherItem, final MismatchType mismatchType,
            final List<String> diffValues)
    {
        RefItem = refItem;
        OtherItem = otherItem;
        MismatchType = mismatchType;
        DiffValues = diffValues;
    }

    public static String header()
    {
        return "Category,MismatchType,Key,RefValues,OtherValues,Differences";
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);

        if(RefItem != null)
            sj.add(RefItem.category().toString());
        else
            sj.add(OtherItem.category().toString());

        sj.add(MismatchType.toString());

        if(RefItem != null)
            sj.add(RefItem.key());
        else
            sj.add(OtherItem.key());

        if(RefItem != null)
            sj.add(diffsStr(RefItem.displayValues()));
        else
            sj.add("");

        if(OtherItem != null)
            sj.add(diffsStr(OtherItem.displayValues()));
        else
            sj.add("");

        sj.add(diffsStr(DiffValues));

        return sj.toString();
    }
}
