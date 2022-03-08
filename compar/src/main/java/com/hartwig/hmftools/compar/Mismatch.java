package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.CommonUtils.DATA_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.diffsStr;

import java.util.List;
import java.util.StringJoiner;

public class Mismatch
{
    public final ComparableItem RefItem;
    public final ComparableItem NewItem;
    public final MismatchType MismatchType;
    public final List<String> DiffValues; // list of the form: field(refValue/otherValue)

    public Mismatch(
            final ComparableItem refItem, final ComparableItem newItem, final MismatchType mismatchType, final List<String> diffValues)
    {
        RefItem = refItem;
        NewItem = newItem;
        MismatchType = mismatchType;
        DiffValues = diffValues;
    }

    public static String commonHeader()
    {
        return "Category,MismatchType,Key";
    }

    public static String header()
    {
        return commonHeader() + ",RefValues,NewValues,Differences";
    }

    public static String commonCsv(final Mismatch mismatch)
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);

        if(mismatch.RefItem != null)
            sj.add(mismatch.RefItem.category().toString());
        else
            sj.add(mismatch.NewItem.category().toString());

        sj.add(mismatch.MismatchType.toString());

        if(mismatch.RefItem != null)
            sj.add(mismatch.RefItem.key());
        else
            sj.add(mismatch.NewItem.key());

        return sj.toString();
    }

    public String toCsv()
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);

        sj.add(commonCsv(this));

        if(RefItem != null)
            sj.add(diffsStr(RefItem.displayValues()));
        else
            sj.add("");

        if(NewItem != null)
            sj.add(diffsStr(NewItem.displayValues()));
        else
            sj.add("");

        sj.add(diffsStr(DiffValues));

        return sj.toString();
    }
}
