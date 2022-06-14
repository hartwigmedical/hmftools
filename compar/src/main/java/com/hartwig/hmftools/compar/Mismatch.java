package com.hartwig.hmftools.compar;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.CommonUtils.DATA_DELIM;
import static com.hartwig.hmftools.compar.CommonUtils.ITEM_DELIM;
import static com.hartwig.hmftools.compar.DiffFunctions.diffsStr;

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
        return "Category,MismatchType,Key,Differences";
    }

    public static String header()
    {
        return commonHeader() + ",AllValues";
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

    public String toCsv(boolean writeFieldValues, final List<String> comparedFieldsNames)
    {
        StringJoiner sj = new StringJoiner(DATA_DELIM);

        sj.add(commonCsv(this));

        sj.add(diffsStr(DiffValues));

        if(writeFieldValues)
        {
            final List<String> refFieldValues = RefItem != null ? RefItem.displayValues() : null;
            final List<String> newFieldValues = NewItem != null ? NewItem.displayValues() : null;
            int fieldCount = refFieldValues != null ? refFieldValues.size() : newFieldValues.size();

            for(int i = 0; i < fieldCount; ++i)
            {
                if(refFieldValues != null)
                    sj.add(refFieldValues.get(i));
                else
                    sj.add("");

                if(newFieldValues != null)
                    sj.add(newFieldValues.get(i));
                else
                    sj.add("");
            }
        }
        else
        {
            ComparableItem item = RefItem != null ? RefItem : NewItem;

            StringJoiner displaySj = new StringJoiner(ITEM_DELIM);

            List<String> itemDisplayValues = item.displayValues();

            for(int i = 0; i < itemDisplayValues.size(); ++i)
            {
                displaySj.add(format("%s=%s", comparedFieldsNames.get(i), itemDisplayValues.get(i)));
            }

            sj.add(displaySj.toString());
        }

        return sj.toString();
    }

    public String toString() { return format("type(%s) item(%) diffs(%d)",
            MismatchType, RefItem != null ? RefItem.key() : NewItem.key(), DiffValues.size()); }
}
