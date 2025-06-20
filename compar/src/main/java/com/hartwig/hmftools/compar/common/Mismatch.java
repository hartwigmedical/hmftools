package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.compar.common.DiffFunctions.diffsStr;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.compar.ComparableItem;

public record Mismatch(ComparableItem RefItem, ComparableItem NewItem, MismatchType MismatchType, List<String> DiffValues)
{
    // DiffValues is list of the form: field(refValue/otherValue)

    public static String commonHeader(boolean includeSampleId, boolean includeCatagory)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        if(includeSampleId)
            sj.add("SampleId");

        if(includeCatagory)
            sj.add("Category");

        sj.add("MismatchType").add("Key").add("Differences");
        return sj.toString();
    }

    public static String header(boolean includeSampleId)
    {
        return commonHeader(includeSampleId, true) + "\tAllValues";
    }

    public static String commonTsv(boolean writeCategory, final Mismatch mismatch)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        if(writeCategory)
        {
            if(mismatch.RefItem != null)
                sj.add(mismatch.RefItem.category().toString());
            else
                sj.add(mismatch.NewItem.category().toString());
        }

        sj.add(mismatch.MismatchType.toString());

        if(mismatch.RefItem != null)
            sj.add(mismatch.RefItem.key());
        else
            sj.add(mismatch.NewItem.key());

        return sj.toString();
    }

    public String toTsv(boolean writeFieldValues, final List<String> comparedFieldsNames)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(commonTsv(!writeFieldValues, this));

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

    public String toString() { return format("type(%s) item(%s) diffs(%d)",
            MismatchType, RefItem != null ? RefItem.key() : NewItem.key(), DiffValues.size()); }
}
