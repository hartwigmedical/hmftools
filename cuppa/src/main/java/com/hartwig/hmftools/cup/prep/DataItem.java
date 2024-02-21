package com.hartwig.hmftools.cup.prep;

public class DataItem
{
    public final DataItemIndex Index;
    public final String Value;

    public static final String FLD_SOURCE = "Source";
    public static final String FLD_CATEGORY = "Category";
    public static final String FLD_KEY = "Key";
    public static final String FLD_VALUE = "Value";

    public DataItem(final DataSource source, final ItemType type, final String key, final String value)
    {
        Index = new DataItemIndex(source, type, key);
        Value = value;
    }
}
