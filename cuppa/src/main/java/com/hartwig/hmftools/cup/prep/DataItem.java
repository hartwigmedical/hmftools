package com.hartwig.hmftools.cup.prep;

public class DataItem
{
    public final Index Index;
    public final String Value;

    public static final String FLD_SOURCE = "Source";
    public static final String FLD_CATEGORY = "Category";
    public static final String FLD_KEY = "Key";
    public static final String FLD_VALUE = "Value";

    public class Index
    {
        public final DataSource Source;
        public final ItemType Type;
        public final String Key;

        public Index(final DataSource source, final ItemType type, final String key)
        {
            Source = source;
            Type = type;
            Key = key;
        }
    }

    public DataItem(final DataSource source, final ItemType type, final String key, final String value)
    {
        Index = new Index(source, type, key);
        Value = value;
    }
}
