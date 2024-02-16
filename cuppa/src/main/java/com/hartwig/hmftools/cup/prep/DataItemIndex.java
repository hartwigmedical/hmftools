package com.hartwig.hmftools.cup.prep;

public class DataItemIndex
{
    public final DataSource Source;
    public final ItemType Type;
    public final String Key;

    public DataItemIndex(final DataSource source, final ItemType type, final String key)
    {
        Source = source;
        Type = type;
        Key = key;
    }
}
