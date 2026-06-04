package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

import java.util.List;

import com.hartwig.hmftools.compar.ComparableItem;

public class Mismatch
{
    public final ComparableItem OldItem;
    public final ComparableItem NewItem;
    public final MismatchType Type;
    public final List<String> DiffValues;

    public Mismatch(final ComparableItem oldItem, final ComparableItem newItem, final MismatchType type, final List<String> diffValues)
    {
        OldItem = oldItem;
        NewItem = newItem;
        Type = type;
        DiffValues = diffValues;
    }

    @Override
    public boolean equals(final Object other)
    {
        if(this == other)
            return true;

        if(!(other instanceof Mismatch))
            return false;

        Mismatch otherMismatch = (Mismatch)other;
        return Type == otherMismatch.Type && OldItem == otherMismatch.OldItem && NewItem == otherMismatch.NewItem;
    }

    public ComparableItem nonNullItem() { return OldItem != null ? OldItem : NewItem; }

    public String toString()
    {
        return format("category(%s) type(%s) item(%s) diffs(%d)",
                nonNullItem().category(), Type, nonNullItem().key(), DiffValues.size());
    }
}
