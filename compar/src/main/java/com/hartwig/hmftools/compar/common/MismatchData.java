package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

import java.util.List;

public class MismatchData
{
    public final CategoryType Category;
    public final MismatchType Type;
    public final String Key;
    public final List<String> DiffValues;

    public MismatchData(final CategoryType category, final MismatchType type, final String key, final List<String> diffValues)
    {
        Category = category;
        Type = type;
        Key = key;
        DiffValues = diffValues;
    }

    public boolean matches(final Mismatch mismatch)
    {
        if(mismatch.Type != Type)
            return false;

        CategoryType category = mismatch.nonNullItem().category();

        if(category != Category)
            return false;

        String key = mismatch.nonNullItem().key();

        if(!key.equals(Key))
            return false;

        if(mismatch.DiffValues.size() != DiffValues.size())
            return false;

        for(int i = 0; i < DiffValues.size(); ++i)
        {
            if(!DiffValues.get(i).equals(mismatch.DiffValues.get(i)))
                return false;
        }

        return true;
    }

    public static MismatchData fromMismatch(final Mismatch mismatch)
    {
        CategoryType category = mismatch.nonNullItem().category();
        String key = mismatch.nonNullItem().key();
        return new MismatchData(category, mismatch.Type, key, mismatch.DiffValues);
    }

    public String toString()
    {
        return format("category(%s) type(%s) key(%s) diffs(%d)", Category, Type, Key, DiffValues.size());
    }
}
