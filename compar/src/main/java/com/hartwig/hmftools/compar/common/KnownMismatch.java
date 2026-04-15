package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class KnownMismatch
{
    public final CategoryType Category;
    public final MismatchType Type;
    public final String Key;
    public final List<String> DiffValues;
    public final CurationInfo Curation;

    public KnownMismatch(
            final CategoryType category, final MismatchType type, final String key, final List<String> diffValues,
            final CurationInfo curation)
    {
        Category = category;
        Type = type;
        Key = key;
        DiffValues = diffValues;
        Curation = curation;
    }

    public static Map<String,CurationInfo> matchCurations(final Mismatch mismatch, final List<KnownMismatch> curations)
    {
        Map<String,CurationInfo> curationsMap = null;

        for(KnownMismatch knownMismatch : curations)
        {
            if(knownMismatch.Curation.Type == CurationType.NONE || knownMismatch.Curation.Type == CurationType.INVALID)
                continue;

            Map<String,CurationInfo> matchedCurations = knownMismatch.matchCurations(mismatch);

            if(matchedCurations.isEmpty())
                continue;

            if(curationsMap == null)
                curationsMap = Maps.newHashMap();

            curationsMap.putAll(matchedCurations);
        }

        return curationsMap == null ? Collections.emptyMap() : curationsMap;
    }

    public Map<String,CurationInfo> matchCurations(final Mismatch mismatch)
    {
        if(mismatch.Type != Type)
            return Collections.emptyMap();

        CategoryType category = mismatch.nonNullItem().category();

        if(category != Category)
            return Collections.emptyMap();

        String key = mismatch.nonNullItem().key();

        if(!key.equals(Key))
            return Collections.emptyMap();

        Map<String,CurationInfo> curationMap = Maps.newHashMap();

        if(mismatch.DiffValues.isEmpty() && DiffValues.isEmpty())
        {
            curationMap.put(Type.toString(), Curation);
        }
        else
        {
            for(String diff : DiffValues)
            {
                if(mismatch.DiffValues.contains(diff))
                    curationMap.put(diff, Curation);
            }
        }

        return curationMap;
    }

    @Deprecated
    public boolean matchesMismatch(final Mismatch mismatch)
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

    public String toString()
    {
        return format("category(%s) type(%s) key(%s) diffs(%d) known(%s)", Category, Type, Key, DiffValues.size(), Curation.Type);
    }
}
