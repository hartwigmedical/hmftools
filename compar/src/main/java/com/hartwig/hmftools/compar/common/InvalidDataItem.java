package com.hartwig.hmftools.compar.common;

import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.ItemComparer;

public class InvalidDataItem extends ComparableItem
{
    private final CategoryType mCategory;

    public InvalidDataItem(final CategoryType category) { mCategory = category; }

    @Override
    public CategoryType category() { return mCategory; }

    @Override
    public String key() { return ""; }

    @Override
    public boolean matches(final ComparableItem other) { return false; }

    @Override
    public Mismatch findMismatch(
            final ItemComparer comparer, final ComparableItem other, final MatchLevel matchLevel, final boolean includeMatches)
    {
        return null;
    }

    @Override
    public boolean isValid()
    {
        return false;
    }
}
