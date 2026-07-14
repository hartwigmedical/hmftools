package com.hartwig.hmftools.compar.common;

import com.hartwig.hmftools.compar.ComparableItem;

public class InvalidDataItem implements ComparableItem
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
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final FieldConfig fieldConfig,
            final boolean includeMatches)
    {
        return null;
    }

    @Override
    public boolean isValid()
    {
        return false;
    }
}
