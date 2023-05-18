package com.hartwig.hmftools.compar;

import java.util.Collections;
import java.util.List;

public class InvalidDataItem implements ComparableItem
{
    private final Category mCategory;

    public InvalidDataItem(final Category category) { mCategory = category; }

    @Override
    public Category category() { return mCategory; }

    @Override
    public String key() { return ""; }

    @Override
    public List<String> displayValues() { return Collections.EMPTY_LIST; }

    @Override
    public boolean reportable()
    {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other) { return false; }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        return null;
    }
}
