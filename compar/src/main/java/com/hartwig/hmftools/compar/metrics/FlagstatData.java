package com.hartwig.hmftools.compar.metrics;

import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public class FlagstatData implements ComparableItem
{
    private final CategoryType mCategory;
    private final BamFlagStats mFlagstat;

    public FlagstatData(final CategoryType category, final BamFlagStats flagstat)
    {
        mCategory = category;
        mFlagstat = flagstat;
    }

    public BamFlagStats flagStats() { return mFlagstat; }

    @Override
    public CategoryType category()
    {
        return mCategory;
    }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }
}
