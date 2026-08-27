package com.hartwig.hmftools.compar.metrics;

import java.util.List;

import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class FlagstatData extends ComparableItem
{
    private final CategoryType mCategory;
    private final BamFlagStats mFlagstat;

    public FlagstatData(final CategoryType category, final BamFlagStats flagstat, final List<FieldInfo> fields)
    {
        mCategory = category;
        mFlagstat = flagstat;

        addDoubleValue(FlagstatComparer.Fields.MappedProportion.toString(), flagstat.mappedProportion(), fields);
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
