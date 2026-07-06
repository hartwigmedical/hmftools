package com.hartwig.hmftools.compar.common.field;

import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

class TestFieldItem<T> implements ComparableItem
{
    public final T Value;

    TestFieldItem(final T value)
    {
        Value = value;
    }

    @Override
    public CategoryType category()
    {
        return null;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        return true;
    }

    @Override
    public String key()
    {
        return "";
    }
}
