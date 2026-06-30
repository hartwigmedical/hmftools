package com.hartwig.hmftools.compar;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public interface ComparableItem
{
    CategoryType category();

    boolean matches(final ComparableItem other);

    Mismatch findMismatch(
            final ComparableItem other, final MatchLevel matchLevel, final FieldConfig fieldConfig, final boolean includeMatches);

    String key();

    List<String> displayValues();

    // values for display only
    default List<String> extraInfoValues() { return Collections.emptyList(); }

    default String geneName() { return ""; }

    default boolean reportable()
    {
        return true;
    }

    default boolean isPass()
    {
        return true;
    }
}
