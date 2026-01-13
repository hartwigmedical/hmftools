package com.hartwig.hmftools.compar;

import java.util.List;

import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public interface ComparableItem
{
    CategoryType category();

    boolean matches(final ComparableItem other);

    Mismatch findMismatch(
            final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds, final boolean includeMatches);

    String key();

    List<String> displayValues();

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
