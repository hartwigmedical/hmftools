package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.CommonUtils.findDiffs;

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

    default Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final FieldConfig fieldConfig,
            final boolean includeMatches)
    {
        final List<String> diffs = findDiffs(this, other, fieldConfig.getFields(category()));
        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }

    String key();

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
