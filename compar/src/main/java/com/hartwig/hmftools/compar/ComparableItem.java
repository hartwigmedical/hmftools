package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.CommonUtils.findDiffs;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.FieldCheckCache;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.field.FieldValue;

public interface ComparableItem
{
    CategoryType category();

    boolean matches(final ComparableItem other);

    default Mismatch findMismatch(
            final ComparableItem other, final MatchLevel matchLevel, final FieldCheckCache fieldConfig, final boolean includeMatches)
    {
        List<String> diffs;

        if(supportTruthsetData())
        {
            diffs = findDiffs(this, other);
        }
        else
        {
            diffs = findDiffs(this, other, fieldConfig.getFields(category()));
        }

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }

    String key();

    default String geneName() { return ""; }

    default boolean reportable()
    {
        return true;
    }

    default boolean isPass()
    {
        return true;
    }

    default boolean isValid()
    {
        return true;
    }

    // TODO: likely temporary
    default boolean supportTruthsetData() { return false; }

    default Map<String,FieldValue> fieldValues() { return Collections.emptyMap(); }
    default List<String> fieldNames() { return Collections.emptyList(); }
}
