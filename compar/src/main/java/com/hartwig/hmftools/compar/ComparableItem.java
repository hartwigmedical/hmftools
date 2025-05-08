package com.hartwig.hmftools.compar;

import java.util.List;

import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public interface ComparableItem
{
    Category category();

    boolean matches(final ComparableItem other);

    boolean reportable();

    boolean isPass();

    Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches);

    String key();

    List<String> displayValues();
}
