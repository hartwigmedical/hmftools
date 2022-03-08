package com.hartwig.hmftools.compar;

import java.util.List;

public interface ComparableItem
{
    Category category();

    boolean matches(final ComparableItem other);

    boolean reportable();

    Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel);

    String key();

    List<String> displayValues();
}
