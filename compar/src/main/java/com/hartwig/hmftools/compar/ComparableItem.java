package com.hartwig.hmftools.compar;

import java.util.List;

public interface ComparableItem
{
    Category category();

    boolean matches(final ComparableItem other);

    List<String> findDifferences(final ComparableItem other);

    String description();

}
