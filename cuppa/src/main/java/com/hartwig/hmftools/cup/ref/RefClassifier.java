package com.hartwig.hmftools.cup.ref;

import com.hartwig.hmftools.common.cuppa.CategoryType;

@Deprecated
public interface RefClassifier
{
    CategoryType categoryType();

    boolean buildRefDataSets();
}
