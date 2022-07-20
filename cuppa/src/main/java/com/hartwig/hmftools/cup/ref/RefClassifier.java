package com.hartwig.hmftools.cup.ref;

import com.hartwig.hmftools.common.cuppa.CategoryType;

public interface RefClassifier
{
    CategoryType categoryType();

    void buildRefDataSets();
}
