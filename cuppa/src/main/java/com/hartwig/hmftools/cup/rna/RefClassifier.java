package com.hartwig.hmftools.cup.rna;

import com.hartwig.hmftools.cup.common.CategoryType;

public interface RefClassifier
{
    CategoryType categoryType();

    void buildRefDataSets();
}
