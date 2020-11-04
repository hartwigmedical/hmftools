package com.hartwig.hmftools.cup.rna;

import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.ref.RefDataConfig;

public interface RefClassifier
{
    CategoryType categoryType();

    void buildRefDataSets();
}
