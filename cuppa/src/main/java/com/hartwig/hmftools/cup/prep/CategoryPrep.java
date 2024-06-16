package com.hartwig.hmftools.cup.prep;

import java.util.List;

import com.hartwig.hmftools.cup.common.CategoryType;

public interface CategoryPrep
{
    CategoryType categoryType();

    List<DataItem> extractSampleData(final String sampleId);
}
