package com.hartwig.hmftools.cup.prep;

import java.util.List;

public interface CategoryPrep
{
    CategoryType categoryType();

    List<DataItem> extractSampleData(final String sampleId);
}
