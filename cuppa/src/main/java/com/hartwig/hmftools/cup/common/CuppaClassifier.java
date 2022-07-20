package com.hartwig.hmftools.cup.common;

import java.util.List;

import com.hartwig.hmftools.common.cuppa.CategoryType;

public interface CuppaClassifier
{
    CategoryType categoryType();

    boolean loadData();

    boolean processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities);

    void close();

}
