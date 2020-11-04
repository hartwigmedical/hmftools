package com.hartwig.hmftools.cup.common;

import java.util.List;

public interface CuppaClassifier
{
    CategoryType categoryType();
    boolean isValid();

    void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities);

}
