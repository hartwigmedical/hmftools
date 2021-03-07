package com.hartwig.hmftools.compar;

import java.util.List;

public interface Comparator
{
    void processSample(final String sampleId, final List<DataMismatch> mismatches);
}
