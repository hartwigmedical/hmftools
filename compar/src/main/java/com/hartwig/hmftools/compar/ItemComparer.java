package com.hartwig.hmftools.compar;

import java.util.List;

public interface ItemComparer
{
    void processSample(final String sampleId, final List<Mismatch> mismatches);
}
