package com.hartwig.hmftools.cobalt.calculations;

import java.util.Collection;

public interface ResultsNormaliser
{
    void recordValue(BamRatio bamRatio);
    void applyNormalisation(BamRatio bamRatio);
}
