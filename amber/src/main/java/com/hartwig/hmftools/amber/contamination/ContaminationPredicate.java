package com.hartwig.hmftools.amber.contamination;

public interface ContaminationPredicate
{
    boolean test(TumorContamination contamination);
}
