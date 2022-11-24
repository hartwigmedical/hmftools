package com.hartwig.hmftools.peach.haplotype;

import com.hartwig.hmftools.peach.event.HaplotypeEvent;

public interface Haplotype
{
    boolean isRelevantFor(HaplotypeEvent event);
}
