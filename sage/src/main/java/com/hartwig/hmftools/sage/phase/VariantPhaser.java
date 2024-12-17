package com.hartwig.hmftools.sage.phase;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;

public interface VariantPhaser
{
    void initialise(final ChrBaseRegion region, final String sample);

    void registeredPhasedVariants(final List<ReadContextCounter> posCounters, final List<ReadContextCounter> negCounters);
}
