package com.hartwig.hmftools.esvee.assembly;

import com.hartwig.hmftools.esvee.common.Direction;

@Deprecated
public class JunctionMetrics
{
    public final String JunctionChromosome;
    public final int JunctionPosition;
    public final Direction JunctionDirection;
    public final PrimaryAssemblerCounters Counters;

    public JunctionMetrics(final String junctionChromosome, final int junctionPosition, final Direction junctionDirection,
            final PrimaryAssemblerCounters counters)
    {
        JunctionChromosome = junctionChromosome;
        JunctionPosition = junctionPosition;
        JunctionDirection = junctionDirection;
        Counters = counters;
    }
}
