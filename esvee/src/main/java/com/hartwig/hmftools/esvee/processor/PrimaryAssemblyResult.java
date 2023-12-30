package com.hartwig.hmftools.esvee.processor;

import java.util.List;

import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.assembly.PrimaryAssemblerCounters;
import com.hartwig.hmftools.esvee.sequence.PrimaryAssembly;

@Deprecated
public class PrimaryAssemblyResult
{
    public final Junction OriginalJunction;
    public final PrimaryAssemblerCounters Counters;
    public final List<PrimaryAssembly> Assemblies;

    public PrimaryAssemblyResult(final Junction junction, final PrimaryAssemblerCounters counters, final List<PrimaryAssembly> assemblies)
    {
        OriginalJunction = junction;
        Counters = counters;
        Assemblies = assemblies;
    }
}
