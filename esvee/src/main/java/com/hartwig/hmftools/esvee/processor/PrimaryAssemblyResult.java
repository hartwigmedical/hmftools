package com.hartwig.hmftools.esvee.processor;

import java.util.List;

import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.common.RegionOfInterest;
import com.hartwig.hmftools.esvee.assembly.PrimaryAssemblerCounters;
import com.hartwig.hmftools.esvee.sequence.PrimaryAssembly;
import com.hartwig.hmftools.esvee.read.Read;

public class PrimaryAssemblyResult
{
    public final com.hartwig.hmftools.esvee.common.Junction Junction;
    public final PrimaryAssemblerCounters Counters;
    public final List<PrimaryAssembly> Assemblies;

    public PrimaryAssemblyResult(final Junction junction, final PrimaryAssemblerCounters counters, final List<PrimaryAssembly> assemblies)
    {
        Junction = junction;
        Counters = counters;
        Assemblies = assemblies;
    }
}
