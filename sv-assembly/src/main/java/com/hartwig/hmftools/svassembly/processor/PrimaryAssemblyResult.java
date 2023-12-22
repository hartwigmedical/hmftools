package com.hartwig.hmftools.svassembly.processor;

import java.util.List;

import com.hartwig.hmftools.svassembly.Junction;
import com.hartwig.hmftools.svassembly.RegionOfInterest;
import com.hartwig.hmftools.svassembly.assembly.PrimaryAssemblerCounters;
import com.hartwig.hmftools.svassembly.models.PrimaryAssembly;
import com.hartwig.hmftools.svassembly.models.Record;

public class PrimaryAssemblyResult
{
    public final Junction Junction;
    public final PrimaryAssemblerCounters Counters;
    public final List<Record> RecordsOfInterest;
    public final List<PrimaryAssembly> Assemblies;
    public final List<RegionOfInterest> InterestingRegions;

    public PrimaryAssemblyResult(final Junction junction, final PrimaryAssemblerCounters counters, final List<Record> records,
            final List<PrimaryAssembly> assemblies, final List<RegionOfInterest> interestingRegions)
    {
        Junction = junction;
        Counters = counters;
        RecordsOfInterest = records;
        Assemblies = assemblies;
        InterestingRegions = interestingRegions;
    }
}
