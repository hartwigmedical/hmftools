package com.hartwig.hmftools.esvee.processor;

import java.util.List;

import com.hartwig.hmftools.esvee.Junction;
import com.hartwig.hmftools.esvee.RegionOfInterest;
import com.hartwig.hmftools.esvee.assembly.PrimaryAssemblerCounters;
import com.hartwig.hmftools.esvee.models.PrimaryAssembly;
import com.hartwig.hmftools.esvee.models.Record;

public class PrimaryAssemblyResult
{
    public final com.hartwig.hmftools.esvee.Junction Junction;
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
