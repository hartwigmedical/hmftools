package com.hartwig.hmftools.esvee.processor;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.Context;
import com.hartwig.hmftools.esvee.common.RegionOfInterest;
import com.hartwig.hmftools.esvee.models.AlignedAssembly;
import com.hartwig.hmftools.esvee.models.Alignment;
import com.hartwig.hmftools.esvee.models.Record;
import com.hartwig.hmftools.esvee.util.Counter;

import org.jetbrains.annotations.Nullable;

public class SupportScanner
{
    private final Context mContext;
    private final Counter mExtraSupportCounter;

    public SupportScanner(final Context context, final Counter counter)
    {
        mContext = context;
        mExtraSupportCounter = counter;
    }

    @Nullable
    public AlignedAssembly tryRescanSupport(final AlignedAssembly assembly)
    {
        try
        {
            return rescanSupport(assembly);
        }
        catch(final Throwable throwable)
        {
            mContext.Problems.add(new Problem("Failed while re-scanning support", throwable, assembly));
            return null;
        }
    }

    public AlignedAssembly rescanSupport(final AlignedAssembly assembly)
    {
        if(assembly.getAlignmentBlocks().size() < 2)
            return assembly; // No variant, don't care about support

        final List<RegionOfInterest> alignmentRegions = assembly.getAlignmentBlocks().stream()
                .filter(Alignment::isMapped)
                .map(block -> new RegionOfInterest(block.Chromosome, block.ReferenceStartPosition,
                        block.ReferenceStartPosition + block.Length))
                .collect(Collectors.toList());
        final List<RegionOfInterest> deduplicated = RegionOfInterest.tryMerge(alignmentRegions);

        final List<Record> potentialSupport = deduplicated.stream()
                .flatMap(region -> mContext.SAMSource.findReadsContaining(region).stream())
                .distinct()
                .collect(Collectors.toList());

        for(final Record record : potentialSupport)
            tryAddSupport(assembly, record);

        return assembly;
    }

    private void tryAddSupport(final AlignedAssembly assembly, final Record record)
    {
        if(assembly.containsSupport(record))
            return;

        if(assembly.tryAddSupport(mContext.SupportChecker, record))
            mExtraSupportCounter.add(1);
        else
        {
            final Record inverted = record.flipStrand();
            if(assembly.tryAddSupport(mContext.SupportChecker, inverted))
                mExtraSupportCounter.add(1);
        }
    }
}
