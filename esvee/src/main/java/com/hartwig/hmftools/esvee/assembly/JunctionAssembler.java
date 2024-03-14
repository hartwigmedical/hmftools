package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.MIN_VARIANT_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READS;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.processIndelJunction;
import static com.hartwig.hmftools.esvee.common.AssemblyOutcome.DUP_SPLIT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.buildFromJunctionReads;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.expandReferenceBases;
import static com.hartwig.hmftools.esvee.read.ReadFilters.recordSoftClipsAtJunction;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadFilters;

public class JunctionAssembler
{
    private final Junction mJunction;
    private final List<Read> mNonJunctionReads;

    public JunctionAssembler(final Junction junction)
    {
        mJunction = junction;
        mNonJunctionReads = Lists.newArrayList();
    }

    public List<Read> nonJunctionReads() { return mNonJunctionReads; }

    public List<JunctionAssembly> processJunction(final List<Read> rawReads)
    {
        if(mJunction.IndelBased)
            return processIndelJunction(mJunction, mNonJunctionReads, rawReads);

        List<Read> junctionReads = Lists.newArrayList();
        List<Read> preciseJunctionReads = Lists.newArrayList();

        for(Read read : rawReads)
        {
            if(!ReadFilters.recordSoftClipsAndCrossesJunction(read, mJunction))
            {
                mNonJunctionReads.add(read);
                continue;
            }

            if(recordSoftClipsAtJunction(read, mJunction))
                preciseJunctionReads.add(read);

            junctionReads.add(read);
        }

        // CHECK: is read realignment to the specific junction required?

        if(preciseJunctionReads.isEmpty() || junctionReads.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return List.of();

        List<JunctionAssembly> initialAssemblies = createInitialAssemblies(junctionReads);

        // filters

        // by soft-clip length
        List<JunctionAssembly> filteredAssemblies = initialAssemblies.stream()
                .filter(x -> x.extensionLength() >= PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH)
                .collect(Collectors.toList());

        filteredAssemblies.forEach(x -> x.buildRepeatInfo());

        // dedup assemblies for this junction based on overlapping read support
        AssemblyDeduper.dedupJunctionAssemblies(filteredAssemblies);

        if(filteredAssemblies.size() > 1)
            filteredAssemblies.forEach(x -> x.setOutcome(DUP_SPLIT));

        return filteredAssemblies;
    }

    private List<JunctionAssembly> createInitialAssemblies(final List<Read> junctionReads)
    {
        JunctionAssembly assembly = buildFromJunctionReads(mJunction, junctionReads, true);

        // CHECK: why keep these if less than the 32 soft-clip length, since will be dropped anyway
        if(assembly == null || assembly.baseLength() < MIN_VARIANT_LENGTH)
            return Collections.emptyList();

        // no filtering of the initial sequence and instead rely on the sequence splitting to do this with all initial mismatches preserved
        AssemblyMismatchSplitter splitter = new AssemblyMismatchSplitter(assembly);
        List<JunctionAssembly> assemblies = splitter.splitOnMismatches(MIN_VARIANT_LENGTH, PRIMARY_ASSEMBLY_SPLIT_MIN_READS);

        // extend these sequences in the direction away from the junction
        assemblies.forEach(x -> expandReferenceBases(x));

        return assemblies;
    }
}
