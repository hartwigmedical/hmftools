package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConstants.MIN_INDEL_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.buildFromJunctionReads;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.expandReferenceBases;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.IndelCoords;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.Direction;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.read.ReadFilters;

public class JunctionAssembler
{
    private final SvConfig mConfig;

    private final Junction mJunction;

    private final List<Read> mNonJunctionReads;
    private final List<Read> mFilteredReads;

    private int mNextAssemblyNumber = 1;

    public JunctionAssembler(final SvConfig config, final Junction junction)
    {
        mConfig = config;
        mJunction = junction;
        mFilteredReads = Lists.newArrayList();
        mNonJunctionReads = Lists.newArrayList();
    }

    public List<Read> nonJunctionReads() { return mNonJunctionReads; }

    public List<JunctionAssembly> processJunction(final List<Read> rawReads)
    {
        if(mJunction.IndelBased)
            return processIndelJunction(rawReads);

        List<Read> junctionReads = Lists.newArrayList();

        for(Read read : rawReads)
        {
            if(!ReadFilters.recordSoftClipsNearJunction(read, mJunction))
            {
                mNonJunctionReads.add(read);
                continue;
            }

            if(!ReadFilters.isAboveBaseQualAvgPastJunctionThreshold(read, mJunction))
            {
                mFilteredReads.add(read);
                continue;
            }

            junctionReads.add(read);
        }

        // CHECK: is read realignment to the specific junction required?

        if(junctionReads.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
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

        return filteredAssemblies;
    }

    private List<JunctionAssembly> createInitialAssemblies(final List<Read> junctionReads)
    {
        JunctionAssembly junctionSequence = buildFromJunctionReads(mJunction, junctionReads, true);

        // CHECK: why keep these if less than the 32 soft-clip length, since will be dropped anyway
        if(junctionSequence.baseLength() < PRIMARY_ASSEMBLY_MIN_LENGTH)
            return Collections.emptyList();

        // no filtering of the initial sequence and instead rely on the sequence splitting to do this with all initial mismatches preserved
        AssemblyMismatchSplitter splitter = new AssemblyMismatchSplitter(junctionSequence);
        List<JunctionAssembly> junctionSequences = splitter.splitOnMismatches(PRIMARY_ASSEMBLY_MIN_LENGTH);

        // extend these sequences in the direction away from the junction
        junctionSequences.forEach(x -> expandReferenceBases(x));

        return junctionSequences;
    }

    public List<JunctionAssembly> processIndelJunction(final List<Read> rawReads)
    {
        List<Read> indelReads = Lists.newArrayList();
        List<Read> shortIndelReads = Lists.newArrayList();
        List<Read> softClippedReads = Lists.newArrayList();

        for(Read read : rawReads)
        {
            IndelCoords indelCoords = read.indelCoords();

            if(indelCoords != null)
            {
                // must match junction exactly to be considered for support
                if(!indelCoords.matchesJunction(mJunction.Position, mJunction.Orientation))
                    continue;

                if(indelCoords.Length >= MIN_INDEL_LENGTH)
                    indelReads.add(read);
                else
                    shortIndelReads.add(read);
            }
            else
            {
                if(!ReadFilters.recordSoftClipsNearJunction(read, mJunction))
                {
                    mNonJunctionReads.add(read);
                    continue;
                }

                if(!ReadFilters.isAboveBaseQualAvgPastJunctionThreshold(read, mJunction))
                {
                    mFilteredReads.add(read);
                    continue;
                }

                softClippedReads.add(read);
            }
        }

        // CHECK: is read realignment to the specific junction required?

        if(indelReads.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return List.of();

        JunctionAssembly assembly = IndelBuilder.buildFromIndelReads(mJunction, indelReads, shortIndelReads, softClippedReads);

        if(assembly.supportCount() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return Collections.emptyList();

        expandReferenceBases(assembly);

        assembly.buildRepeatInfo();

        return Lists.newArrayList(assembly);
    }

    private String nextAssemblyName()
    {
        // consider naming based on initial length and support? try to make typically deterministic
        return String.format("%s:%s%s:%s", mJunction.Chromosome, mJunction.Position,
                mJunction.direction() == Direction.FORWARDS ? "F" : "R", mNextAssemblyNumber++);
    }
}
