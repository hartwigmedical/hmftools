package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

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
import com.hartwig.hmftools.esvee.SvConstants;
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
        List<Read> junctionReads = Lists.newArrayList();

        for(Read read : rawReads)
        {
            if(!ReadFilters.hasAcceptableMapQ(read, SvConstants.READ_FILTER_MIN_JUNCTION_MAPQ))
            {
                // CHECK: any use for these, eg for extension?
                mFilteredReads.add(read);
                continue;
            }

            if(!ReadFilters.recordSoftClipsNearJunction(read, mJunction))
            {
                mNonJunctionReads.add(read);
                continue;
            }

            if(!ReadFilters.isRecordAverageQualityPastJunctionAbove(read, mJunction, SvConstants.AVG_BASE_QUAL_THRESHOLD))
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

    private String nextAssemblyName()
    {
        // consider naming based on initial length and support? try to make typically deterministic
        return String.format("%s:%s%s:%s", mJunction.Chromosome, mJunction.Position,
                mJunction.direction() == Direction.FORWARDS ? "F" : "R", mNextAssemblyNumber++);
    }

    private List<JunctionAssembly> createInitialAssemblies(final List<Read> junctionReads)
    {
        JunctionAssembly junctionSequence = buildFromJunctionReads(mJunction, junctionReads, true);

        if(junctionSequence.baseLength() < PRIMARY_ASSEMBLY_MIN_LENGTH)
            return Collections.emptyList();

        // no filtering of the initial sequence and instead rely on the sequence splitting to do this with all initial mismatches preserved
        AssemblyMismatchSplitter splitter = new AssemblyMismatchSplitter(junctionSequence);
        List<JunctionAssembly> junctionSequences = splitter.splitOnMismatches(PRIMARY_ASSEMBLY_MIN_LENGTH);

        // extend these sequences in the direction away from the junction
        junctionSequences.forEach(x -> expandReferenceBases(x));

        return junctionSequences;
    }
}
