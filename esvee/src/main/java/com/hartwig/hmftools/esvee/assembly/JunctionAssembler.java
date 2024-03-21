package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConstants.INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.SvConstants.MIN_VARIANT_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_CONSENSUS_MISMATCH;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READS;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_SUPPORT_MISMATCH;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.processIndelJunction;
import static com.hartwig.hmftools.esvee.common.AssemblyOutcome.DUP_SPLIT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.buildFromJunctionReads;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.expandReferenceBases;
import static com.hartwig.hmftools.esvee.read.ReadFilters.readJunctionExtensionLength;
import static com.hartwig.hmftools.esvee.read.ReadFilters.recordSoftClipsAtJunction;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.common.AssemblySupport;
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

        // find prominent reads to establish the extension sequence, taking any read meeting min soft-clip lengths
        // and repetitive indels

        int consensusMismatch = PRIMARY_ASSEMBLY_CONSENSUS_MISMATCH;

        List<Read> junctionReads = Lists.newArrayList();
        List<Read> extensionReads = Lists.newArrayList();

        Map<Integer,List<Read>> indelLengthReads = Maps.newHashMap();

        for(Read read : rawReads)
        {
            if(!ReadFilters.recordSoftClipsAndCrossesJunction(read, mJunction))
            {
                mNonJunctionReads.add(read);
                continue;
            }

            if(recordSoftClipsAtJunction(read, mJunction) && readJunctionExtensionLength(read, mJunction) >= PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH)
                extensionReads.add(read);

            if((mJunction.isForward() && read.indelImpliedAlignmentEnd() > 0)
            || (mJunction.isReverse() && read.indelImpliedAlignmentStart() > 0))
            {
                buildIndelFrequencies(indelLengthReads, read);
            }

            junctionReads.add(read);
        }

        List<Read> dominantIndelReads = findMaxFrequencyIndelReads(indelLengthReads);

        if(extensionReads.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT && dominantIndelReads.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return Collections.emptyList();

        extensionReads.addAll(dominantIndelReads);

        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads, consensusMismatch);

        if(!extensionSeqBuilder.isValid())
            return Collections.emptyList();

        List<AssemblySupport> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        if(assemblySupport.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return Collections.emptyList();

        JunctionAssembly assembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualitiies(), assemblySupport,
                extensionSeqBuilder.repeatInfo());

        // test other reads against this new assembly
        int supportMismatch = PRIMARY_ASSEMBLY_SUPPORT_MISMATCH;

        int mismatchReads = extensionSeqBuilder.mismatches();

        for(Read read : junctionReads)
        {
            if(assemblySupport.stream().anyMatch(x -> x.read() == read))
                continue;

            if(assembly.checkReadMatches(read, supportMismatch))
                assembly.addJunctionRead(read, false);
            else
                ++mismatchReads;
        }

        assembly.addMismatchReadCount(mismatchReads);

        // deal with mismatch reads by forming a new assembly if they are significant

        expandReferenceBases(assembly);

        return Lists.newArrayList(assembly);
    }

    private void buildIndelFrequencies(final Map<Integer,List<Read>> indelLengthReads, final Read read)
    {
        int maxIndelLength = read.cigarElements().stream()
                .filter(x -> x.getOperator().isIndel())
                .mapToInt(x -> x.getLength()).max().orElse(0);

        if(maxIndelLength >= INDEL_TO_SC_MIN_SIZE_SOFTCLIP)
        {
            List<Read> lengthReads = indelLengthReads.get(maxIndelLength);
            if(lengthReads == null)
            {
                lengthReads = Lists.newArrayList();
                indelLengthReads.put(maxIndelLength, lengthReads);
            }

            lengthReads.add(read);
        }
    }

    private static List<Read> findMaxFrequencyIndelReads(final Map<Integer,List<Read>> indelLengthReads)
    {
        if(indelLengthReads.isEmpty())
            return Collections.emptyList();

        int maxFrequency = 0;
        int indelLength = 0;

        for(Map.Entry<Integer,List<Read>> entry : indelLengthReads.entrySet())
        {
            if(entry.getValue().size() > maxFrequency)
            {
                indelLength = entry.getKey();
                maxFrequency = entry.getValue().size();
            }
        }

        return indelLengthReads.get(indelLength);
    }

    private List<JunctionAssembly> processJunctionOld(final List<Read> rawReads)
    {
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

        if(assembly == null || assembly.baseLength() < MIN_VARIANT_LENGTH)
            return Collections.emptyList();

        int mismatchReads = (int)assembly.support().stream().filter(x -> x.mismatchCount() > PRIMARY_ASSEMBLY_CONSENSUS_MISMATCH).count();

        if(mismatchReads == 0)
        {
            expandReferenceBases(assembly);
            return List.of(assembly);
        }

        // no filtering of the initial sequence and instead rely on the sequence splitting to do this with all initial mismatches preserved
        AssemblyMismatchSplitter splitter = new AssemblyMismatchSplitter(assembly);
        List<JunctionAssembly> assemblies = splitter.splitOnMismatches(MIN_VARIANT_LENGTH, PRIMARY_ASSEMBLY_SPLIT_MIN_READS);

        // extend these sequences in the direction away from the junction
        assemblies.forEach(x -> expandReferenceBases(x));

        return assemblies;
    }
}
