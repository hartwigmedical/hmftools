package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.bam.CigarUtils.maxIndelLength;
import static com.hartwig.hmftools.esvee.AssemblyConstants.INDEL_TO_SC_MIN_SIZE_SOFTCLIP;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_CONSENSUS_MISMATCH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_SUPPORT_MISMATCH;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findIndelExtensionReads;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.expandReferenceBases;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.readJunctionExtensionLength;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.recordSoftClipsAtJunction;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadFilters;

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
        // find prominent reads to establish the extension sequence, taking any read meeting min soft-clip lengths
        // and repetitive indels

        int consensusMismatch = PRIMARY_ASSEMBLY_CONSENSUS_MISMATCH;

        List<Read> junctionReads = Lists.newArrayList();
        List<Read> extensionReads = Lists.newArrayList();

        if(mJunction.IndelBased)
        {
            findIndelExtensionReads(mJunction, rawReads, extensionReads, junctionReads, mNonJunctionReads);
        }
        else
        {
            Map<Integer, List<Read>> indelLengthReads = Maps.newHashMap();

            // the only difference for indel-based junctions is that only the long indels are used to build the consensus extension

            for(Read read : rawReads)
            {
                if(!ReadFilters.recordSoftClipsAndCrossesJunction(read, mJunction))
                {
                    mNonJunctionReads.add(read);
                    continue;
                }

                if(recordSoftClipsAtJunction(read, mJunction)
                && readJunctionExtensionLength(read, mJunction) >= PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH)
                {
                    extensionReads.add(read);
                }

                if((mJunction.isForward() && read.indelImpliedAlignmentEnd() > 0)
                || (mJunction.isReverse() && read.indelImpliedAlignmentStart() > 0))
                {
                    buildIndelFrequencies(indelLengthReads, read);
                }

                junctionReads.add(read);
            }

            List<Read> dominantIndelReads = findMaxFrequencyIndelReads(indelLengthReads);

            extensionReads.addAll(dominantIndelReads);
        }

        if(extensionReads.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return Collections.emptyList();

        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads, consensusMismatch);

        if(!extensionSeqBuilder.isValid())
            return Collections.emptyList();

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        if(assemblySupport.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return Collections.emptyList();

        JunctionAssembly assembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualitiies(), assemblySupport,
                extensionSeqBuilder.repeatInfo());

        // test other reads against this new assembly
        int supportMismatch = PRIMARY_ASSEMBLY_SUPPORT_MISMATCH;

        int mismatchCount = extensionSeqBuilder.mismatches();

        for(Read read : junctionReads)
        {
            if(assemblySupport.stream().anyMatch(x -> x.cachedRead() == read))
                continue;

            if(!assembly.checkAddJunctionRead(read, supportMismatch))
                ++mismatchCount;
        }

        assembly.addMismatchReadCount(mismatchCount);

        // deal with mismatch reads by forming a new assembly if they are significant

        // dedup assemblies for this junction based on overlapping read support
        // AssemblyDeduper.dedupJunctionAssemblies(filteredAssemblies);

        expandReferenceBases(assembly);

        assembly.buildRepeatInfo();

        return Lists.newArrayList(assembly);
    }

    private void buildIndelFrequencies(final Map<Integer,List<Read>> indelLengthReads, final Read read)
    {
        int maxIndelLength = maxIndelLength(read.cigarElements());

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
}
