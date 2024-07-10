package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.buildIndelFrequencies;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findIndelExtensionReads;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findMaxFrequencyIndelReads;
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

        List<Read> junctionReads = Lists.newArrayList();
        List<Read> extensionReads = Lists.newArrayList();

        if(mJunction.IndelBased)
        {
            findIndelExtensionReads(mJunction, rawReads, extensionReads, junctionReads, mNonJunctionReads);
        }
        else
        {
            Map<Integer,List<Read>> indelLengthReads = Maps.newHashMap();

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

        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads);

        if(!extensionSeqBuilder.isValid())
            return Collections.emptyList();

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        if(assemblySupport.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return Collections.emptyList();

        JunctionAssembly firstAssembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualities(), assemblySupport,
                extensionSeqBuilder.repeatInfo());

        List<JunctionAssembly> assemblies = Lists.newArrayList(firstAssembly);

        // test for a second well-supported, alternative assembly at the same junction
        JunctionAssembly secondAssembly = checkSecondAssembly(extensionSeqBuilder.mismatchReads(), firstAssembly);

        if(secondAssembly != null)
            assemblies.add(secondAssembly);

        for(JunctionAssembly assembly : assemblies)
        {
            int mismatchCount = 0;

            // test other reads against this new assembly
            for(Read read : junctionReads)
            {
                if(assembly.support().stream().anyMatch(x -> x.cachedRead() == read)) // skip those already added
                    continue;

                if(!assembly.checkAddJunctionRead(read))
                    ++mismatchCount;
            }

            assembly.addMismatchReadCount(mismatchCount);

            expandReferenceBases(assembly);

            assembly.buildRepeatInfo();
        }

        return assemblies;
    }

    private JunctionAssembly checkSecondAssembly(final List<Read> extensionReads, final JunctionAssembly firstAssembly)
    {
        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads);

        if(!extensionSeqBuilder.isValid())
            return null;

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        if(assemblySupport.size() < PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT)
            return null;

        JunctionAssembly newAssembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualities(), assemblySupport,
                extensionSeqBuilder.repeatInfo());

        if(newAssembly.extensionLength() < PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH)
            return null;

        // perform a final sequence comparison check with more liberal comparison tests
        boolean closeMatch = SequenceCompare.matchedAssemblySequences(firstAssembly, newAssembly);
        return !closeMatch ? newAssembly : null;
    }

    private void expandReferenceBases(final JunctionAssembly assembly)
    {
        // find the longest length of aligned reference bases extending back from the junction
        int minAlignedPosition = assembly.minAlignedPosition();
        int maxAlignedPosition = assembly.maxAlignedPosition();

        int junctionPosition = assembly.junction().Position;
        boolean junctionIsForward = assembly.junction().isForward();

        int maxDistanceFromJunction = 0;

        SupportRead topSupportRead = null;

        for(SupportRead support : assembly.support())
        {
            Read read = support.cachedRead();
            int readJunctionIndex = read.getReadIndexAtReferencePosition(junctionPosition, true);

            // for positive orientations, if read length is 10, and junction index is at 4, then extends with indices 0-3 ie 4
            // for negative orientations, if read length is 10, and junction index is at 6, then extends with indices 7-9 ie 4
            int readExtensionDistance;

            if(junctionIsForward)
            {
                minAlignedPosition = min(minAlignedPosition, read.alignmentStart());
                readExtensionDistance = max(readJunctionIndex - read.leftClipLength(), 0);
            }
            else
            {
                maxAlignedPosition = max(maxAlignedPosition, read.alignmentEnd());
                readExtensionDistance = max(read.basesLength() - readJunctionIndex - 1 - read.rightClipLength(), 0);
            }

            assembly.checkAddRefSideSoftClip(read);

            maxDistanceFromJunction = max(maxDistanceFromJunction, readExtensionDistance);

            if(topSupportRead == null)
            {
                topSupportRead = support; // will be the initial
            }
            else
            {
                // select the read with the fewest SNVs in the aligned bases that also has the equal least number of junction mismatches
                if(support.junctionMismatches() <= topSupportRead.junctionMismatches()
                && support.junctionMatches() >= PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH
                && read.snvCount() < topSupportRead.cachedRead().snvCount())
                {
                    topSupportRead = support;
                }
            }
        }

        assembly.extendRefBases(maxDistanceFromJunction, minAlignedPosition, maxAlignedPosition, null);

        if(topSupportRead != null)
        {
            assembly.extendRefBasesWithJunctionRead(topSupportRead.cachedRead(), topSupportRead);
        }

        for(SupportRead support : assembly.support())
        {
            if(support == topSupportRead)
                continue;

            assembly.extendRefBasesWithJunctionRead(support.cachedRead(), support);
        }
    }
}
