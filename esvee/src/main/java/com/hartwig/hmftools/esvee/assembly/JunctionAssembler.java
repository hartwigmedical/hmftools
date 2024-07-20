package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_READ_MAX_MISMATCH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.basesMatch;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.mismatchesPerComparisonLength;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.buildIndelFrequencies;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findIndelExtensionReads;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findMaxFrequencyIndelReads;
import static com.hartwig.hmftools.esvee.assembly.SequenceCompare.calcReadSequenceMismatches;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.readJunctionExtensionLength;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.recordSoftClipsAtJunction;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.INVALID_INDEX;
import static com.hartwig.hmftools.esvee.assembly.types.ReadAssemblyIndices.getJunctionReadExtensionIndices;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.assembly.types.ReadAssemblyIndices;
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
            int mismatchReadCount = 0;

            // test other reads against this new assembly
            for(Read read : junctionReads)
            {
                if(assembly.support().stream().anyMatch(x -> x.cachedRead() == read)) // skip those already added
                    continue;

                if(!canAddJunctionRead(assembly, read))
                    ++mismatchReadCount;
            }

            assembly.addMismatchReadCount(mismatchReadCount);

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

        int secondSupport = assemblySupport.size();
        double secondSupportPerc = secondSupport / (double)firstAssembly.supportCount();

        if(secondSupport < PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT || secondSupportPerc < PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC)
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

    private boolean canAddJunctionRead(final JunctionAssembly assembly, final Read read)
    {
        int readJunctionIndex = read.getReadIndexAtReferencePosition(mJunction.Position, true);

        if(readJunctionIndex == INVALID_INDEX)
            return false;

        ReadAssemblyIndices readAssemblyIndices = getJunctionReadExtensionIndices(
                assembly.junction(), assembly.junctionIndex(), read, readJunctionIndex);

        int assemblyIndex = readAssemblyIndices.AssemblyIndexStart;

        if(assemblyIndex < 0)
            return false;

        int mismatchCount = 0;
        int highQualMatchCount = 0;
        int checkedBaseCount = 0;
        int assemblyBaseLength = assembly.baseLength();

        final byte[] assemblyBases = assembly.bases();
        final byte[] assemblyBaseQuals = assembly.baseQuals();

        for(int i = readAssemblyIndices.ReadIndexStart; i <= readAssemblyIndices.ReadIndexEnd; ++i, ++assemblyIndex)
        {
            if(assemblyIndex < 0)
                continue;

            if(assemblyIndex >= assemblyBaseLength)
                break;

            byte qual = read.getBaseQuality()[i];
            ++checkedBaseCount;

            if(basesMatch(read.getBases()[i], assemblyBases[assemblyIndex], qual, assemblyBaseQuals[assemblyIndex], LOW_BASE_QUAL_THRESHOLD))
            {
                if(qual >= LOW_BASE_QUAL_THRESHOLD)
                    ++highQualMatchCount;
            }
            else
            {
                ++mismatchCount;

                if(mismatchCount > PRIMARY_ASSEMBLY_READ_MAX_MISMATCH)
                    break;
            }
        }

        int permittedMismatches = mismatchesPerComparisonLength(checkedBaseCount);

        if(mismatchCount > permittedMismatches)
        {
            // test again taking repeats into consideration
            mismatchCount = calcReadSequenceMismatches(
                    mJunction.isForward(), assemblyBases, assemblyBaseQuals, assembly.repeatInfo(), read, readJunctionIndex, permittedMismatches);
        }

        if(mismatchCount > permittedMismatches)
            return false;

        assembly.addSupport(read, JUNCTION, readJunctionIndex, highQualMatchCount, mismatchCount, 0);
        return true;
    }

    private void expandReferenceBases(final JunctionAssembly assembly)
    {
        // find the longest length of aligned reference bases extending back from the junction
        int newRefBasePosition = assembly.junction().Position;

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
                newRefBasePosition = min(newRefBasePosition, read.alignmentStart());
                readExtensionDistance = max(readJunctionIndex - read.leftClipLength(), 0);
            }
            else
            {
                newRefBasePosition = max(newRefBasePosition, read.alignmentEnd());
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

        assembly.extendRefBases(newRefBasePosition, Collections.emptyList(), null);

        if(topSupportRead != null)
        {
            extendRefBasesWithJunctionRead(assembly, topSupportRead.cachedRead(), topSupportRead);
        }

        for(SupportRead support : assembly.support())
        {
            if(support == topSupportRead)
                continue;

            extendRefBasesWithJunctionRead(assembly, support.cachedRead(), support);
        }
    }

    private void extendRefBasesWithJunctionRead(final JunctionAssembly assembly, final Read read, final SupportRead existingSupport)
    {
        ReadAssemblyIndices readAssemblyIndices = ReadAssemblyIndices.getJunctionReadRefBaseIndices(
                mJunction, assembly.junctionIndex(), read, existingSupport.junctionReadIndex());

        assembly.addRead(read, readAssemblyIndices, existingSupport.type(), existingSupport);
    }
}
