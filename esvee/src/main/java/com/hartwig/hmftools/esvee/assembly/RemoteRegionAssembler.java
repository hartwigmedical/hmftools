package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.floor;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_LINK_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyLinker.MATCH_SUBSEQUENCE_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyLinker.formLink;
import static com.hartwig.hmftools.esvee.common.CommonUtils.createByteArray;
import static com.hartwig.hmftools.esvee.common.CommonUtils.initialise;
import static com.hartwig.hmftools.esvee.common.SvConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;
import static com.hartwig.hmftools.esvee.read.ReadUtils.isDiscordantFragment;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.read.BamReader;
import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.types.AssemblyLink;
import com.hartwig.hmftools.esvee.types.AssemblySupport;
import com.hartwig.hmftools.esvee.types.Junction;
import com.hartwig.hmftools.esvee.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.types.JunctionSequence;
import com.hartwig.hmftools.esvee.types.RemoteRegion;
import com.hartwig.hmftools.esvee.types.SupportType;

import org.checkerframework.checker.units.qual.A;

import htsjdk.samtools.SAMRecord;

public class RemoteRegionAssembler
{
    private final AssemblyConfig mConfig;
    private final BamReader mBamReader;

    private RemoteRegion mRemoteRegion;
    private final Map<String,Read> mSourceReads;
    private final List<Read> mMatchedRemoteReads;

    public RemoteRegionAssembler(final AssemblyConfig config, final BamReader bamReader)
    {
        mConfig = config;
        mBamReader = bamReader;

        mRemoteRegion = null;
        mSourceReads = Maps.newHashMap();
        mMatchedRemoteReads = Lists.newArrayList();
    }

    public static boolean isExtensionCandidateAssembly(final JunctionAssembly assembly)
    {
        if(assembly.refBaseTrimLength() < MIN_VARIANT_LENGTH)
            return false;

        int maxExtBaseMatchCount = 0;
        int secondExtBaseMatchCount = 0;
        int remoteJuncMates = 0;

        for(AssemblySupport support : assembly.support())
        {
            if(support.read().isReference()) // for now no reference support
                return false;

            Read read = support.read();

            if(support.type() != SupportType.JUNCTION)
                continue;

            if(support.junctionMatches() > maxExtBaseMatchCount)
            {
                secondExtBaseMatchCount = maxExtBaseMatchCount; // promote the second highest
                maxExtBaseMatchCount = support.junctionMatches();
            }
            else if(support.junctionMatches() > secondExtBaseMatchCount)
            {
                secondExtBaseMatchCount = support.junctionMatches();
            }

            if(read.isSupplementary() || read.isMateUnmapped())
                continue;

            boolean matePastJunction = (read.orientation() == POS_ORIENT) == assembly.isForwardJunction();

            if(isDiscordantFragment(read))
            {
                if(matePastJunction)
                    ++remoteJuncMates;
            }
            else
            {
                if((assembly.isForwardJunction() && read.mateAlignmentStart() > assembly.junction().Position)
                || (!assembly.isForwardJunction() && read.mateAlignmentEnd() < assembly.junction().Position))
                {
                    return false;
                }
            }
        }

        if(secondExtBaseMatchCount < MIN_VARIANT_LENGTH)
            return false;

        if(remoteJuncMates < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return false;

        return true;
    }

    public static boolean assemblyOverlapsRemoteRegion(final JunctionAssembly assembly, final RemoteRegion remoteRegion)
    {
        return positionsOverlap(assembly.minAlignedPosition(), assembly.maxAlignedPosition(), remoteRegion.start(), remoteRegion.end());
    }

    // private static final int REMOTE_REF_BASE_LENGTH_MAX = 300;
    // private static final int REMOTE_REF_BASE_EXTENSION_LENGTH = 100;

    public AssemblyLink tryRemoteAssemblyLink(final JunctionAssembly assembly, final RemoteRegion remoteRegion, final List<Read> sourceReads)
    {
        mRemoteRegion = remoteRegion;

        SV_LOGGER.trace("remote region({}) slice", mRemoteRegion);

        mMatchedRemoteReads.clear();
        sourceReads.forEach(x -> mSourceReads.put(x.id(), x));

        mBamReader.sliceBam(mRemoteRegion.Chromosome, mRemoteRegion.start(), mRemoteRegion.end(), this::processRecord);

        SV_LOGGER.debug("remote region({}) sourcedReads(matched={} unmatched={})",
                mRemoteRegion, mMatchedRemoteReads.size(), mSourceReads.size());

        if(mMatchedRemoteReads.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return null;

        // form a remote ref-based assembly from these reads but without a specific junction
        int remoteRegionStart = mMatchedRemoteReads.stream().mapToInt(x -> x.alignmentStart()).min().orElse(0);
        int remoteRegionEnd = mMatchedRemoteReads.stream().mapToInt(x -> x.alignmentEnd()).max().orElse(0);

        byte[] refGenomeBases = mConfig.RefGenome.getBases(remoteRegion.Chromosome, remoteRegionStart, remoteRegionEnd);

        AssemblyLink assemblyLink = tryAssemblyRemoteRefOverlap(assembly, remoteRegionStart, remoteRegionEnd, refGenomeBases);

        if(assemblyLink == null)
            return null;

        SV_LOGGER.debug("assembly({}) links with remote region({}) matchedReads({})",
                assembly, remoteRegion.toString(), mMatchedRemoteReads.size());

        return assemblyLink;
    }

    private void processRecord(final SAMRecord record)
    {
        Read sourceRead = mSourceReads.remove(record.getReadName());

        if(sourceRead == null)
            return;

        Read remoteRead = new Read(record);

        if(mBamReader.currentIsReferenceSample())
            remoteRead.markReference();

        mMatchedRemoteReads.add(remoteRead);
    }

    public AssemblyLink tryAssemblyRemoteRefOverlap(
            final JunctionAssembly assembly, final int remoteRegionStart, final int remoteRegionEnd, final byte[] refGenomeBases)
    {
        byte[] refBaseQuals = createMinBaseQuals(refGenomeBases.length);

        boolean assemblyReversed = false;
        boolean remoteReversed = false;

        if(assembly.junction().Orientation == mRemoteRegion.orientation())
        {
            if(assembly.junction().isForward())
                remoteReversed = true;
            else
                assemblyReversed = true;
        }

        JunctionSequence assemblySeq = new JunctionSequence(assembly, assemblyReversed, -1);

        JunctionSequence remoteRefSeq = new JunctionSequence(refGenomeBases, refBaseQuals, mRemoteRegion.orientation(), remoteReversed);

        // start with a simple comparison looking for the first sequence around its junction in the second
        String firstJunctionSequence = assemblySeq.junctionSequence();
        int firstJunctionSeqLength = firstJunctionSequence.length();

        // take a smaller sections of the first's junction sequence and try to find their start index in the second sequence
        int juncSeqStartIndex = 0;
        List<int[]> alternativeIndexStarts = Lists.newArrayList();

        int subsequenceLength = MATCH_SUBSEQUENCE_LENGTH;
        int subSeqIterations = (int)floor(firstJunctionSeqLength / subsequenceLength);

        for(int i = 0; i < subSeqIterations; ++i)
        {
            juncSeqStartIndex = i * subsequenceLength;
            int juncSeqEndIndex = juncSeqStartIndex + subsequenceLength;

            if(juncSeqEndIndex >= firstJunctionSeqLength)
                break;

            String firstSubSequence = firstJunctionSequence.substring(juncSeqStartIndex, juncSeqStartIndex + subsequenceLength);

            int secondSubSeqIndex = remoteRefSeq.FullSequence.indexOf(firstSubSequence);

            if(secondSubSeqIndex < 0)
                continue;

            alternativeIndexStarts.add(new int[] {juncSeqStartIndex, secondSubSeqIndex});

            secondSubSeqIndex = remoteRefSeq.FullSequence.indexOf(firstSubSequence, secondSubSeqIndex + subsequenceLength);

            while(secondSubSeqIndex >= 0)
            {
                alternativeIndexStarts.add(new int[] {juncSeqStartIndex, secondSubSeqIndex});
                secondSubSeqIndex = remoteRefSeq.FullSequence.indexOf(firstSubSequence, secondSubSeqIndex + subsequenceLength);
            }
        }

        // now perform a full junction sequence search in the second using the sequence matching logic
        Set<Integer> testedOffsets = Sets.newHashSet();

        int minOverlapLength = min(assembly.extensionLength(), ASSEMBLY_LINK_OVERLAP_BASES);

        for(int[] indexStarts : alternativeIndexStarts)
        {
            // find a comparison range that falls within both sequence's index range around the junction
            int firstJuncSeqMatchIndex = indexStarts[0];
            int secondMatchIndex = indexStarts[1];

            int matchOffset = secondMatchIndex - firstJuncSeqMatchIndex;

            if(testedOffsets.contains(matchOffset))
                continue;

            testedOffsets.add(matchOffset);

            int secondIndexStart = secondMatchIndex - firstJuncSeqMatchIndex;
            int secondIndexEnd = secondIndexStart + firstJunctionSeqLength - 1;

            int firstJuncIndexStart = 0;
            int firstJuncIndexEnd = firstJunctionSeqLength - 1;

            if(secondIndexStart < 0)
            {
                firstJuncIndexStart += -(secondIndexStart);
                secondIndexStart = 0;
            }

            // discount this match if the implied end of the match in the second sequence is past its ref base end
            if(secondIndexEnd >= remoteRefSeq.BaseLength)
                continue;

            int firstIndexStart = firstJuncIndexStart + assemblySeq.junctionSeqStartIndex();
            int firstIndexEnd = min(firstJuncIndexEnd + assemblySeq.junctionSeqStartIndex(), assemblySeq.BaseLength - 1);

            if(secondIndexEnd - secondIndexStart  < minOverlapLength || firstIndexEnd - firstIndexStart  < minOverlapLength)
                continue;

            int mismatchCount = SequenceCompare.compareSequences(
                    assemblySeq.bases(), assemblySeq.baseQuals(), firstIndexStart, firstIndexEnd, assemblySeq.repeatInfo(),
                    remoteRefSeq.bases(), remoteRefSeq.baseQuals(), secondIndexStart, secondIndexEnd, remoteRefSeq.repeatInfo(),
                    PRIMARY_ASSEMBLY_MERGE_MISMATCH);

            if(mismatchCount > PRIMARY_ASSEMBLY_MERGE_MISMATCH)
                continue;

            // now that the index in the remote ref sequence has a match and it is clear where this is in the assembly's extension sequence,
            // the implied junction position in the remote can be determined
            return formLinkWithRemote(
                    assembly, assemblySeq, remoteRefSeq, refGenomeBases, refBaseQuals, remoteRegionStart, remoteRegionEnd,
                    firstIndexStart, secondIndexStart);
        }

        return null;
    }

    private AssemblyLink formLinkWithRemote(
            final JunctionAssembly assembly, final JunctionSequence assemblySeq, final JunctionSequence initialRemoteRefSeq,
            final byte[] refGenomeBases, final byte[] refBaseQuals, final int remoteRegionStart, final  int remoteRegionEnd,
            int firstIndexStart, int secondIndexStart)
    {
        // use the start position of the match to infer where the junction may be in the remote location despite there not being any junction
        // spanning reads at that position
        int assemblyMatchJunctionOffset;
        if(assemblySeq.Reversed)
        {
            int assemblyJuncIndex = assembly.baseLength() - assembly.junctionIndex();
            assemblyMatchJunctionOffset = firstIndexStart - assemblyJuncIndex;
        }
        else
        {
            assemblyMatchJunctionOffset = assembly.junctionIndex() - firstIndexStart;
        }

        int remoteJunctionPosition, remoteJunctionIndex;

        int inferredRemoteJunctionPosition;

        if(mRemoteRegion.isForward())
        {
            remoteJunctionPosition = remoteRegionEnd;
            remoteJunctionIndex = refGenomeBases.length - 1;
            inferredRemoteJunctionPosition = remoteRegionEnd + assemblyMatchJunctionOffset;
        }
        else
        {
            remoteJunctionPosition = remoteRegionStart;
            remoteJunctionIndex = 0;
            inferredRemoteJunctionPosition = remoteRegionStart - assemblyMatchJunctionOffset;
        }

        boolean hasSpanningReads = mMatchedRemoteReads.stream()
                .anyMatch(x -> positionWithin(inferredRemoteJunctionPosition, x.unclippedStart(), x.unclippedEnd()));

        JunctionSequence remoteRefSeq = initialRemoteRefSeq;

        if(assemblyMatchJunctionOffset > 0 && !hasSpanningReads)
        {
            remoteJunctionPosition = inferredRemoteJunctionPosition;

            // suggests that the real junction location is further back into the ref bases
            // test if these bases match the assembly's extension sequence
            int inferredStart, inferredEnd;

            if(mRemoteRegion.isForward())
            {
                inferredStart = remoteRegionEnd;
                inferredEnd = remoteRegionEnd + assemblyMatchJunctionOffset;

            }
            else
            {
                inferredStart = remoteRegionStart - assemblyMatchJunctionOffset;
                inferredEnd = remoteRegionStart;
            }

            String inferredRefBases = mConfig.RefGenome.getBaseString(mRemoteRegion.Chromosome, inferredStart, inferredEnd);
            // byte[] inferredRefBaseQuals = createMinBaseQuals(inferredRefBasees.length);

            // check for a simple match at the assembly's junction
            String assemblyJunctionSequence = assemblySeq.junctionSequence();
            // int assemblyJunctionSeqLength = firstJunctionSequence.length();

            int secondIndexInFirst = assemblyJunctionSequence.indexOf(inferredRefBases);

            if(secondIndexInFirst >= 0)
            {
                firstIndexStart = assemblySeq.junctionSeqStartIndex();

                if(!mRemoteRegion.isForward())
                    secondIndexStart = 0;

                int adjustedRemoteStart = mRemoteRegion.isForward() ? remoteRegionStart : inferredRemoteJunctionPosition;
                int adjustedRemoteEnd = mRemoteRegion.isForward() ? inferredRemoteJunctionPosition : remoteRegionEnd;

                byte[] remoteRefBases = mConfig.RefGenome.getBases(mRemoteRegion.Chromosome, adjustedRemoteStart, adjustedRemoteEnd);
                byte[] remoteRefBaseQuals = createMinBaseQuals(remoteRefBases.length);

                remoteRefSeq = new JunctionSequence(remoteRefBases, remoteRefBaseQuals, mRemoteRegion.orientation(), initialRemoteRefSeq.Reversed);
            }
        }
        else
        {
            // suggests that the break is further into the initially selected ref bases, ie they need to be truncated
            // but those bases are still valid extension bases
            // remoteJunctionIndex needs adjusting?

        }

        List<AssemblySupport> remoteSupport = Lists.newArrayList();

        for(Read read : mMatchedRemoteReads)
        {
            boolean spansJunction = positionWithin(remoteJunctionPosition, read.unclippedStart(), read.unclippedEnd());

            int matchLength = read.basesLength();

            AssemblySupport support = new AssemblySupport(
                    read, spansJunction ? SupportType.JUNCTION : SupportType.DISCORDANT,
                    0, 0, matchLength, 0);

            remoteSupport.add(support);
        }

        Junction remoteJunction = new Junction(mRemoteRegion.Chromosome, remoteJunctionPosition, mRemoteRegion.orientation());

        JunctionAssembly remoteAssembly = new JunctionAssembly(
                remoteJunction, refGenomeBases, refBaseQuals, remoteSupport, Lists.newArrayList());

        remoteAssembly.setJunctionIndex(remoteJunctionIndex);

        remoteAssembly.buildRepeatInfo();

        return formLink(assembly, remoteAssembly, assemblySeq, remoteRefSeq, firstIndexStart, secondIndexStart, false);

    }

    private static byte[] createMinBaseQuals(final int length) { return createByteArray(length, (byte) (LOW_BASE_QUAL_THRESHOLD + 1)); }
}
