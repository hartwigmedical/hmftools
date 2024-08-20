package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.Arrays.reverseArray;
import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_LINK_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.AssemblyConstants.MATCH_SUBSEQUENCE_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.calcTrimmedExtensionBaseLength;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.createLowBaseQuals;
import static com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker.findBestSequenceMatch;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.createMinBaseQuals;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.esvee.assembly.SequenceCompare;
import com.hartwig.hmftools.esvee.assembly.read.BamReader;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.JunctionSequence;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

import htsjdk.samtools.SAMRecord;

public class RemoteRegionAssembler
{
    private final RefGenomeInterface mRefGenome;
    private final BamReader mBamReader;

    private RemoteRegion mRemoteRegion;
    private final Set<String> mSourceReadIds;
    private final List<Read> mMatchedRemoteReads;

    private int mTotalRemoteReadsSearch;
    private int mTotalRemoteReadsMatched;

    public RemoteRegionAssembler(final RefGenomeInterface refGenome, final BamReader bamReader)
    {
        mRefGenome = refGenome;
        mBamReader = bamReader;

        mRemoteRegion = null;
        mSourceReadIds = Sets.newHashSet();
        mMatchedRemoteReads = Lists.newArrayList();

        mTotalRemoteReadsSearch = 0;
        mTotalRemoteReadsMatched = 0;
    }

    public int totalRemoteReadsSearch() { return mTotalRemoteReadsSearch; }
    public int totalRemoteReadsMatched() { return mTotalRemoteReadsMatched; }

    public static boolean isExtensionCandidateAssembly(final JunctionAssembly assembly)
    {
        // apply some filters to limit the number of assemblies which attempt to find a remote discordant match
        if(assembly.refBaseTrimLength() < MIN_VARIANT_LENGTH)
            return false;

        if(assembly.stats().SoftClipSecondMaxLength < MIN_VARIANT_LENGTH)
            return false;

        if(assembly.stats().JuncMateDiscordantRemote < ASSEMBLY_MIN_READ_SUPPORT)
            return false;

        // check for sufficient diversity in the extension bases
        int trimmedExtBaseLength = calcTrimmedExtensionBaseLength(assembly);

        if(trimmedExtBaseLength < MIN_VARIANT_LENGTH)
            return false;

        return true;
    }

    public static boolean assemblyOverlapsRemoteRegion(final JunctionAssembly assembly, final RemoteRegion remoteRegion)
    {
        if(!assembly.junction().Chromosome.equals(remoteRegion.Chromosome))
            return false;

        if(assembly.isForwardJunction())
            return positionsOverlap(assembly.refBasePosition(), assembly.junction().Position, remoteRegion.start(), remoteRegion.end());
        else
            return positionsOverlap(assembly.junction().Position, assembly.refBasePosition(), remoteRegion.start(), remoteRegion.end());
    }

    public List<Read> extractRemoteReads(final RemoteRegion remoteRegion)
    {
        mRemoteRegion = remoteRegion;

        mSourceReadIds.clear();
        mSourceReadIds.addAll(remoteRegion.readIds());

        mTotalRemoteReadsSearch += remoteRegion.readIds().size();

        if(mBamReader != null)
        {
            mMatchedRemoteReads.clear();

            SV_LOGGER.trace("remote region({}) slice", mRemoteRegion);

            mBamReader.sliceBam(mRemoteRegion.Chromosome, mRemoteRegion.start(), mRemoteRegion.end(), this::processRecord);

            SV_LOGGER.trace("remote region({}) sourcedReads(matched={} unmatched={})",
                    mRemoteRegion, mMatchedRemoteReads.size(), mSourceReadIds.size());

            mTotalRemoteReadsMatched += mMatchedRemoteReads.size();
        }

        // ignore supplementaries since their bases provide no new assembly sequence information
        return mMatchedRemoteReads.stream().filter(x -> !x.isSupplementary()).collect(Collectors.toList());
    }

    public AssemblyLink tryRemoteAssemblyLink(
            final JunctionAssembly assembly, final RemoteRegion remoteRegion, final Set<String> sourceReadIds)
    {
        if(assembly.extensionLength() > remoteRegion.length())
            return null;

        mRemoteRegion = remoteRegion;

        mSourceReadIds.clear();
        mSourceReadIds.addAll(sourceReadIds);

        mTotalRemoteReadsSearch += sourceReadIds.size();

        if(mBamReader != null)
        {
            mMatchedRemoteReads.clear();

            SV_LOGGER.trace("remote region({}) slice", mRemoteRegion);

            mBamReader.sliceBam(mRemoteRegion.Chromosome, mRemoteRegion.start(), mRemoteRegion.end(), this::processRecord);

            SV_LOGGER.trace("remote region({}) sourcedReads(matched={} unmatched={})",
                    mRemoteRegion, mMatchedRemoteReads.size(), mSourceReadIds.size());

            mTotalRemoteReadsMatched += mMatchedRemoteReads.size();
        }

        if(mMatchedRemoteReads.size() < ASSEMBLY_MIN_READ_SUPPORT)
            return null;

        // form a remote ref-based assembly from these reads but without a specific junction
        int remoteRegionStart = mMatchedRemoteReads.stream().mapToInt(x -> x.alignmentStart()).min().orElse(0);
        int remoteRegionEnd = mMatchedRemoteReads.stream().mapToInt(x -> x.alignmentEnd()).max().orElse(0);

        // require a region at least as long as the extension sequence
        // if(assembly.extensionLength() > remoteRegionEnd - remoteRegionStart + 1)
        //    return null;

        byte[] refGenomeBases = mRefGenome.getBases(remoteRegion.Chromosome, remoteRegionStart, remoteRegionEnd);

        AssemblyLink assemblyLink = tryAssemblyRemoteRefOverlap(assembly, remoteRegionStart, remoteRegionEnd, refGenomeBases);

        if(assemblyLink == null)
            return null;

        SV_LOGGER.trace("assembly({}) links with remote region({}) matchedReads({})",
                assembly, remoteRegion.toString(), mMatchedRemoteReads.size());

        mSourceReadIds.clear();
        mMatchedRemoteReads.clear();

        return assemblyLink;
    }

    private void processRecord(final SAMRecord record)
    {
        // the read IDs have been trimmed, so has to match on what has been kept
        boolean containedRead = mSourceReadIds.stream().anyMatch(x -> record.getReadName().contains(x));

        if(!containedRead)
            return;

        Read remoteRead = new Read(record);

        if(mBamReader.currentIsReferenceSample())
            remoteRead.markReference();

        mMatchedRemoteReads.add(remoteRead);
    }

    public AssemblyLink tryAssemblyRemoteRefOverlap(
            final JunctionAssembly assembly, final int remoteRegionStart, final int remoteRegionEnd, final byte[] refGenomeBases)
    {
        String remoteRefBases = new String(refGenomeBases);

        String assemblyExtensionBases = assembly.formJunctionSequence();

        // test in orientation in turn
        // JunctionSequence assemblySeq = JunctionSequence.formFullExtensionMatchSequence(assembly, false);

        AssemblyLink assemblyLink = tryAssemblyRemoteRefOverlap(
                assembly, assemblyExtensionBases, false, remoteRegionStart, remoteRegionEnd, remoteRefBases);

        if(assemblyLink != null)
            return assemblyLink;

        assemblyExtensionBases = Nucleotides.reverseComplementBases(assemblyExtensionBases);

        return tryAssemblyRemoteRefOverlap(
                assembly, assemblyExtensionBases, true, remoteRegionStart, remoteRegionEnd, remoteRefBases);
    }

    private AssemblyLink tryAssemblyRemoteRefOverlap(
            final JunctionAssembly assembly, String assemblyExtBases, boolean assemblyReversed,
            final int remoteRegionStart, final int remoteRegionEnd, final String remoteRefBases)
    {
        Orientation remoteOrientation = assembly.isForwardJunction() == (!assemblyReversed) ? REVERSE : FORWARD;

        // start with a simple comparison looking for the first sequence around its junction in the second
        int assemblyExtBaseLength = assemblyExtBases.length();

        // first a simple match where the extension bases are fully contained
        int extBaseIndexStartInRef = remoteRefBases.indexOf(assemblyExtBases);

        if(extBaseIndexStartInRef >= 0)
        {
            return formLinkWithRemote(
                    assembly, assemblyReversed, remoteRefBases, remoteRegionStart, remoteRegionEnd, remoteOrientation,
                    0, extBaseIndexStartInRef);
        }

        // take a smaller sections of the assembly's extension sequence and try to find their start index in the second sequence
        int subsequenceLength = MATCH_SUBSEQUENCE_LENGTH;
        int subSeqIterations = (int)floor(assemblyExtBaseLength / subsequenceLength);

        byte[] assemblyExtBaseQuals = null;
        List<RepeatInfo> assemblyRepeats = null;

        for(int i = 0; i < subSeqIterations; ++i)
        {
            int extBaseSeqStartIndex = i * subsequenceLength;
            int extBaseSeqEndIndex = extBaseSeqStartIndex + subsequenceLength;

            if(extBaseSeqEndIndex > assemblyExtBaseLength)
                break;

            String assemblySubSequence = assemblyExtBases.substring(extBaseSeqStartIndex, extBaseSeqEndIndex);

            int extBaseSubSeqIndexInRef = remoteRefBases.indexOf(assemblySubSequence);

            if(extBaseSubSeqIndexInRef < 0)
                continue;

            // expand the ext base sequence to the start of the remote ref bases and check for a match
            int extBaseStartOffset = extBaseSubSeqIndexInRef - extBaseSeqStartIndex;
            int inferredRegionPosStart = remoteRegionStart + extBaseStartOffset;
            int inferredRegionPosEnd = inferredRegionPosStart + assemblyExtBaseLength - 1;

            String extendedRemoteRefBases = remoteRefBases;
            boolean useAdjustedStart = false;

            if(inferredRegionPosStart < remoteRegionStart)
            {
                int inferredStart = inferredRegionPosStart;
                int inferredEnd = remoteRegionStart - 1;

                String inferredRefBases = mRefGenome.getBaseString(mRemoteRegion.Chromosome, inferredStart, inferredEnd);
                extendedRemoteRefBases = inferredRefBases + remoteRefBases;
                useAdjustedStart = true;
                extBaseIndexStartInRef = extendedRemoteRefBases.indexOf(assemblyExtBases); // exact match test
            }
            else if(inferredRegionPosEnd > remoteRegionEnd)
            {
                int inferredStart = remoteRegionEnd + 1;
                int inferredEnd = inferredRegionPosEnd;

                String inferredRefBases = mRefGenome.getBaseString(mRemoteRegion.Chromosome, inferredStart, inferredEnd);
                extendedRemoteRefBases = remoteRefBases + inferredRefBases;
                extBaseIndexStartInRef = extendedRemoteRefBases.indexOf(assemblyExtBases);
            }

            if(extBaseIndexStartInRef >= 0)
            {
                // exact match using the extended ref bases ie beyond the original remote region bounds
                return formLinkWithRemote(
                        assembly, assemblyReversed, extendedRemoteRefBases,
                        min(inferredRegionPosStart, remoteRegionStart), max(inferredRegionPosEnd, remoteRegionEnd),
                        remoteOrientation, 0, extBaseIndexStartInRef);
            }

            // finally attempt a sequence comparison allowing mismatches - extract the same length ref base sequence to compare
            int remoteRefBaseIndexStart = useAdjustedStart ? 0 : inferredRegionPosStart - remoteRegionStart;

            byte[] remoteRefBaseBytes = extendedRemoteRefBases.substring(
                    remoteRefBaseIndexStart, remoteRefBaseIndexStart + assemblyExtBaseLength).getBytes();

            byte[] remoteRefBaseQuals = createMinBaseQuals(assemblyExtBaseLength);

            if(assemblyExtBaseQuals == null)
            {
                // assemblyRepeats = RepeatInfo.findRepeats(assemblyExtBases.getBytes());
                int extBaseIndexStart, extBaseIndexEnd;

                if(assembly.isForwardJunction())
                {
                    extBaseIndexStart = assembly.junctionIndex() + 1;
                    extBaseIndexEnd = assembly.baseLength() - 1;
                }
                else
                {
                    extBaseIndexStart = 0;
                    extBaseIndexEnd = assembly.junctionIndex() - 1;
                }

                if(assemblyReversed)
                    assemblyExtBaseQuals = reverseArray(subsetArray(assembly.baseQuals(), extBaseIndexStart, extBaseIndexEnd));
                else
                    assemblyExtBaseQuals = subsetArray(assembly.baseQuals(), extBaseIndexStart, extBaseIndexEnd);
            }

            int compareIndexEnd = assemblyExtBaseLength - 1;

            int mismatchCount = SequenceCompare.compareSequences(
                    assemblyExtBases.getBytes(), assemblyExtBaseQuals, 0, compareIndexEnd, Collections.emptyList(),
                    remoteRefBaseBytes, remoteRefBaseQuals, 0, compareIndexEnd, Collections.emptyList(),
                    PRIMARY_ASSEMBLY_MERGE_MISMATCH);

            if(mismatchCount <= PRIMARY_ASSEMBLY_MERGE_MISMATCH)
            {
                // a full match with no inserted bases
                return formLinkWithRemote(
                        assembly, assemblyReversed, extendedRemoteRefBases,
                        min(inferredRegionPosStart, remoteRegionStart), max(inferredRegionPosEnd, remoteRegionEnd),
                        remoteOrientation, 0, remoteRefBaseIndexStart);
            }
        }

        return null;
    }

    private AssemblyLink formLinkWithRemote(
            final JunctionAssembly assembly, boolean assemblyReversed, final String remoteRefBases,
            final int remoteRegionStart, final  int remoteRegionEnd, final Orientation remoteOrientation,
            int extBaseMatchIndexStart, int extBaseIndexStartInRef)
    {
        // work out the remote region implied position where the assembly extension joins, factoring in any inserted bases

        int remoteExtBaseLength = 0; // remains undefined or some overlap from the assembly?

        int assemblyExtBaseLength = assembly.extensionLength();

        int remoteJunctionPosition, remoteRefPosition, remoteJunctionIndex, remoteRefLength, remoteRefBaseIndexStart;

        if(remoteOrientation.isForward())
        {
            remoteJunctionPosition = remoteRegionStart + extBaseIndexStartInRef + assemblyExtBaseLength - 1;
            remoteRefPosition = remoteRegionStart;
            remoteJunctionIndex = remoteJunctionPosition - remoteRefPosition + 1;
            remoteRefLength = remoteJunctionIndex;
            remoteRefBaseIndexStart = 0;
        }
        else
        {
            remoteJunctionPosition = remoteRegionStart + extBaseIndexStartInRef;
            remoteRefPosition = remoteRegionEnd;
            remoteJunctionIndex = 0; // no extension bases for now
            remoteRefLength = remoteRefPosition - remoteJunctionPosition + 1;
            remoteRefBaseIndexStart = extBaseIndexStartInRef;
        }

        Junction remoteJunction = new Junction(mRemoteRegion.Chromosome, remoteJunctionPosition, remoteOrientation);

        byte[] remoteRefBaseQuals = createLowBaseQuals(remoteRefLength); // since not supported by reads yet
        byte[] remoteRefTrimmedBases = remoteRefBases.substring(remoteRefBaseIndexStart, remoteRefBaseIndexStart + remoteRefLength).getBytes();

        List<SupportRead> remoteSupport = Lists.newArrayList();

        for(Read read : mMatchedRemoteReads)
        {
            boolean spansJunction = positionWithin(remoteJunctionPosition, read.unclippedStart(), read.unclippedEnd());

            int matchLength = read.basesLength();

            int junctionReadStartDistance = remoteJunctionPosition - read.unclippedStart();

            SupportRead support = new SupportRead(
                    read, spansJunction ? SupportType.JUNCTION : SupportType.DISCORDANT, junctionReadStartDistance, matchLength, 0);

            remoteSupport.add(support);
        }

        JunctionAssembly remoteAssembly = new JunctionAssembly(
                remoteJunction, remoteRefTrimmedBases, remoteRefBaseQuals, remoteSupport, Lists.newArrayList());

        remoteAssembly.setJunctionIndex(remoteJunctionIndex);

        remoteAssembly.buildRepeatInfo();

        // switch if first is -ve orientation as per normal link testing
        if(assemblyReversed)
            return new AssemblyLink(remoteAssembly, assembly, LinkType.SPLIT, "", "");
        else
            return new AssemblyLink(assembly, remoteAssembly, LinkType.SPLIT, "", "");
    }

    private boolean hasSpanningReads(int remotePosition)
    {
        return mMatchedRemoteReads.stream()
                .anyMatch(x -> positionWithin(remotePosition, x.unclippedStart(), x.unclippedEnd()));
    }

    public AssemblyLink tryAssemblyRemoteRefOverlapOld(
            final JunctionAssembly assembly, final int remoteRegionStart, final int remoteRegionEnd, final byte[] refGenomeBases)
    {
        AssemblyLink assemblyLink = tryAssemblyRemoteRefOverlapOld(assembly, remoteRegionStart, remoteRegionEnd, refGenomeBases, REVERSE);

        if(assemblyLink == null)
            assemblyLink = tryAssemblyRemoteRefOverlapOld(assembly, remoteRegionStart, remoteRegionEnd, refGenomeBases, FORWARD);

        return assemblyLink;
    }

    private AssemblyLink tryAssemblyRemoteRefOverlapOld(
            final JunctionAssembly assembly, final int remoteRegionStart, final int remoteRegionEnd, final byte[] refGenomeBases,
            final Orientation remoteOrientation)
    {
        byte[] refBaseQuals = createMinBaseQuals(refGenomeBases.length);

        boolean assemblyReversed = false;
        boolean remoteReversed = false;

        if(assembly.junction().Orient == remoteOrientation)
        {
            if(assembly.junction().isForward())
                remoteReversed = true;
            else
                assemblyReversed = true;
        }

        JunctionSequence assemblySeq = JunctionSequence.formFullExtensionMatchSequence(assembly, assemblyReversed);

        JunctionSequence remoteRefSeq = new JunctionSequence(refGenomeBases, refBaseQuals, remoteOrientation, remoteReversed);

        // start with a simple comparison looking for the first sequence around its junction in the second
        String firstMatchSequence = assemblySeq.matchSequence();
        int firstMatchSeqLength = firstMatchSequence.length();

        // first a simple local match
        int remoteSeqIndexInRef = remoteRefSeq.FullSequence.indexOf(firstMatchSequence);

        if(remoteSeqIndexInRef >= 0)
        {
            return formLinkWithRemoteOld(
                    assembly, assemblySeq, remoteRefSeq, refGenomeBases, remoteRegionStart, remoteRegionEnd, remoteOrientation,
                    assemblySeq.matchSeqStartIndex(), remoteSeqIndexInRef);
        }

        // take a smaller sections of the first's junction sequence and try to find their start index in the second sequence
        int matchSeqStartIndex = 0;
        List<int[]> alternativeIndexStarts = Lists.newArrayList();

        int subsequenceLength = MATCH_SUBSEQUENCE_LENGTH;
        int subSeqIterations = (int)floor(firstMatchSeqLength / subsequenceLength);

        for(int i = 0; i < subSeqIterations; ++i)
        {
            matchSeqStartIndex = i * subsequenceLength;
            int matchSeqEndIndex = matchSeqStartIndex + subsequenceLength;

            if(matchSeqEndIndex >= firstMatchSeqLength)
                break;

            String firstSubSequence = firstMatchSequence.substring(matchSeqStartIndex, matchSeqStartIndex + subsequenceLength);

            int secondSubSeqIndex = remoteRefSeq.FullSequence.indexOf(firstSubSequence);

            if(secondSubSeqIndex < 0)
                continue;

            alternativeIndexStarts.add(new int[] {matchSeqStartIndex, secondSubSeqIndex});

            secondSubSeqIndex = remoteRefSeq.FullSequence.indexOf(firstSubSequence, secondSubSeqIndex + subsequenceLength);

            while(secondSubSeqIndex >= 0)
            {
                alternativeIndexStarts.add(new int[] {matchSeqStartIndex, secondSubSeqIndex});
                secondSubSeqIndex = remoteRefSeq.FullSequence.indexOf(firstSubSequence, secondSubSeqIndex + subsequenceLength);
            }
        }

        if(alternativeIndexStarts.isEmpty())
            return null;

        // now perform a full junction sequence search in the second using the sequence matching logic
        int minOverlapLength = min(assembly.extensionLength(), ASSEMBLY_LINK_OVERLAP_BASES);

        int[] topMatchIndices = findBestSequenceMatch(assemblySeq, remoteRefSeq, minOverlapLength, alternativeIndexStarts);

        if(topMatchIndices != null)
        {
            int firstIndexStart = topMatchIndices[0];
            int secondIndexStart = topMatchIndices[1];

            // now that the index in the remote ref sequence has a match and it is clear where this is in the assembly's extension sequence,
            // the implied junction position in the remote can be determined
            return formLinkWithRemoteOld(
                    assembly, assemblySeq, remoteRefSeq, refGenomeBases, remoteRegionStart, remoteRegionEnd, remoteOrientation,
                    firstIndexStart, secondIndexStart);
        }

        return null;
    }

    private AssemblyLink formLinkWithRemoteOld(
            final JunctionAssembly assembly, final JunctionSequence assemblySeq, final JunctionSequence initialRemoteRefSeq,
            final byte[] refGenomeBases, final int remoteRegionStart, final  int remoteRegionEnd, final Orientation remoteOrientation,
            int firstMatchSeqIndexStart, int remoteMatchIndexStart)
    {
        // use the start position of the match to infer where the junction may be in the remote location despite there not being any junction
        // spanning reads at that position
        int assemblyMatchJunctionOffset = firstMatchSeqIndexStart - assemblySeq.matchSeqStartIndex();

        int remoteJunctionPosition, remoteJunctionIndex;

        int inferredRemoteJunctionPosition;

        if(remoteOrientation.isForward())
        {
            remoteJunctionPosition = remoteRegionEnd;
            remoteJunctionIndex = refGenomeBases.length - 1;
            inferredRemoteJunctionPosition = remoteRegionEnd + assemblyMatchJunctionOffset;
        }
        else
        {
            remoteJunctionPosition = remoteRegionStart + remoteMatchIndexStart;
            remoteJunctionIndex = 0;
            inferredRemoteJunctionPosition = remoteRegionStart - assemblyMatchJunctionOffset + remoteMatchIndexStart;
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

            if(remoteOrientation.isForward())
            {
                inferredStart = remoteRegionEnd;
                inferredEnd = remoteRegionEnd + assemblyMatchJunctionOffset;

            }
            else
            {
                inferredStart = remoteRegionStart - assemblyMatchJunctionOffset;
                inferredEnd = remoteRegionStart;
            }

            String inferredRefBases = mRefGenome.getBaseString(mRemoteRegion.Chromosome, inferredStart, inferredEnd);

            // check for a simple match at the assembly's junction
            String assemblyExtBases = assemblySeq.matchSequence();

            int secondIndexInFirst = assemblyExtBases.indexOf(inferredRefBases);

            if(secondIndexInFirst >= 0)
            {
                int adjustedRemoteStart = remoteOrientation.isForward() ? remoteRegionStart : inferredRemoteJunctionPosition;
                int adjustedRemoteEnd = remoteOrientation.isForward() ? inferredRemoteJunctionPosition : remoteRegionEnd;

                byte[] remoteRefBases = mRefGenome.getBases(mRemoteRegion.Chromosome, adjustedRemoteStart, adjustedRemoteEnd);
                byte[] remoteRefBaseQuals = createMinBaseQuals(remoteRefBases.length);

                remoteRefSeq = new JunctionSequence(remoteRefBases, remoteRefBaseQuals, remoteOrientation, initialRemoteRefSeq.Reversed);

                if(remoteOrientation.isForward())
                    remoteJunctionIndex = remoteRefBases.length - 1;;
            }
        }
        else
        {
            // suggests that the break is further into the initially selected ref bases, ie they need to be truncated
            // but those bases are still valid extension bases
            // remoteJunctionIndex needs adjusting?
        }

        List<SupportRead> remoteSupport = Lists.newArrayList();

        for(Read read : mMatchedRemoteReads)
        {
            boolean spansJunction = positionWithin(remoteJunctionPosition, read.unclippedStart(), read.unclippedEnd());

            int matchLength = read.basesLength();

            int junctionReadStartDistance = remoteJunctionPosition - read.unclippedStart();

            SupportRead support = new SupportRead(
                    read, spansJunction ? SupportType.JUNCTION : SupportType.DISCORDANT, junctionReadStartDistance, matchLength, 0);

            remoteSupport.add(support);
        }

        // TODO: consider adjusting the start pos and sequence to the inferred remote junction, or let reads from other remote locations
        // fill in this gap?

        Junction remoteJunction = new Junction(mRemoteRegion.Chromosome, remoteJunctionPosition, remoteOrientation);

        JunctionAssembly remoteAssembly = new JunctionAssembly(
                remoteJunction, remoteRefSeq.originalBases(), remoteRefSeq.originalBaseQuals(), remoteSupport, Lists.newArrayList());

        remoteAssembly.setJunctionIndex(remoteJunctionIndex);

        remoteAssembly.buildRepeatInfo();

        // switch if first is -ve orientation as per normal link testing
        if(assemblySeq.Reversed)
            return new AssemblyLink(remoteAssembly, assembly, LinkType.SPLIT, "", "");
        else
            return new AssemblyLink(assembly, remoteAssembly, LinkType.SPLIT, "", "");
    }

    @VisibleForTesting
    public void addMatchedReads(final List<Read> reads, final RemoteRegion remoteRegion)
    {
        mRemoteRegion = remoteRegion;
        mMatchedRemoteReads.clear();
        mMatchedRemoteReads.addAll(reads);
    }
}
