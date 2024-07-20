package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.floor;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.AssemblyConstants.ASSEMBLY_LINK_OVERLAP_BASES;
import static com.hartwig.hmftools.esvee.AssemblyConstants.MATCH_SUBSEQUENCE_LENGTH;
import static com.hartwig.hmftools.esvee.AssemblyConstants.PRIMARY_ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.calcTrimmedExtensionBaseLength;
import static com.hartwig.hmftools.esvee.assembly.phase.AssemblyLinker.findBestSequenceMatch;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.createMinBaseQuals;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.read.BamReader;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
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

        if(assembly.stats().JuncMateDiscordantRemote < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return false;

        // check for sufficient diversity in the extension bases
        int trimmedExtBaseLength = calcTrimmedExtensionBaseLength(assembly);

        if(trimmedExtBaseLength < MIN_VARIANT_LENGTH)
            return false;

        return true;
    }

    public static boolean assemblyOverlapsRemoteRegion(final JunctionAssembly assembly, final RemoteRegion remoteRegion)
    {
        if(assembly.isForwardJunction())
            return positionsOverlap(assembly.refBasePosition(), assembly.junction().Position, remoteRegion.start(), remoteRegion.end());
        else
            return positionsOverlap(assembly.junction().Position, assembly.refBasePosition(), remoteRegion.start(), remoteRegion.end());
    }

    public AssemblyLink tryRemoteAssemblyLink(
            final JunctionAssembly assembly, final RemoteRegion remoteRegion, final List<String> sourceReadIds)
    {
        mRemoteRegion = remoteRegion;

        SV_LOGGER.trace("remote region({}) slice", mRemoteRegion);

        mMatchedRemoteReads.clear();
        mSourceReadIds.clear();

        mTotalRemoteReadsSearch += sourceReadIds.size();

        mSourceReadIds.addAll(sourceReadIds);

        mBamReader.sliceBam(mRemoteRegion.Chromosome, mRemoteRegion.start(), mRemoteRegion.end(), this::processRecord);

        SV_LOGGER.trace("remote region({}) sourcedReads(matched={} unmatched={})",
                mRemoteRegion, mMatchedRemoteReads.size(), mSourceReadIds.size());

        mTotalRemoteReadsMatched += mMatchedRemoteReads.size();

        if(mMatchedRemoteReads.size() < PRIMARY_ASSEMBLY_MIN_READ_SUPPORT)
            return null;

        // form a remote ref-based assembly from these reads but without a specific junction
        int remoteRegionStart = mMatchedRemoteReads.stream().mapToInt(x -> x.alignmentStart()).min().orElse(0);
        int remoteRegionEnd = mMatchedRemoteReads.stream().mapToInt(x -> x.alignmentEnd()).max().orElse(0);

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
        byte[] refBaseQuals = createMinBaseQuals(refGenomeBases.length);

        boolean assemblyReversed = false;
        boolean remoteReversed = false;

        if(assembly.junction().Orient == mRemoteRegion.orientation())
        {
            if(assembly.junction().isForward())
                remoteReversed = true;
            else
                assemblyReversed = true;
        }

        JunctionSequence assemblySeq = new JunctionSequence(assembly, assemblyReversed, 0, -1);

        JunctionSequence remoteRefSeq = new JunctionSequence(refGenomeBases, refBaseQuals, mRemoteRegion.orientation(), remoteReversed);

        // start with a simple comparison looking for the first sequence around its junction in the second
        String firstJunctionSequence = assemblySeq.junctionSequence();
        int firstJunctionSeqLength = firstJunctionSequence.length();

        // first a simple local match
        int remoteSeqIndexInRef = remoteRefSeq.FullSequence.indexOf(firstJunctionSequence);

        if(remoteSeqIndexInRef >= 0)
        {
            return formLinkWithRemote(
                    assembly, assemblySeq, remoteRefSeq, refGenomeBases, remoteRegionStart, remoteRegionEnd,
                    assemblySeq.junctionSeqStartIndex(), remoteSeqIndexInRef);
        }

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
        int minOverlapLength = min(assembly.extensionLength(), ASSEMBLY_LINK_OVERLAP_BASES);

        int[] topMatchIndices = findBestSequenceMatch(assemblySeq, remoteRefSeq, minOverlapLength, alternativeIndexStarts);

        if(topMatchIndices != null)
        {
            int firstIndexStart = topMatchIndices[0];
            int secondIndexStart = topMatchIndices[1];

            // now that the index in the remote ref sequence has a match and it is clear where this is in the assembly's extension sequence,
            // the implied junction position in the remote can be determined
            return formLinkWithRemote(
                    assembly, assemblySeq, remoteRefSeq, refGenomeBases, remoteRegionStart, remoteRegionEnd,
                    firstIndexStart, secondIndexStart);
        }

        return null;
    }

    private AssemblyLink formLinkWithRemote(
            final JunctionAssembly assembly, final JunctionSequence assemblySeq, final JunctionSequence initialRemoteRefSeq,
            final byte[] refGenomeBases, final int remoteRegionStart, final  int remoteRegionEnd,
            int firstJuncSeqIndexStart, int secondIndexStart)
    {
        // use the start position of the match to infer where the junction may be in the remote location despite there not being any junction
        // spanning reads at that position
        int assemblyMatchJunctionOffset = firstJuncSeqIndexStart - assemblySeq.junctionSeqStartIndex();

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

            String inferredRefBases = mRefGenome.getBaseString(mRemoteRegion.Chromosome, inferredStart, inferredEnd);

            // check for a simple match at the assembly's junction
            String assemblyJunctionSequence = assemblySeq.junctionSequence();

            int secondIndexInFirst = assemblyJunctionSequence.indexOf(inferredRefBases);

            if(secondIndexInFirst >= 0)
            {
                firstJuncSeqIndexStart = assemblySeq.junctionSeqStartIndex();

                if(!mRemoteRegion.isForward())
                    secondIndexStart = 0;

                int adjustedRemoteStart = mRemoteRegion.isForward() ? remoteRegionStart : inferredRemoteJunctionPosition;
                int adjustedRemoteEnd = mRemoteRegion.isForward() ? inferredRemoteJunctionPosition : remoteRegionEnd;

                byte[] remoteRefBases = mRefGenome.getBases(mRemoteRegion.Chromosome, adjustedRemoteStart, adjustedRemoteEnd);
                byte[] remoteRefBaseQuals = createMinBaseQuals(remoteRefBases.length);

                remoteRefSeq = new JunctionSequence(remoteRefBases, remoteRefBaseQuals, mRemoteRegion.orientation(), initialRemoteRefSeq.Reversed);

                if(mRemoteRegion.isForward())
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

            SupportRead support = new SupportRead(
                    read, spansJunction ? SupportType.JUNCTION : SupportType.DISCORDANT, 0, matchLength, 0);

            remoteSupport.add(support);
        }

        // TODO: consider adjusting the start pos and sequence to the inferred remote junction, or let reads from other remote locations
        // fill in this gap?

        Junction remoteJunction = new Junction(mRemoteRegion.Chromosome, remoteJunctionPosition, mRemoteRegion.orientation());

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
