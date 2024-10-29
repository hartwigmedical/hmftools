package com.hartwig.hmftools.esvee.assembly.phase;

import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.Arrays.reverseArray;
import static com.hartwig.hmftools.common.utils.Arrays.subsetArray;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MATCH_SUBSEQUENCE_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PRIMARY_ASSEMBLY_MERGE_MISMATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.calcTrimmedExtensionBaseLength;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.createLowBaseQuals;
import static com.hartwig.hmftools.esvee.assembly.AssemblyUtils.createMinBaseQuals;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.Collections;
import java.util.Comparator;
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
    private int mRemoteReadSlices;

    public RemoteRegionAssembler(final RefGenomeInterface refGenome, final BamReader bamReader)
    {
        mRefGenome = refGenome;
        mBamReader = bamReader;

        mRemoteRegion = null;
        mSourceReadIds = Sets.newHashSet();
        mMatchedRemoteReads = Lists.newArrayList();

        mTotalRemoteReadsSearch = 0;
        mTotalRemoteReadsMatched = 0;
        mRemoteReadSlices = 0;
    }

    public int remoteReadSlices() { return mRemoteReadSlices; }
    public int remoteReadsSearch() { return mTotalRemoteReadsSearch; }
    public int remoteReadsMatched() { return mTotalRemoteReadsMatched; }

    public static boolean isExtensionCandidateAssembly(final JunctionAssembly assembly)
    {
        // apply some filters to limit the number of assemblies which attempt to find a remote discordant match
        if(assembly.refBaseTrimLength() < MIN_VARIANT_LENGTH)
            return false;

        if(assembly.stats().SoftClipSecondMaxLength < MIN_VARIANT_LENGTH)
            return false;

        // check for sufficient diversity in the extension bases
        int trimmedExtBaseLength = calcTrimmedExtensionBaseLength(assembly);

        if(trimmedExtBaseLength < MIN_VARIANT_LENGTH)
            return false;

        return true;
    }

    public static List<RemoteRegion> collectCandidateRemoteRegions(
            final JunctionAssembly assembly, final List<JunctionAssembly> phasedAssemblies, boolean checkAssemblyMatches)
    {
        List<RemoteRegion> combinedRemoteRegions = Lists.newArrayList();

        // ignore remote regions which overlap an assembly in this phase group with any matching reads
        for(RemoteRegion remoteRegion : assembly.remoteRegions())
        {
            if(remoteRegion.isSuppOnlyRegion())
                continue;

            if(!checkAssemblyMatches)
            {
                combinedRemoteRegions.add(remoteRegion);
                continue;
            }

            boolean matchesAssembly = false;

            for(JunctionAssembly otherAssembly : phasedAssemblies)
            {
                if(otherAssembly == assembly)
                    continue;

                if(otherAssembly.discordantOnly()) // allow overlaps with discordant junctions since these don't have precise junctions
                    continue;

                if(!remoteRegion.overlapsAssembly(otherAssembly))
                    continue;

                // otherwise ignore if the overlapping assembly has these remote reads already
                // note the max reads checked to avoid checking all reads in very large groups
                int readsChecked = 0;
                for(String readId : remoteRegion.readIds())
                {
                    if(otherAssembly.support().stream().filter(x -> x.type().isSplitSupport()).anyMatch(x -> x.id().equals(readId)))
                    {
                        matchesAssembly = true;
                        break;
                    }

                    if(otherAssembly.candidateSupport().stream().anyMatch(x -> x.id().equals(readId)))
                    {
                        matchesAssembly = true;
                        break;
                    }

                    ++readsChecked;
                    if(readsChecked > MAX_MATCHED_READ_CHECK)
                        break;
                }
            }

            if(!matchesAssembly)
                combinedRemoteRegions.add(remoteRegion);
        }

        return combinedRemoteRegions;
    }

    private static final int HIGH_REMOTE_REGION_TOTAL_READ_COUNT = 1000;
    private static final int HIGH_REMOTE_REGION_READ_COUNT = 200;
    private static final int MAX_MATCHED_READ_CHECK = 100;

    public void extractRemoteRegionReads(
            final int phaseGroupId, final List<RemoteRegion> remoteRegions, final List<Read> remoteRegionReads, boolean applyThresholds)
    {
        // slices the BAM at the specified locations for mates of reads in one or more assemblies
        // first merges overlapping reads to avoid repeated reads
        RemoteRegion.mergeRegions(remoteRegions);

        int minRemoteRegionReadCount = 1;
        int maxRemoteRegionReadCount = 0;

        int totalRemoteReads = remoteRegions.stream().mapToInt(x -> x.readCount()).sum();

        // impose restrictions on which remote regions are used - excluding those with 1 or very few reads relatively, eg for 1000 then
        // a region will be required to have 2+ reads
        if(applyThresholds && totalRemoteReads >= HIGH_REMOTE_REGION_TOTAL_READ_COUNT)
        {
            minRemoteRegionReadCount = (int)floor(log10(totalRemoteReads));
            Collections.sort(remoteRegions, Comparator.comparingInt(x -> -x.readCount()));

            maxRemoteRegionReadCount = HIGH_REMOTE_REGION_TOTAL_READ_COUNT;

            int cappedTotalRemoteReads = remoteRegions.stream().mapToInt(x -> min(x.readCount(), HIGH_REMOTE_REGION_READ_COUNT)).sum();

            if(cappedTotalRemoteReads > HIGH_REMOTE_REGION_TOTAL_READ_COUNT)
            {
                double cappedFraction = HIGH_REMOTE_REGION_TOTAL_READ_COUNT/(double)cappedTotalRemoteReads;
                maxRemoteRegionReadCount = (int)round(cappedFraction * maxRemoteRegionReadCount);
            }
        }

        for(RemoteRegion remoteRegion : remoteRegions)
        {
            if(minRemoteRegionReadCount > 1 && remoteRegion.readCount() < minRemoteRegionReadCount)
                break;

            List<Read> remoteReads = extractRemoteReads(phaseGroupId, remoteRegion);

            if(maxRemoteRegionReadCount > 0 && remoteReads.size() > maxRemoteRegionReadCount)
            {
                remoteReads = remoteReads.subList(0, maxRemoteRegionReadCount);
            }

            // finally as a way of down-sample, restrict a single region's read contribution
            remoteRegionReads.addAll(remoteReads);
        }
    }

    public List<Read> extractRemoteReads(final int phaseGroupId, final RemoteRegion remoteRegion)
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

            // SV_LOGGER.trace("remote region({}) sourcedReads(matched={} unmatched={})",
            //        mRemoteRegion, mMatchedRemoteReads.size(), mSourceReadIds.size());

            mTotalRemoteReadsMatched += mMatchedRemoteReads.size();
        }

        ++mRemoteReadSlices;

        if((mRemoteReadSlices % 1000) == 0 && mTotalRemoteReadsMatched > 100000)
        {
            SV_LOGGER.debug("pgId({}) remote region read extraction: slices({}) matched({}) last region({}:{}-{})",
                    phaseGroupId, mRemoteReadSlices, mTotalRemoteReadsMatched,
                    mRemoteRegion.Chromosome, mRemoteRegion.start(), mRemoteRegion.end());
        }

        // ignore supplementaries since their bases provide no new assembly sequence information
        return mMatchedRemoteReads.stream().filter(x -> !x.isSupplementary()).collect(Collectors.toList());
    }

    public AssemblyLink tryRemoteAssemblyLink(
            final JunctionAssembly assembly, final RemoteRegion remoteRegion, final Set<String> sourceReadIds)
    {
        // require a region at least as long as the extension sequence
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

        // form a reference sequence from these reads but without a specific junction
        int remoteRegionStart = -1;
        int remoteRegionEnd = -1;

        for(Read read : mMatchedRemoteReads)
        {
            if(remoteRegionStart < 0 || read.alignmentStart() < remoteRegionStart)
                remoteRegionStart = read.alignmentStart();

            remoteRegionEnd = max(remoteRegionEnd, read.alignmentEnd());
        }

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
        List<RepeatInfo> assemblyRepeats = null; // would establish these to aid with matching but for not only expect SNVs

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
            int remoteRefBaseIndexEnd = remoteRefBaseIndexStart + assemblyExtBaseLength;

            if(remoteRefBaseIndexStart >= extendedRemoteRefBases.length()|| remoteRefBaseIndexEnd > extendedRemoteRefBases.length())
                return null;

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
            remoteRefLength = remoteJunctionPosition - remoteRefPosition + 1;
            remoteJunctionIndex = remoteRefLength - 1;
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

    @VisibleForTesting
    public void addMatchedReads(final List<Read> reads, final RemoteRegion remoteRegion)
    {
        mRemoteRegion = remoteRegion;
        mMatchedRemoteReads.clear();
        mMatchedRemoteReads.addAll(reads);
    }
}
