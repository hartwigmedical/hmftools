package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_DISTINCT_FRAGS;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MAX_OBSERVED_CONCORDANT_FRAG_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findIndelExtensionReads;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.hasIndelJunctionReads;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.isLineWithLocalAlignedInsert;
import static com.hartwig.hmftools.esvee.assembly.RemoteRegionFinder.addOrCreateMateRemoteRegion;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.readJunctionExtensionLength;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.readSoftClipsAndCrossesJunction;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.recordSoftClipsAtJunction;
import static com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip.checkSupportVsRefSideSoftClip;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteRegion.mergeRegions;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_EXTENSION_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_HOTSPOT_JUNCTION_SUPPORT;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.read.Read;

public class JunctionAssembler
{
    private Junction mJunction;
    private final RefGenomeInterface mRefGenome;
    private final List<Read> mNonJunctionReads;

    public JunctionAssembler(final Junction junction, final RefGenomeInterface refGenome)
    {
        mJunction = junction;
        mRefGenome = refGenome;
        mNonJunctionReads = Lists.newArrayList();
    }

    public List<Read> nonJunctionReads() { return mNonJunctionReads; }

    public List<JunctionAssembly> processJunction(final List<Read> rawReads)
    {
        // find prominent reads to establish the extension sequence, taking any read meeting min soft-clip lengths
        // and repetitive indels

        List<Read> junctionReads = Lists.newArrayList();
        List<Read> extensionReads = Lists.newArrayList();

        if(!mJunction.indelBased() && hasIndelJunctionReads(mJunction, rawReads))
        {
            // fall-back in case Prep didn't set this state or junctions are loaded from config
            mJunction.markAsIndel();
        }

        boolean hasMinLengthSoftClipRead = false;

        if(mJunction.indelBased())
        {
            findIndelExtensionReads(mJunction, rawReads, extensionReads, junctionReads, mNonJunctionReads);
            hasMinLengthSoftClipRead = !extensionReads.isEmpty();
        }
        else if(mJunction.DiscordantOnly)
        {
            // look for a common soft-clip position, otherwise take the min variant length back from the inner most read as the junction
            assessDiscordantJunction(rawReads, extensionReads, junctionReads);
            hasMinLengthSoftClipRead = !extensionReads.isEmpty();
        }
        else
        {
            // the only difference for indel-based junctions is that only the long indels are used to build the consensus extension
            for(Read read : rawReads)
            {
                if(!readSoftClipsAndCrossesJunction(read, mJunction, mRefGenome))
                {
                    mNonJunctionReads.add(read);
                    continue;
                }

                if(recordSoftClipsAtJunction(read, mJunction))
                {
                    int softClipJunctionExtension = readJunctionExtensionLength(read, mJunction);

                    if(read.hasLineTail())
                    {
                        hasMinLengthSoftClipRead |= softClipJunctionExtension >= LINE_MIN_EXTENSION_LENGTH;

                        if(softClipJunctionExtension >= LINE_MIN_SOFT_CLIP_SECONDARY_LENGTH)
                            extensionReads.add(read);
                    }
                    else
                    {

                        hasMinLengthSoftClipRead |= softClipJunctionExtension >= ASSEMBLY_MIN_SOFT_CLIP_LENGTH;

                        if(softClipJunctionExtension >= ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH)
                            extensionReads.add(read);
                    }
                }

                junctionReads.add(read);
            }
        }

        if(!hasMinLengthSoftClipRead || !aboveMinReadThreshold(extensionReads))
            return Collections.emptyList();

        ExtensionSeqBuilderOld extensionSeqBuilder = new ExtensionSeqBuilderOld(mJunction, extensionReads);

        int reqExtensionLength = extensionSeqBuilder.hasLineSequence() ? LINE_MIN_EXTENSION_LENGTH : ASSEMBLY_MIN_SOFT_CLIP_LENGTH;

        if(!extensionSeqBuilder.isValid() || extensionSeqBuilder.extensionLength() < reqExtensionLength)
            return Collections.emptyList();

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        if(!meetsMinSupportThreshold(assemblySupport))
            return Collections.emptyList();

        JunctionAssembly firstAssembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualities(), assemblySupport,
                extensionSeqBuilder.repeatInfo());

        // filter LINE source-type sites marked by opposition orientation poly A/T sequences
        if(!firstAssembly.indel() && LineUtils.hasLineSourceSequence(firstAssembly))
            return Collections.emptyList();

        firstAssembly.setExtBaseBuildInfo(extensionSeqBuilder.buildInformation());

        if(extensionSeqBuilder.hasLineSequence())
            firstAssembly.markLineSequence();

        List<JunctionAssembly> assemblies = Lists.newArrayList(firstAssembly);

        addJunctionReads(firstAssembly, extensionSeqBuilder, junctionReads);

        if(firstAssembly.hasLineSequence())
        {
            if(isLineWithLocalAlignedInsert(firstAssembly))
            {
                firstAssembly.unmarkLineSequence();

                if(firstAssembly.extensionLength() < ASSEMBLY_MIN_SOFT_CLIP_LENGTH)
                    return Collections.emptyList();
            }
        }

        // test for a second well-supported, alternative assembly at the same junction
        JunctionAssembly secondAssembly = checkSecondAssembly(extensionSeqBuilder.mismatchReads(), firstAssembly, junctionReads);

        if(secondAssembly != null)
            assemblies.add(secondAssembly);

        for(JunctionAssembly assembly : assemblies)
        {
            // call routine to purge reads likely add to multiple assemblies and breaching a dominant RSSC
            // NOTE: this could be done for all assemblies, not just split ones
            if(assemblies.size() > 1)
                checkSupportVsRefSideSoftClip(assembly);

            RefBaseSeqBuilder refBaseSeqBuilder = new RefBaseSeqBuilder(assembly);
            assembly.setRefBases(refBaseSeqBuilder);

            assembly.buildRepeatInfo();
        }

        return assemblies;
    }

    private void assessDiscordantJunction(final List<Read> rawReads, final List<Read> extensionReads, final List<Read> junctionReads)
    {
        if(rawReads.size() < ASSEMBLY_MIN_READ_SUPPORT)
            return;

        int originalJuncPosition = mJunction.Position;

        // find the applicable remote region from amongst
        List<Read> candidateJunctionReads = rawReads.stream()
                .filter(x -> (mJunction.isForward() && x.alignmentEnd() == originalJuncPosition)
                        || (mJunction.isReverse() && x.alignmentStart() == originalJuncPosition))
                .collect(Collectors.toList());

        if(candidateJunctionReads.isEmpty())
            return;

        List<RemoteRegion> remoteRegions = Lists.newArrayList();
        candidateJunctionReads.forEach(x -> addOrCreateMateRemoteRegion(remoteRegions, x, true));

        if(remoteRegions.isEmpty())
            return;

        if(remoteRegions.size() > 1)
        {
            mergeRegions(remoteRegions);
            Collections.sort(remoteRegions, Comparator.comparingInt(x -> -x.readCount()));
        }

        RemoteRegion mainRemoteRegion = remoteRegions.get(0);

        // take the inner-most read position and its remote region, then only use other reads with matching remote regions
        List<Integer> extensionJuncPositions = Lists.newArrayListWithCapacity(rawReads.size());

        Set<String> candidateReadIds = Sets.newHashSet();

        int minReadPosStart, maxReadPosEnd;

        if(mJunction.isForward())
        {
            minReadPosStart = max(1, originalJuncPosition - MAX_OBSERVED_CONCORDANT_FRAG_LENGTH);
            maxReadPosEnd = originalJuncPosition;
        }
        else
        {
            minReadPosStart = originalJuncPosition;
            maxReadPosEnd = originalJuncPosition + MAX_OBSERVED_CONCORDANT_FRAG_LENGTH;
        }

        for(Read read : rawReads)
        {
            if(read.mappingQuality() < ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY)
                continue;

            // check that the read maps to the same remote region to qualify as an junction candidate
            if(!mainRemoteRegion.overlaps(read.mateChromosome(), read.mateAlignmentStart(), read.mateAlignmentEnd()))
                continue;

            // ensure reads aren't past the inner-most discordant read, nor outside valid bounds
            if(read.alignmentStart() < minReadPosStart || read.alignmentEnd() > maxReadPosEnd)
                continue;

            candidateReadIds.add(read.id());

            if(mJunction.isForward())
            {
                if(read.alignmentEnd() > originalJuncPosition)
                    continue;

                extensionJuncPositions.add(read.unclippedEnd());
            }
            else
            {
                if(read.alignmentStart() < originalJuncPosition)
                    continue;

                extensionJuncPositions.add(read.unclippedStart());
            }
        }

        if(extensionJuncPositions.size() < ASSEMBLY_MIN_READ_SUPPORT)
            return;

        if(mJunction.isForward())
            Collections.sort(extensionJuncPositions, Comparator.reverseOrder());
        else
            Collections.sort(extensionJuncPositions);

        int outerJuncPosition = extensionJuncPositions.get(0);
        int secondJuncPosition = extensionJuncPositions.get(ASSEMBLY_MIN_READ_SUPPORT - 1);
        int adjustedJuncPosition;

        if(mJunction.isForward())
        {
            outerJuncPosition -= ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
            secondJuncPosition -= ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
            adjustedJuncPosition = min(outerJuncPosition, secondJuncPosition);
        }
        else
        {
            outerJuncPosition += ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
            secondJuncPosition += ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
            adjustedJuncPosition = max(outerJuncPosition, secondJuncPosition);
        }

        mJunction = new Junction(
                mJunction.Chromosome, adjustedJuncPosition, mJunction.Orient, true, false, false);

        mJunction.setRawDiscordantPosition(originalJuncPosition);

        for(Read read : rawReads)
        {
            if(read.alignmentStart() < minReadPosStart || read.alignmentEnd() > maxReadPosEnd)
                continue;

            // skip indels supporting discordant-only junctions
            if(mJunction.isForward())
            {
                if(read.hasIndelImpliedUnclippedEnd())
                    continue;
            }
            else
            {
                if(read.hasIndelImpliedUnclippedStart())
                    continue;
            }

            if(!candidateReadIds.contains(read.id()) || read.mappingQuality() < ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY)
            {
                mNonJunctionReads.add(read);
                continue;
            }

            int extensionLength;

            if(mJunction.isForward())
            {
                if(read.alignmentEnd() > originalJuncPosition)
                    continue;

                extensionLength = read.alignmentEnd() >= adjustedJuncPosition ? read.unclippedEnd() - adjustedJuncPosition : -1;
            }
            else
            {
                if(read.alignmentStart() < originalJuncPosition)
                    continue;

                extensionLength = read.alignmentStart() <= adjustedJuncPosition ? adjustedJuncPosition - read.unclippedStart() : -1;
            }

            if(extensionLength >= ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH)
                extensionReads.add(read);
            else if(extensionLength > 0)
                junctionReads.add(read);
            else
                mNonJunctionReads.add(read);
        }
    }

    private JunctionAssembly checkSecondAssembly(
            final List<Read> extensionReads, final JunctionAssembly firstAssembly, final List<Read> junctionReads)
    {
        if(extensionReads.isEmpty() || mJunction.DiscordantOnly)
            return null;

        if(firstAssembly.hasLineSequence())
            return null;

        ExtensionSeqBuilderOld extensionSeqBuilder = new ExtensionSeqBuilderOld(mJunction, extensionReads);

        if(!extensionSeqBuilder.isValid())
            return null;

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        int secondSupport = assemblySupport.size();
        double secondSupportPerc = secondSupport / (double)firstAssembly.supportCount();

        if(secondSupport < ASSEMBLY_SPLIT_MIN_READ_SUPPORT || secondSupportPerc < PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC)
            return null;

        if(!passDistinctFragmentsFilter(assemblySupport))
            return null;

        JunctionAssembly newAssembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualities(), assemblySupport,
                extensionSeqBuilder.repeatInfo());

        if(newAssembly.extensionLength() < ASSEMBLY_MIN_SOFT_CLIP_LENGTH)
            return null;

        // perform a final sequence comparison check with more liberal comparison tests
        boolean closeMatch = SequenceCompare.matchedAssemblySequences(firstAssembly, newAssembly);

        if(closeMatch)
            return null;

        if(LineUtils.hasLineSourceSequence(newAssembly))
            return null;

        newAssembly.setExtBaseBuildInfo(extensionSeqBuilder.buildInformation());

        addJunctionReads(newAssembly, extensionSeqBuilder, junctionReads);

        return newAssembly;
    }

    private void addJunctionReads(
            final JunctionAssembly assembly, final ExtensionSeqBuilderOld extensionSeqBuilder, final List<Read> junctionReads)
    {
        int mismatchReadCount = 0;

        // test other junction-spanning reads against this new assembly
        for(Read read : junctionReads)
        {
            if(assembly.support().stream().anyMatch(x -> x.cachedRead() == read)) // skip those already added
                continue;

            if(!canAddJunctionRead(assembly, extensionSeqBuilder, read))
                ++mismatchReadCount;
        }

        assembly.addMismatchReadCount(mismatchReadCount);
    }

    private boolean canAddJunctionRead(final JunctionAssembly assembly, final ExtensionSeqBuilderOld extensionSeqBuilder, final Read read)
    {
        ExtReadParseState readParseState = extensionSeqBuilder.checkAddJunctionRead(read);

        if(readParseState == null)
            return false;

        if(readParseState.exceedsMaxMismatches() || !extensionSeqBuilder.sufficientHighQualMatches(readParseState))
            return false;

        assembly.addSupport(read, JUNCTION, readParseState.junctionIndex(), readParseState.matchedBases(), readParseState.mismatches());
        return true;
    }

    private void addJunctionReads(
            final JunctionAssembly assembly, final ExtensionSeqBuilder extensionSeqBuilder, final List<Read> junctionReads)
    {
        int mismatchReadCount = 0;

        // test other junction-spanning reads against this new assembly
        for(Read read : junctionReads)
        {
            if(assembly.support().stream().anyMatch(x -> x.cachedRead() == read)) // skip those already added
                continue;

            ReadParseState readInfo = extensionSeqBuilder.checkAddJunctionRead(read);
            if(readInfo == null || readInfo.mismatched())
            {
                ++mismatchReadCount;
            }
            else
            {
                assembly.addSupport(read, JUNCTION, readInfo.startIndex(), readInfo.matchedBases(), readInfo.mismatchCount());
            }
        }

        assembly.addMismatchReadCount(mismatchReadCount);
    }

    private boolean meetsMinSupportThreshold(final List<SupportRead> support)
    {
        if(!aboveMinSupportThreshold(support))
            return false;

        return passDistinctFragmentsFilter(support);
    }

    private boolean aboveMinSupportThreshold(final List<SupportRead> support)
    {
        int minRequiredReadCount = minReadThreshold(mJunction);

        // account for overlapping fragments
        if(support.size() >= minRequiredReadCount * 2)
            return true;

        Set<String> uniqueReadIds = Sets.newHashSet();
        support.forEach(x -> uniqueReadIds.add(x.id()));

        return uniqueReadIds.size() >= minRequiredReadCount;
    }

    private boolean aboveMinReadThreshold(final List<Read> reads)
    {
        int minRequiredReadCount = minReadThreshold(mJunction);

        if(reads.size() >= minRequiredReadCount * 2)
            return true;

        Set<String> uniqueReadIds = Sets.newHashSet();
        reads.forEach(x -> uniqueReadIds.add(x.id()));

        return uniqueReadIds.size() >= minRequiredReadCount;
    }

    protected static int minReadThreshold(final Junction junction)
    {
        return junction.Hotspot ? MIN_HOTSPOT_JUNCTION_SUPPORT : ASSEMBLY_MIN_READ_SUPPORT;
    }

    private boolean passDistinctFragmentsFilter(final List<SupportRead> support)
    {
        if(mJunction.Hotspot)
            return true;

        Set<Integer> readPositions = Sets.newHashSet();
        Set<Integer> readEndPositions = null;

        for(SupportRead read : support)
        {
            if(read.isPairedRead())
            {
                readPositions.add(read.cachedRead().fivePrimeFragmentPosition());

                if(readPositions.size() >= ASSEMBLY_MIN_DISTINCT_FRAGS)
                    return true;
            }
            else
            {
                if(readEndPositions == null)
                    readEndPositions = Sets.newHashSet();

                readPositions.add(read.unclippedStart());
                readEndPositions.add(read.unclippedEnd());

                if(readPositions.size() >= ASSEMBLY_MIN_DISTINCT_FRAGS && readEndPositions.size() >= ASSEMBLY_MIN_DISTINCT_FRAGS)
                    return true;
            }
        }

        return false;
    }

    @VisibleForTesting
    public JunctionAssembler(final Junction junction)
    {
        this(junction, null);
    }
}
