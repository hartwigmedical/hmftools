package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.buildIndelFrequencies;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findIndelExtensionReads;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findMaxFrequencyIndelReads;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.hasIndelJunctionReads;
import static com.hartwig.hmftools.esvee.assembly.RemoteRegionFinder.addOrCreateMateRemoteRegion;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.readJunctionExtensionLength;
import static com.hartwig.hmftools.esvee.assembly.read.ReadFilters.recordSoftClipsAtJunction;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteRegion.mergeRegions;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_EXTENSION_LENGTH;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.assembly.types.RemoteRegion;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.read.ReadFilters;

public class JunctionAssembler
{
    private Junction mJunction;
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
            Map<Integer,List<Read>> indelLengthReads = Maps.newHashMap();

            // the only difference for indel-based junctions is that only the long indels are used to build the consensus extension

            for(Read read : rawReads)
            {
                if(!ReadFilters.recordSoftClipsAndCrossesJunction(read, mJunction))
                {
                    mNonJunctionReads.add(read);
                    continue;
                }

                if(recordSoftClipsAtJunction(read, mJunction))
                {
                    int softClipJunctionExtension = readJunctionExtensionLength(read, mJunction);

                    hasMinLengthSoftClipRead |= softClipJunctionExtension >= ASSEMBLY_MIN_SOFT_CLIP_LENGTH
                            || (read.hasLineTail() && softClipJunctionExtension >= LINE_MIN_EXTENSION_LENGTH);

                    if(softClipJunctionExtension >= ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH)
                    {
                        extensionReads.add(read);
                    }
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

        if(!hasMinLengthSoftClipRead || !aboveMinReadThreshold(extensionReads))
            return Collections.emptyList();

        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads);

        int reqExtensionLength = extensionSeqBuilder.hasLineSequence() ? LINE_MIN_EXTENSION_LENGTH : ASSEMBLY_MIN_SOFT_CLIP_LENGTH;

        if(!extensionSeqBuilder.isValid() || extensionSeqBuilder.extensionLength() < reqExtensionLength)
            return Collections.emptyList();

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        if(!aboveMinSupportThreshold(assemblySupport))
            return Collections.emptyList();

        JunctionAssembly firstAssembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualities(), assemblySupport,
                extensionSeqBuilder.repeatInfo());

        if(extensionSeqBuilder.hasLineSequence())
            firstAssembly.markLineSequence();

        List<JunctionAssembly> assemblies = Lists.newArrayList(firstAssembly);

        addJunctionReads(firstAssembly, extensionSeqBuilder, junctionReads);

        // test for a second well-supported, alternative assembly at the same junction
        JunctionAssembly secondAssembly = checkSecondAssembly(extensionSeqBuilder.mismatchReads(), firstAssembly, junctionReads);

        if(secondAssembly != null)
            assemblies.add(secondAssembly);

        for(JunctionAssembly assembly : assemblies)
        {
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

        for(Read read : rawReads)
        {
            if(read.mappingQuality() < ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY)
                continue;

            // check that the read maps to the same remote region to qualify as an junction candidate
            if(!mainRemoteRegion.overlaps(read.mateChromosome(), read.mateAlignmentStart(), read.mateAlignmentEnd()))
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

        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads);

        if(!extensionSeqBuilder.isValid())
            return null;

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        int secondSupport = assemblySupport.size();
        double secondSupportPerc = secondSupport / (double)firstAssembly.supportCount();

        if(secondSupport < ASSEMBLY_SPLIT_MIN_READ_SUPPORT || secondSupportPerc < PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC)
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

        addJunctionReads(newAssembly, extensionSeqBuilder, junctionReads);

        return newAssembly;
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

            if(!canAddJunctionRead(assembly, extensionSeqBuilder, read))
                ++mismatchReadCount;
        }

        assembly.addMismatchReadCount(mismatchReadCount);
    }

    private boolean canAddJunctionRead(final JunctionAssembly assembly, final ExtensionSeqBuilder extensionSeqBuilder, final Read read)
    {
        ExtReadParseState readParseState = extensionSeqBuilder.checkAddJunctionRead(read);

        if(readParseState == null)
            return false;

        if(readParseState.exceedsMaxMismatches() || readParseState.highQualMatches() < ASSEMBLY_MIN_EXTENSION_READ_HIGH_QUAL_MATCH)
            return false;

        assembly.addSupport(read, JUNCTION, readParseState.junctionIndex(), readParseState.matchedBases(), readParseState.mismatches());
        return true;
    }

    private static boolean aboveMinSupportThreshold(final List<SupportRead> assemblySupport)
    {
        // account for overlapping fragments
        if(assemblySupport.size() >= ASSEMBLY_MIN_READ_SUPPORT * 2)
            return true;

        Set<String> uniqueReadIds = Sets.newHashSet();
        assemblySupport.forEach(x -> uniqueReadIds.add(x.id()));

        return uniqueReadIds.size() >= ASSEMBLY_MIN_READ_SUPPORT;
    }

    private static boolean aboveMinReadThreshold(final List<Read> reads)
    {
        if(reads.size() >= ASSEMBLY_MIN_READ_SUPPORT * 2)
            return true;

        Set<String> uniqueReadIds = Sets.newHashSet();
        reads.forEach(x -> uniqueReadIds.add(x.id()));

        return uniqueReadIds.size() >= ASSEMBLY_MIN_READ_SUPPORT;
    }
}
