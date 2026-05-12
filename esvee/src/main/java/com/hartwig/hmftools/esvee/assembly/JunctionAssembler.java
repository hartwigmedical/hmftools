package com.hartwig.hmftools.esvee.assembly;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_DISCORDANT_MIN_MAP_QUALITY;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH_LOWER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH_LOWER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_DISTINCT_FRAGS;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.MAX_OBSERVED_CONCORDANT_FRAG_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findIndelExtensionReads;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.hasIndelJunctionReads;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.isLineWithLocalAlignedInsert;
import static com.hartwig.hmftools.esvee.assembly.RemoteRegionFinder.addOrCreateMateRemoteRegion;
import static com.hartwig.hmftools.esvee.assembly.SeqTechUtils.findSbxPossibleDuplicates;
import static com.hartwig.hmftools.esvee.assembly.SeqTechUtils.passSbxDistinctPrimePositionsFilter;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.readJunctionExtensionLength;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.readSoftClipsAndCrossesJunction;
import static com.hartwig.hmftools.esvee.assembly.read.ReadUtils.recordSoftClipsAtJunction;
import static com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip.checkSupportVsRefSideSoftClip;
import static com.hartwig.hmftools.esvee.assembly.types.RemoteRegion.mergeRegions;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_EXTENSION_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.LINE_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.isSbx;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_HOTSPOT_JUNCTION_SUPPORT;

import java.util.ArrayList;
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
import com.hartwig.hmftools.esvee.common.SagaMatcher;

import org.jetbrains.annotations.Nullable;

public class JunctionAssembler
{
    private Junction mJunction;
    private final RefGenomeInterface mRefGenome;
    @Nullable
    private final SagaMatcher mSagaMatcher;
    private final List<Read> mNonJunctionReads;

    public JunctionAssembler(final Junction junction, final RefGenomeInterface refGenome, @Nullable final SagaMatcher sagaMatcher)
    {
        mJunction = junction;
        mRefGenome = refGenome;
        mSagaMatcher = sagaMatcher;
        mNonJunctionReads = Lists.newArrayList();
    }

    public List<Read> nonJunctionReads() { return mNonJunctionReads; }

    public List<JunctionAssembly> processJunction(final List<Read> rawReads)
    {
        // find prominent reads to establish the extension sequence, taking any read meeting min soft-clip lengths
        // and repetitive indels

        if(!mJunction.indelBased() && hasIndelJunctionReads(mJunction, rawReads))
        {
            // fall-back in case Prep didn't set this state or junctions are loaded from config
            mJunction.markAsIndel();
        }

        List<Read> junctionReads;
        List<Read> extensionReads;
        boolean hasMinLengthSoftClipRead;
        boolean usedRelaxedFilters;

        if(mJunction.indelBased())
        {
            junctionReads = new ArrayList<>();
            extensionReads = new ArrayList<>();
            findIndelExtensionReads(mJunction, rawReads, extensionReads, junctionReads, mNonJunctionReads);
            hasMinLengthSoftClipRead = !extensionReads.isEmpty();
            usedRelaxedFilters = false;
        }
        else if(mJunction.DiscordantOnly)
        {
            // look for a common soft-clip position, otherwise take the min variant length back from the inner most read as the junction
            junctionReads = new ArrayList<>();
            extensionReads = new ArrayList<>();
            assessDiscordantJunction(rawReads, extensionReads, junctionReads);
            hasMinLengthSoftClipRead = !extensionReads.isEmpty();
            usedRelaxedFilters = false;
        }
        else
        {
            // First, try with the regular soft clip length requirement.
            AssessJunctionReadsResult result = assessSoftClipJunction(rawReads, false);
            if(!checkJunctionReadExtension(result.hasMinLengthSoftClipRead, result.extensionReads))
            {
                // If that failed, try again with a relaxed requirement - later this will be permitted if the junction sequence matches SAGA.
                result = assessSoftClipJunction(rawReads, true);
            }
            junctionReads = result.junctionReads();
            extensionReads = result.extensionReads();
            mNonJunctionReads.addAll(result.nonJunctionReads());
            hasMinLengthSoftClipRead = result.hasMinLengthSoftClipRead();
            usedRelaxedFilters = result.usedRelaxedLimits();
        }

        if(!checkJunctionReadExtension(hasMinLengthSoftClipRead, extensionReads))
        {
            SV_LOGGER.trace("filter stage=junctionAssembly reason=\"read extension\" data=junction({})", mJunction);
            if(!mJunction.isSagaMatched())
            {
                return Collections.emptyList();
            }
        }

        if (usedRelaxedFilters)
        {
            SV_LOGGER.trace("junction({}) passed relaxed filters", mJunction);
        }

        List<Read> duplicateLongExtensionReads = isSbx() ? findSbxPossibleDuplicates(mJunction, extensionReads) : Collections.emptyList();
        duplicateLongExtensionReads.forEach(x -> extensionReads.remove(x));
        duplicateLongExtensionReads.forEach(x -> junctionReads.remove(x));

        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads, usedRelaxedFilters);

        int reqExtensionLength = extensionSeqBuilder.hasLineSequence() ? LINE_MIN_EXTENSION_LENGTH : (usedRelaxedFilters ? ASSEMBLY_MIN_SOFT_CLIP_LENGTH_LOWER :  ASSEMBLY_MIN_SOFT_CLIP_LENGTH);

        if(!extensionSeqBuilder.isValid() || extensionSeqBuilder.extensionLength() < reqExtensionLength)
        {
            SV_LOGGER.trace("filter stage=junctionAssembly reason=\"extension sequence\" data=junction({})", mJunction);
            if(!mJunction.isSagaMatched())
            {
                return Collections.emptyList();
            }
        }

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        JunctionAssembly firstAssembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualities(), assemblySupport,
                extensionSeqBuilder.repeats());

        if(!meetsMinSupportThreshold(assemblySupport))
        {
            SV_LOGGER.trace("filter stage=junctionAssembly reason=\"min support\" data=junction({})", mJunction);
            if(!mJunction.isSagaMatched())
            {
                return Collections.emptyList();
            }
        }

        // filter LINE source-type sites marked by opposition orientation poly A/T sequences
        if(!firstAssembly.indel() && LineUtils.hasLineSourceSequence(firstAssembly))
        {
            SV_LOGGER.trace("filter stage=junctionAssembly reason=\"LINE source site\" data=assembly({})", firstAssembly);
            if(!mJunction.isSagaMatched())
            {
                return Collections.emptyList();
            }
        }

        firstAssembly.setExtBaseBuildInfo(extensionSeqBuilder.buildInformation());

        if(extensionSeqBuilder.hasLineSequence())
            firstAssembly.markLineSequence();

        int initialAssemblySupport = assemblySupport.size();
        addJunctionReads(firstAssembly, extensionSeqBuilder, junctionReads);

        if(firstAssembly.hasLineSequence())
        {
            if(isLineWithLocalAlignedInsert(firstAssembly))
            {
                firstAssembly.unmarkLineSequence();

                if(firstAssembly.extensionLength() < ASSEMBLY_MIN_SOFT_CLIP_LENGTH)
                {
                    SV_LOGGER.trace("filter stage=junctionAssembly reason=\"LINE with local aligned insert SC\" data=assembly({})", firstAssembly);
                    if(!mJunction.isSagaMatched())
                    {
                        return Collections.emptyList();
                    }
                }
            }
        }

        // test for a second well-supported, alternative assembly at the same junction
        JunctionAssembly secondAssembly = checkSecondAssembly(extensionSeqBuilder.mismatchReads(), firstAssembly, junctionReads);

        List<JunctionAssembly> assemblies = Lists.newArrayList(firstAssembly);
        if(secondAssembly != null)
        {
            assemblies.add(secondAssembly);

            if(!keepSecondAssembly(secondAssembly.supportCount(), initialAssemblySupport))
                assemblies.remove(firstAssembly);
        }

        if (firstAssembly.bases().length == 0)
        {
            assemblies.remove(firstAssembly);
        }

        for(JunctionAssembly assembly : assemblies)
        {
            // call routine to purge reads likely add to multiple assemblies and breaching a dominant RSSC
            // NOTE: this could be done for all assemblies, not just split ones
            if(assemblies.size() > 1)
                checkSupportVsRefSideSoftClip(assembly);

            RefBaseSeqBuilder refBaseSeqBuilder = new RefBaseSeqBuilder(assembly);
            assembly.setRefBases(refBaseSeqBuilder);

            assembly.buildRepeatInfo();

            // Now we have the assembly sequence, try to match to SAGA if necessary.
            // Note matching to SAGA via coordinate is not enough; we need to know that the variant sequence is actually the same.
            boolean sagaMatched = tryMatchToSaga(assembly);
            if (usedRelaxedFilters)
            {
                if (sagaMatched)
                {
                    SV_LOGGER.trace("assembly({}) recovered with SAGA", assembly);
                    assembly.mSagaRecovered = true;
                }
                else
                {
                    // Failed the regular filters and couldn't recover by matching to SAGA.
                    SV_LOGGER.trace("filter stage=junctionAssembly reason=\"could not recover with SAGA\" data=assembly({})", assembly);
                    return Collections.emptyList();
                }
            }
        }

        return assemblies;
    }

    private record AssessJunctionReadsResult(
            List<Read> junctionReads,
            List<Read> extensionReads,
            List<Read> nonJunctionReads,
            boolean hasMinLengthSoftClipRead,
            boolean usedRelaxedLimits
    )
    {}

    private AssessJunctionReadsResult assessSoftClipJunction(final List<Read> rawReads, boolean useRelaxedLimits)
    {
        // the only difference for indel-based junctions is that only the long indels are used to build the consensus extension

        int minSoftClip = useRelaxedLimits ? ASSEMBLY_MIN_SOFT_CLIP_LENGTH_LOWER : ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
        int minSoftClipSecondary = useRelaxedLimits ? ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH_LOWER : ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;

        List<Read> extensionReads = new ArrayList<>();
        List<Read> junctionReads = new ArrayList<>();
        List<Read> nonJunctionReads = new ArrayList<>();
        boolean hasMinLengthSoftClipRead = false;

        for(Read read : rawReads)
        {
            if(!readSoftClipsAndCrossesJunction(read, mJunction, mRefGenome))
            {
                nonJunctionReads.add(read);
                continue;
            }

            if(recordSoftClipsAtJunction(read, mJunction))
            {
                int softClipJunctionExtension = readJunctionExtensionLength(read, mJunction);

                if(read.hasLineTail(mJunction.isForward()))
                {
                    hasMinLengthSoftClipRead |= softClipJunctionExtension >= LINE_MIN_EXTENSION_LENGTH;

                    if(softClipJunctionExtension >= LINE_MIN_SOFT_CLIP_SECONDARY_LENGTH)
                        extensionReads.add(read);
                }
                else
                {
                    hasMinLengthSoftClipRead |= softClipJunctionExtension >= minSoftClip;

                    if(softClipJunctionExtension >= minSoftClipSecondary)
                        extensionReads.add(read);
                }
            }

            junctionReads.add(read);
        }

        return new AssessJunctionReadsResult(junctionReads, extensionReads, nonJunctionReads, hasMinLengthSoftClipRead, useRelaxedLimits);
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

        SV_LOGGER.trace("junction({}) discordant position adjust: {}", mJunction, adjustedJuncPosition);

        mJunction = new Junction(
                mJunction.Chromosome, adjustedJuncPosition, mJunction.Orient, true, false, false, mJunction.sagaMatch());

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

    private boolean checkJunctionReadExtension(boolean hasMinLengthSoftClipRead, final List<Read> extensionReads)
    {
        return hasMinLengthSoftClipRead && aboveMinReadThreshold(extensionReads);
    }

    private boolean tryMatchToSaga(JunctionAssembly assembly)
    {
        if(mSagaMatcher == null)
        {
            return false;
        }
        else
        {
            int junctionOffset = mJunction.isForward() ? assembly.refBaseLength() : assembly.baseLength() - assembly.refBaseLength();
            SagaMatcher.MatchBySequence sagaSeqMatch = mSagaMatcher.matchBySequence(assembly.bases(), List.of(junctionOffset));
            SV_LOGGER.trace("assembly({}) SAGA sequence match {}", assembly, sagaSeqMatch);
            assembly.setSagaMatch(sagaSeqMatch);
            return sagaSeqMatch != null;
        }
    }

    private JunctionAssembly checkSecondAssembly(
            final List<Read> extensionReads, final JunctionAssembly firstAssembly, final List<Read> junctionReads)
    {
        if(extensionReads.isEmpty() || mJunction.DiscordantOnly)
            return null;

        if(firstAssembly.hasLineSequence())
            return null;

        int secondSupport = extensionReads.size();
        double secondSupportPerc = secondSupport / (double)firstAssembly.supportCount();

        if(secondSupport < ASSEMBLY_SPLIT_MIN_READ_SUPPORT || secondSupportPerc < PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC)
            return null;

        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads);

        if(!extensionSeqBuilder.isValid())
            return null;

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        // test min support again from actual supporting reads
        if(!keepSecondAssembly(firstAssembly.supportCount(), assemblySupport.size()))
            return null;

        if(!passDistinctFragmentsFilter(assemblySupport))
            return null;

        JunctionAssembly newAssembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualities(), assemblySupport,
                extensionSeqBuilder.repeats());

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

    public boolean keepSecondAssembly(final int firstSupportCount, final int secondSupportCount)
    {
        double secondSupportPerc = secondSupportCount / (double)firstSupportCount;

        return secondSupportCount >= ASSEMBLY_SPLIT_MIN_READ_SUPPORT && secondSupportPerc >= PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;
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
            else if(extensionSeqBuilder.sufficientQualMatches(readInfo))
            {
                assembly.addJunctionSupport(read, JUNCTION, readInfo.startIndex(), readInfo);
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

        if(isSbx() && !passSbxDistinctPrimePositionsFilter(support))
            return false;

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
        this(junction, null, null);
    }
}
