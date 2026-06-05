package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MAX_JUNCTION_READ_LOW_QUAL_MISMATCH_PERC;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_DISTINCT_FRAGS;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.findIndelExtensionReads;
import static com.hartwig.hmftools.esvee.assembly.IndelBuilder.hasIndelJunctionReads;
import static com.hartwig.hmftools.esvee.assembly.LineUtils.isLineWithLocalAlignedInsert;
import static com.hartwig.hmftools.esvee.assembly.SeqTechUtils.findSbxPossibleDuplicates;
import static com.hartwig.hmftools.esvee.assembly.SeqTechUtils.passSbxDistinctPrimePositionsFilter;
import static com.hartwig.hmftools.esvee.assembly.types.RefSideSoftClip.checkSupportVsRefSideSoftClip;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.EXTENSION;
import static com.hartwig.hmftools.esvee.assembly.types.SupportType.JUNCTION;
import static com.hartwig.hmftools.esvee.common.SvConstants.isSbx;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.MIN_HOTSPOT_JUNCTION_SUPPORT;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.esvee.assembly.read.Read;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.common.saga.SagaMatchBySequence;
import com.hartwig.hmftools.esvee.common.saga.SagaSequenceMatcher;

import org.jetbrains.annotations.Nullable;

public class JunctionAssembler
{
    private Junction mJunction;
    private final SagaSequenceMatcher mSagaMatcher;
    private final List<Read> mNonJunctionReads;

    public JunctionAssembler(final Junction junction, @Nullable final SagaSequenceMatcher sagaMatcher)
    {
        mJunction = junction;
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

        List<Read> junctionReads = Lists.newArrayList();
        List<Read> extensionReads = Lists.newArrayList();

        JunctionThresholdState juncThresholdState = new JunctionThresholdState();

        // set defaults which may be overridden for each type of junction
        juncThresholdState.MinRequiredReads = minReadThreshold(mJunction);
        juncThresholdState.MinExtensionLength = ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
        juncThresholdState.MinSecondExtensionLength = ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;

        boolean isValidAssembly = false;

        if(mJunction.indelBased())
        {
            findIndelExtensionReads(mJunction, rawReads, extensionReads, junctionReads, mNonJunctionReads);

            if(extensionReads.isEmpty())
                return Collections.emptyList();

            juncThresholdState.ExtensionLengthValid = true;
            juncThresholdState.SecondExtensionLengthValid = true;
            // juncValidity.ExtensionLengthValid = !extensionReads.isEmpty();
        }
        else if(mJunction.DiscordantOnly)
        {
            // look for a common soft-clip position, otherwise take the min variant length back from the inner most read as the junction
            DiscordantJunctionBuilder discordantJunctionBuilder = new DiscordantJunctionBuilder(mJunction);
            discordantJunctionBuilder.assessDiscordantJunction(rawReads, extensionReads, junctionReads);

            if(extensionReads.isEmpty())
                return Collections.emptyList();

            if(discordantJunctionBuilder.hasNewJunction())
                mJunction = discordantJunctionBuilder.newJunction();

            juncThresholdState.ExtensionLengthValid = true;
            juncThresholdState.SecondExtensionLengthValid = true;
            // juncValidity.ExtensionLengthValid = !extensionReads.isEmpty();
        }
        else
        {
            JunctionReadTypes juncReadTypes = new JunctionReadTypes(mJunction, rawReads);
            junctionReads.addAll(juncReadTypes.junctionReads());

            juncReadTypes.establishExtensionReads(extensionReads, juncThresholdState);

            if(!juncThresholdState.ExtensionLengthValid)
                return Collections.emptyList();

            if(!juncThresholdState.SecondExtensionLengthValid)
            {
                if(juncThresholdState.IsLINE)
                    return Collections.emptyList();

                juncThresholdState.RequireNonExtensionSupport = true;
            }

            mNonJunctionReads.addAll(juncReadTypes.nonJunctionReads());

            /*
            // First, try with the regular soft clip length requirement.
            AssessJunctionReadsResult result = assessSoftClipJunction(rawReads, false);
            if(!meetsMinExtensionConditions(result.hasMinLengthSoftClipRead, result.extensionReads))
            {
                // If that failed, try again with a relaxed requirement - later this will be permitted if the junction sequence matches SAGA.
                result = assessSoftClipJunction(rawReads, true);
            }

            junctionReads = result.junctionReads();
            extensionReads = result.extensionReads();
            mNonJunctionReads.addAll(result.nonJunctionReads());
            hasMinLengthSoftClipRead = result.hasMinLengthSoftClipRead();
            usedRelaxedFilters = result.usedRelaxedLimits();

            // int minSoftClip = useRelaxedLimits ? ASSEMBLY_MIN_SOFT_CLIP_LENGTH_LOWER : ASSEMBLY_MIN_SOFT_CLIP_LENGTH;

            // int minSoftClipSecondary = useRelaxedLimits
            //        ? ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH_LOWER : ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;

            // hasMinLengthRead = juncReadTypes.hasMinLengthExtensionRead(ASSEMBLY_MIN_SOFT_CLIP_LENGTH);
            // usedRelaxedFilters = false;
            */
        }

        juncThresholdState.MinReadsValid = extensionReads.size() >= juncThresholdState.MinRequiredReads;

        if(!juncThresholdState.MinReadsValid && !juncThresholdState.RequireNonExtensionSupport)
        {
            return Collections.emptyList();
        }

        /* unique reads by ID is now checked later on
        if(!meetsMinExtensionConditions(hasMinLengthRead, extensionReads))
        {
            return Collections.emptyList();
        }
        */

        List<Read> duplicateLongExtensionReads = isSbx() ? findSbxPossibleDuplicates(mJunction, extensionReads) : Collections.emptyList();
        duplicateLongExtensionReads.forEach(x -> extensionReads.remove(x));
        duplicateLongExtensionReads.forEach(x -> junctionReads.remove(x));

        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads);

        extensionSeqBuilder.checkValidity(juncThresholdState);

        if(!juncThresholdState.ExtensionLengthValid)
            return Collections.emptyList();

        if((!juncThresholdState.SecondExtensionLengthValid || !juncThresholdState.MinReadsValid)
        && !juncThresholdState.RequireNonExtensionSupport)
        {
            return Collections.emptyList();
        }

        List<SupportRead> assemblySupport = extensionSeqBuilder.formAssemblySupport();

        /*
        int reqExtensionLength = extensionSeqBuilder.hasLineSequence()
                ? LINE_MIN_EXTENSION_LENGTH
                : (usedRelaxedFilters ? ASSEMBLY_MIN_SOFT_CLIP_LENGTH_LOWER : ASSEMBLY_MIN_SOFT_CLIP_LENGTH);

        if(!extensionSeqBuilder.isValid() || extensionSeqBuilder.extensionLength() < reqExtensionLength)
            return Collections.emptyList();


        if(!meetsMinSupportThreshold(assemblySupport))
            return Collections.emptyList();
        */

        JunctionAssembly firstAssembly = new JunctionAssembly(
                mJunction, extensionSeqBuilder.extensionBases(), extensionSeqBuilder.baseQualities(), assemblySupport,
                extensionSeqBuilder.repeats());

        // filter LINE source-type sites marked by opposition orientation poly A/T sequences
        if(!firstAssembly.indel() && LineUtils.hasLineSourceSequence(firstAssembly))
            return Collections.emptyList();

        firstAssembly.setExtBaseBuildInfo(extensionSeqBuilder.buildInformation());

        if(extensionSeqBuilder.hasLineSequence())
            firstAssembly.markLineSequence();

        int initialAssemblySupport = assemblySupport.size();
        addJunctionReads(firstAssembly, extensionSeqBuilder, junctionReads);

        if(juncThresholdState.RequireNonExtensionSupport)
        {
            // evaluate if any high-match junction reads can provide min support
            for(SupportRead read : firstAssembly.support())
            {
                if(assemblySupport.contains(read))
                    continue;

                if(read.extensionBaseMismatches() > 0)
                    continue;

                int baseMatches = read.extensionBaseMatches();
                int extensionLength = read.extensionLength(mJunction.Orient);
                double matchPerc = baseMatches / (double)extensionLength;

                if(matchPerc < ASSEMBLY_MAX_JUNCTION_READ_LOW_QUAL_MISMATCH_PERC)
                    continue;

                if(extensionLength >= juncThresholdState.MinSecondExtensionLength)
                {
                    juncThresholdState.SecondExtensionLengthValid = true;
                    assemblySupport.add(read);
                }
            }

            juncThresholdState.MinReadsValid = assemblySupport.size() >= juncThresholdState.MinRequiredReads;
        }

        if(!meetsMinSupportThreshold(assemblySupport, juncThresholdState.MinRequiredReads))
            return Collections.emptyList();

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

        List<JunctionAssembly> assemblies = Lists.newArrayList(firstAssembly);
        if(secondAssembly != null)
        {
            assemblies.add(secondAssembly);

            if(!keepSecondAssembly(secondAssembly.supportCount(), initialAssemblySupport))
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

            // check for a SAGA match by sequence, even if a coord match was previously found
            boolean sagaMatched = matchJunctionAssemblyToSaga(assembly);

            if(juncThresholdState.UsesLowerSagaLimits)
            {
                if(sagaMatched)
                {
                    SV_LOGGER.trace("assembly({}) recovered with SAGA", assembly);
                }
                else
                {
                    return Collections.emptyList();
                }
            }
        }

        return assemblies;
    }

    private JunctionAssembly checkSecondAssembly(
            final List<Read> extensionReads, final JunctionAssembly firstAssembly, final List<Read> junctionReads)
    {
        if(extensionReads.isEmpty() || mJunction.DiscordantOnly)
            return null;

        if(firstAssembly.hasLineSequence())
            return null;

        int secondSupport = extensionReads.size();
        double secondSupportPerc = secondSupport / (double) firstAssembly.supportCount();

        JunctionThresholdState juncThresholdState = new JunctionThresholdState();
        juncThresholdState.MinRequiredReads = ASSEMBLY_SPLIT_MIN_READ_SUPPORT;
        juncThresholdState.MinExtensionLength = ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
        juncThresholdState.MinSecondExtensionLength = ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;

        if(secondSupport < juncThresholdState.MinRequiredReads || secondSupportPerc < PRIMARY_ASSEMBLY_SPLIT_MIN_READ_SUPPORT_PERC)
            return null;

        ExtensionSeqBuilder extensionSeqBuilder = new ExtensionSeqBuilder(mJunction, extensionReads);
        extensionSeqBuilder.checkValidity(juncThresholdState);

        if(!juncThresholdState.isValid())
            return null;

        // if(!extensionSeqBuilder.isValid())
        //    return null;

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
        double secondSupportPerc = secondSupportCount / (double) firstSupportCount;

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

    private boolean meetsMinSupportThreshold(final List<SupportRead> support, int minRequiredReadCount)
    {
        if(!aboveMinSupportThreshold(support, minRequiredReadCount))
            return false;

        return passDistinctFragmentsFilter(support);
    }

    private boolean aboveMinSupportThreshold(final List<SupportRead> support, int minRequiredReadCount)
    {
        // account for overlapping fragments
        if(support.size() >= minRequiredReadCount * 2)
            return true;

        Set<String> uniqueReadIds = Sets.newHashSet();
        support.forEach(x -> uniqueReadIds.add(x.id()));

        return uniqueReadIds.size() >= minRequiredReadCount;
    }

    @Deprecated
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

    private boolean matchJunctionAssemblyToSaga(final JunctionAssembly assembly)
    {
        if(mSagaMatcher == null)
            return false;

        int junctionOffset = mJunction.isForward() ? assembly.refBaseLength() : assembly.baseLength() - assembly.refBaseLength();
        SagaMatchBySequence match = mSagaMatcher.matchBySequence(assembly.bases(), List.of(junctionOffset));
        assembly.setSagaMatch(match);
        return match != null;
    }

    @VisibleForTesting
    public JunctionAssembler(final Junction junction)
    {
        this(junction, null);
    }

    /*
    private record AssessJunctionReadsResult(
            List<Read> junctionReads,
            List<Read> extensionReads,
            List<Read> nonJunctionReads,
            boolean hasMinLengthSoftClipRead,
            boolean usedRelaxedLimits
    )
    {
    }

    private AssessJunctionReadsResult assessSoftClipJunction(final List<Read> rawReads, boolean useRelaxedLimits)
    {
        // the only difference for indel-based junctions is that only the long indels are used to build the consensus extension
        int minSoftClip = useRelaxedLimits ? ASSEMBLY_MIN_SOFT_CLIP_LENGTH_LOWER : ASSEMBLY_MIN_SOFT_CLIP_LENGTH;

        int minSoftClipSecondary = useRelaxedLimits
                ? ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH_LOWER : ASSEMBLY_MIN_SOFT_CLIP_SECONDARY_LENGTH;

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
    */

}
