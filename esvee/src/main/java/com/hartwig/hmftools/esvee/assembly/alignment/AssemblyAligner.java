package com.hartwig.hmftools.esvee.assembly.alignment;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.SvUtils.isShortLocalDelDupIns;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ALIGNMENT_REQUERY_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_MIN_SOFT_CLIP_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.ASSEMBLY_UNLINKED_WEAK_ASSEMBLY_EXTENSION_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.WEAK_ASSEMBLY_DISC_RATE;
import static com.hartwig.hmftools.esvee.common.SvConstants.WEAK_ASSEMBLY_LONG_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.WEAK_ASSEMBLY_MIN_INS_LENGTH;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.WEAK_ASSEMBLY_MIN_REF_REPEAT;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.WEAK_ASSEMBLY_UNPAIRED_LONG_EXT_FACTOR;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.WEAK_ASSEMBLY_UNPAIRED_LONG_EXT_FACTOR_LOWER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.WEAK_ASSEMBLY_UNPAIRED_MAX_READS;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.WEAK_ASSEMBLY_UNPAIRED_READ_FACTOR;
import static com.hartwig.hmftools.esvee.assembly.SeqTechUtils.SBX_MEDIUM_QUAL_DESYNC_COUNT;
import static com.hartwig.hmftools.esvee.assembly.SeqTechUtils.SBX_PRIME_POSITION_RANGE_THRESHOLD;
import static com.hartwig.hmftools.esvee.assembly.alignment.Alignment.writeAssemblyData;
import static com.hartwig.hmftools.esvee.assembly.alignment.AssemblyAlignment.isLocalIndelAssembly;
import static com.hartwig.hmftools.esvee.common.SvConstants.INV_SHORT_MAX_HOMOLOGY_HIGHER;
import static com.hartwig.hmftools.esvee.common.SvConstants.isIllumina;

import java.util.Collections;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bwa.IBwaMemAligner;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.perf.TaskQueue;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.assembly.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.output.AlignmentWriter;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.RepeatInfo;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.assembly.types.ThreadTask;
import com.hartwig.hmftools.esvee.common.saga.SagaAlignment;
import com.hartwig.hmftools.esvee.common.saga.SagaJunctionInfo;
import com.hartwig.hmftools.esvee.common.saga.SagaSequenceMatch;
import com.hartwig.hmftools.esvee.common.saga.SagaSequenceMatcher;
import com.hartwig.hmftools.esvee.common.saga.SagaMatcherFactory;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.Nullable;

public class AssemblyAligner extends ThreadTask
{
    private final AssemblyConfig mConfig;
    private final IBwaMemAligner mAligner;
    private final AlignmentWriter mWriter;

    private final TaskQueue<AssemblyAlignment> mAssemblyAlignments;
    private int mRequeriedSuppCount;
    private int mRequeriedSoftClipCount;

    @Nullable
    private final SagaSequenceMatcher mSagaMatcher;

    public AssemblyAligner(
            final AssemblyConfig config, final IBwaMemAligner aligner, @Nullable final SagaMatcherFactory sagaMatcherFactory,
            final AlignmentWriter writer, final TaskQueue<AssemblyAlignment> assemblyAlignments)
    {
        super("AssemblerAlignment");
        mConfig = config;
        mAligner = aligner;
        mWriter = writer;
        mAssemblyAlignments = assemblyAlignments;
        mRequeriedSuppCount = 0;
        mRequeriedSoftClipCount = 0;
        mSagaMatcher = sagaMatcherFactory == null ? null : sagaMatcherFactory.createSequenceMatcher();
    }

    public int requeriedSuppCount()
    {
        return mRequeriedSuppCount;
    }
    public int requeriedSoftClipCount()
    {
        return mRequeriedSoftClipCount;
    }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                mPerfCounter.start();

                AssemblyAlignment assemblyAlignment = mAssemblyAlignments.removeItem();

                processAssembly(assemblyAlignment);

                stopCheckLog(
                        format("alignment count(%d) assemblies(%s)", assemblyAlignment.assemblies().size(), assemblyAlignment.info()),
                        mConfig.PerfLogTime);
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all alignment tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    @VisibleForTesting
    public void processAssembly(final AssemblyAlignment assemblyAlignment)
    {
        if(!assemblyAlignment.isValid())
        {
            SV_LOGGER.warn("assembly alignment({}) invalid, skipping", assemblyAlignment);
            return;
        }

        if(assemblyAlignment.isMerged())
        {
            writeAssemblyData(mWriter, mConfig, assemblyAlignment, Collections.emptyList(), Collections.emptyList());
            return;
        }

        AlignmentInfo alignmentInfo = null;
        // Try to form breakends from the matched SAGA variant first. Otherwise, create breakends from alignment to ref genome.
        if(tryMatchAssemblyAlignmentToSaga(assemblyAlignment))
        {
            alignmentInfo = processSagaMatchedAssembly(assemblyAlignment);
        }
        if(alignmentInfo == null)
        {
            alignmentInfo = alignAssembly(assemblyAlignment);
        }

        AlignmentFragments alignmentFragments = new AlignmentFragments(assemblyAlignment, mConfig.combinedSampleIds());
        alignmentFragments.allocateBreakendSupport();

        writeAssemblyData(mWriter, mConfig, assemblyAlignment, alignmentInfo.alignments(), alignmentInfo.requeriedAlignments());
    }

    private boolean tryMatchAssemblyAlignmentToSaga(AssemblyAlignment assemblyAlignment)
    {
        if(mSagaMatcher == null)
        {
            return false;
        }

        // Only use SAGA for the simple case of 0 or 1 links. Chained assemblies require more care, not implemented for now.
        boolean isChained = assemblyAlignment.phaseSet() != null && assemblyAlignment.phaseSet().assemblyLinks().size() > 1;
        if(isChained)
        {
            return false;
        }

        List<SagaJunctionInfo> junctionInfos = assemblyAlignment.linkIndices().stream().map(SagaJunctionInfo::new).toList();
        // If an assembly is a LINE site, then the extension sequence requirement is lower.
        // In which case need to lower the requirement for the junction overlap, or else there can be false negative matches.
        boolean lowerJunctionOverlap = assemblyAlignment.assemblies().stream().anyMatch(JunctionAssembly::hasLineSequence);
        SagaSequenceMatch sagaMatch = mSagaMatcher.matchBySequence(assemblyAlignment.fullSequence()
                .getBytes(), junctionInfos, lowerJunctionOverlap, false);
        assemblyAlignment.setSagaMatch(sagaMatch);

        return sagaMatch != null;
    }

    private record AlignmentInfo(
            List<AlignData> alignments,
            List<AlignData> requeriedAlignments
    )
    {
    }

    @Nullable
    private AlignmentInfo processSagaMatchedAssembly(final AssemblyAlignment assemblyAlignment)
    {
        SagaSequenceMatch sagaMatch = assemblyAlignment.sagaMatch();
        if(sagaMatch == null)
        {
            return null;
        }

        // FIXME? validate that the saga match is nearby to the assembly.

        SagaAlignment sagaAlignment = sagaMatch.alignment();
        AlignData alignData = AlignData.fromSaga(sagaAlignment);
        alignData.setFullSequenceData(assemblyAlignment.fullSequence(), assemblyAlignment.fullSequenceLength());

        BreakendBuilder breakendBuilder = new BreakendBuilder(mConfig.RefGenome, assemblyAlignment);
        boolean success = breakendBuilder.formBreakendsFromSaga(sagaAlignment, alignData);
        if(!success)
        {
            return null;
        }

        return new AlignmentInfo(List.of(alignData), Collections.emptyList());
    }

    private AlignmentInfo alignAssembly(final AssemblyAlignment assemblyAlignment)
    {
        List<BwaMemAlignment> bwaAlignments = mAligner.alignSequence(assemblyAlignment.fullSequence().getBytes());

        List<AlignData> alignments = bwaAlignments.stream()
                .map(x -> AlignData.fromRef(x, mConfig.RefGenVersion))
                .filter(x -> x != null).collect(Collectors.toList());

        // set the orientation-adjusted sequence coordinates - done here since used in requery logic
        String fullSequence = assemblyAlignment.fullSequence();
        alignments.forEach(x -> x.setFullSequenceData(fullSequence, assemblyAlignment.fullSequenceLength()));

        List<AlignData> requeriedAlignments = Lists.newArrayList();

        alignments = requerySoftClipAlignments(assemblyAlignment, alignments);

        alignments = requerySupplementaryAlignments(assemblyAlignment, alignments, requeriedAlignments);

        if(!alignments.isEmpty())
        {
            BreakendBuilder breakendBuilder = new BreakendBuilder(mConfig.RefGenome, assemblyAlignment);
            breakendBuilder.formBreakends(alignments);

            // final filters on assembly and alignment results
            if(isWeakSingleReadExtensionAssembly(assemblyAlignment))
            {
                assemblyAlignment.breakends().clear();
            }
        }

        return new AlignmentInfo(alignments, requeriedAlignments);
    }

    private List<AlignData> requerySoftClipAlignments(final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments)
    {
        // re-align long soft-clipped sequences and attach them to the original alignment
        String softClipBases;
        boolean isLeftClip;
        boolean firstBasesMissing;

        AlignData relevantAlignment;
        int fullLength = assemblyAlignment.fullSequenceLength();
        int minSeqStart = alignments.stream().map(AlignData::sequenceStart).min(Integer::compare).orElse(0);
        int maxSeqEnd = alignments.stream().map(AlignData::sequenceEnd).max(Integer::compare).orElse(fullLength);
        int missingBasesOnEnd = fullLength - maxSeqEnd - 1;

        if(minSeqStart >= ALIGNMENT_REQUERY_SOFT_CLIP_LENGTH && minSeqStart >= missingBasesOnEnd)
        {
            relevantAlignment = alignments.stream().filter(x -> x.sequenceStart() == minSeqStart).findFirst().orElse(null);
            isLeftClip = relevantAlignment.orientation().isForward();
            firstBasesMissing = true;

        }
        else if(missingBasesOnEnd >= ALIGNMENT_REQUERY_SOFT_CLIP_LENGTH)
        {
            relevantAlignment = alignments.stream().filter(x -> x.sequenceEnd() == maxSeqEnd).findFirst().orElse(null);
            isLeftClip = relevantAlignment.orientation().isReverse();
            firstBasesMissing = false;
        }
        else
        {
            return alignments;
        }

        int softClipLength = isLeftClip ? relevantAlignment.leftSoftClipLength() : relevantAlignment.rightSoftClipLength();
        if(firstBasesMissing)
        {
            softClipBases = assemblyAlignment.fullSequence().substring(0, softClipLength);
        }
        else
        {
            softClipBases = assemblyAlignment.fullSequence().substring(fullLength - softClipLength);
        }

        ++mRequeriedSoftClipCount;

        List<BwaMemAlignment> requeryBwaAlignments = mAligner.alignSequence(softClipBases.getBytes());

        List<AlignData> newAlignments = requeryBwaAlignments.stream()
                .map(x -> AlignData.fromRef(x, mConfig.RefGenVersion))
                .filter(x -> x != null).collect(Collectors.toList());

        if(newAlignments.isEmpty())
            return alignments;

        SV_LOGGER.trace("assembly({}) requeried single alignment({}) with long soft-clip", assemblyAlignment, relevantAlignment);

        int softClipSeqIndexStart;

        if(isLeftClip == relevantAlignment.orientation().isForward())
        {
            softClipSeqIndexStart = max(relevantAlignment.sequenceStart() - softClipLength, 0);
        }
        else
        {
            softClipSeqIndexStart = relevantAlignment.sequenceEnd() + 1;
        }

        for(AlignData rqAlignment : newAlignments)
        {
            // adjust values to be in terms of the original sequence, and then ordered within the soft-clip bounds
            int adjSequenceStart, adjSequenceEnd;

            if(rqAlignment.orientation().isForward())
            {
                adjSequenceStart = softClipSeqIndexStart + rqAlignment.sequenceStart();
                adjSequenceEnd = softClipSeqIndexStart + rqAlignment.sequenceEnd();
            }
            else
            {
                int newSequenceStart = (softClipLength - 1) - (rqAlignment.rawSequenceEnd() - 1);
                int newSequenceEnd = newSequenceStart + rqAlignment.segmentLength() - 1;

                adjSequenceStart = softClipSeqIndexStart + newSequenceStart;
                adjSequenceEnd = softClipSeqIndexStart + newSequenceEnd;
            }

            rqAlignment.setRequeriedSequenceCoords(adjSequenceStart, adjSequenceEnd);
        }

        if(isLeftClip)
        {
            for(AlignData alignment : alignments)
            {
                newAlignments.add(alignment);
            }
        }
        else
        {
            for(int i = 0; i < alignments.size(); i++)
            {
                newAlignments.add(i, alignments.get(i));
            }
        }

        return newAlignments;
    }

    private List<AlignData> requerySupplementaryAlignments(
            final AssemblyAlignment assemblyAlignment, final List<AlignData> alignments, final List<AlignData> requeriedAlignments)
    {
        // re-align supplementaries to get a more reliable map quality
        if(alignments.stream().noneMatch(x -> x.isSupplementary()) || alignments.stream().allMatch(x -> x.isSupplementary()))
        {
            return alignments;
        }

        List<AlignData> newAlignments = Lists.newArrayList();

        for(AlignData alignData : alignments)
        {
            if(!alignData.isSupplementary())
            {
                newAlignments.add(alignData);
                continue;
            }

            requeriedAlignments.add(alignData);
            newAlignments.addAll(requerySupplementaryAlignment(assemblyAlignment, alignData));
        }

        return newAlignments;
    }

    private List<AlignData> requerySupplementaryAlignment(final AssemblyAlignment assemblyAlignment, final AlignData alignData)
    {
        ++mRequeriedSuppCount;

        String fullSequence = assemblyAlignment.fullSequence();

        alignData.setFullSequenceData(fullSequence, assemblyAlignment.fullSequenceLength());

        String alignmentSequence = fullSequence.substring(alignData.sequenceStart(), alignData.sequenceEnd() + 1);

        List<BwaMemAlignment> requeryBwaAlignments = mAligner.alignSequence(alignmentSequence.getBytes());

        List<AlignData> requeryAlignments = requeryBwaAlignments.stream()
                .map(x -> AlignData.fromRef(x, mConfig.RefGenVersion))
                .filter(x -> x != null).collect(Collectors.toList());

        List<AlignData> convertedAlignments = Lists.newArrayList();

        for(AlignData rqAlignment : requeryAlignments)
        {
            rqAlignment.setFullSequenceData(alignmentSequence, alignmentSequence.length());

            // eg:
            // alignData = {AlignData@3240} "10:2543491-2543563 72S73M fwd seq(72-145 adj=72-144) score(58) flags(2048) mapQual(55 align=73 adj=73)"
            // rqAlignment = {AlignData@3246} "10:2543809-2543878 3S70M fwd seq(3-73 adj=3-72) score(65) flags(0) mapQual(17 align=70 adj=70)"

            AlignData convertedAlignment = new AlignData(
                    rqAlignment.refLocation(),
                    rqAlignment.rawSequenceStart(),
                    rqAlignment.rawSequenceEnd(),
                    rqAlignment.mapQual(), rqAlignment.score(), rqAlignment.flags(), rqAlignment.cigar(), rqAlignment.nMatches(),
                    rqAlignment.xaTag(), rqAlignment.mdTag());

            // restore values to be in terms of the original sequence
            int rqSeqOffsetStart = rqAlignment.sequenceStart();
            int adjSequenceStart = alignData.sequenceStart() + rqSeqOffsetStart;
            int rqSeqOffsetEnd = alignmentSequence.length() - 1 - rqAlignment.sequenceEnd();
            int adjSequenceEnd = alignData.sequenceEnd() - rqSeqOffsetEnd;
            convertedAlignment.setRequeriedSequenceCoords(adjSequenceStart, adjSequenceEnd);

            convertedAlignment.setSoftClipLengths(
                    convertedAlignment.leftSoftClipLength() + alignData.leftSoftClipLength(),
                    convertedAlignment.rightSoftClipLength() + alignData.rightSoftClipLength());

            convertedAlignments.add(convertedAlignment);

            alignData.markDroppedOnRequery();
        }

        return convertedAlignments;
    }

    private static boolean isWeakSingleReadExtensionAssembly(final AssemblyAlignment assemblyAlignment)
    {
        if(assemblyAlignment.breakends().isEmpty())
            return false;

        if(assemblyAlignment.assemblies().size() != 1)
        {
            // a single assembly including if it locally aligned
            if(!isLocalIndelAssembly(assemblyAlignment))
                return false;
        }

        // check not a short indel
        if(assemblyAlignment.breakends().size() == 2)
        {
            Breakend breakend = assemblyAlignment.breakends().get(0);
            if(isShortLocalDelDupIns(breakend.svType(), breakend.svLength()))
            {
                return false;
            }
        }

        JunctionAssembly assembly = assemblyAlignment.assemblies().get(0);

        if(assembly.stats().MaxExtBaseMatchCount < ASSEMBLY_UNLINKED_WEAK_ASSEMBLY_EXTENSION_LENGTH)
        {
            return false;
        }

        if(assembly.hasLineSequence() || assembly.junction().DiscordantOnly || assembly.junction().indelBased()) // simple split junction
        {
            return false;
        }

        boolean hasSuspectExtension = hasSuspectExtension(assemblyAlignment);

        boolean testMinorityExtensions = false;

        if(isIllumina())
        {
            if(assembly.stats().SoftClipSecondMaxLength < ASSEMBLY_MIN_SOFT_CLIP_LENGTH) // only 1 read above the min length
                return true;

            if(!hasSuspectExtension)
                return false;
        }
        else
        {
            testMinorityExtensions = true;
        }

        if(testMinorityExtensions)
        {
            // see if there is a minority of longer reads being
            return hasLongerMinorityExtensions(assembly.support(), assembly.junction().Orient, hasSuspectExtension);
        }

        return false;
    }

    @VisibleForTesting
    public static boolean hasSuspectExtension(final AssemblyAlignment assemblyAlignment)
    {
        Breakend breakend = assemblyAlignment.breakends().get(0);

        // the variant is a BND or has over 100k length
        if(breakend.svType() != StructuralVariantType.BND && breakend.svLength() < WEAK_ASSEMBLY_LONG_LENGTH)
            return false;

        // has 6+ bases of inexact homology AND the sample has a discordant fragment rate above 2%
        if(breakend.Homology.length() >= INV_SHORT_MAX_HOMOLOGY_HIGHER && AssemblyConfig.SampleDiscordantRate > WEAK_ASSEMBLY_DISC_RATE)
            return true;

        // OR has at least 16 bases of insert sequence
        if(breakend.InsertedBases.length() >= WEAK_ASSEMBLY_MIN_INS_LENGTH)
            return true;

        // OR is adjacent to a ref repeat of 7+ bp)
        for(JunctionAssembly assembly : assemblyAlignment.assemblies())
        {
            String refBases = assembly.formRefBaseSequence();
            String adjacentRefBases;
            int minRefLength = WEAK_ASSEMBLY_MIN_REF_REPEAT * 3; // so will only allow repeats of 3mers or less

            if(refBases.length() <= minRefLength)
            {
                adjacentRefBases = refBases;
            }
            else if(assembly.isForwardJunction())
            {
                adjacentRefBases = refBases.substring(refBases.length() - minRefLength);
            }
            else
            {
                adjacentRefBases = refBases.substring(0, minRefLength);
            }

            List<RepeatInfo> refRepeats = RepeatInfo.findRepeats(adjacentRefBases.getBytes(), WEAK_ASSEMBLY_MIN_REF_REPEAT);

            if(assembly.isForwardJunction())
                return refRepeats.stream().anyMatch(x -> x.lastIndex() >= adjacentRefBases.length() - x.repeatLength());
            else
                return refRepeats.stream().anyMatch(x -> x.Index == 0);
        }

        return false;
    }

    @VisibleForTesting
    public static boolean hasLongerMinorityExtensions(
            final List<SupportRead> reads, final Orientation assemblyOrientation, boolean useIncreasedLimits)
    {
        // in these cases:
        // - add one extra read to maxLongReadCount (still capped at 3)
        // - lower WEAK_ASSEMBLY_UNPAIRED_LONG_EXT_FACTOR to 2 for all comparisons before the last allowed one (so if maxLongReadCount = 2,
        //   the first check uses 2, the second uses 2.5)

        List<Integer> extensionLengths = Lists.newArrayListWithCapacity(reads.size());

        for(SupportRead read : reads)
        {
            if(read.type() != SupportType.JUNCTION)
                continue;

            int extensionLength = read.junctionExtensionLength(assemblyOrientation);

            extensionLengths.add(extensionLength);
        }

        Collections.sort(extensionLengths, Collections.reverseOrder());

        int minLongReadCount = 1;

        if(extensionLengths.size() > 2)
        {
            int secondLongestLength = extensionLengths.get(1);

            List<SupportRead> topReads = reads.stream()
                    .filter(x -> x.type() == SupportType.JUNCTION && x.junctionExtensionLength(assemblyOrientation) >= secondLongestLength)
                    .collect(Collectors.toList());

            SupportRead read1 = topReads.get(0);
            SupportRead read2 = topReads.get(1);

            if(readsCloseMatch(read1, read2)
            && (read1.mediumQualCount() >= SBX_MEDIUM_QUAL_DESYNC_COUNT || read2.mediumQualCount() >= SBX_MEDIUM_QUAL_DESYNC_COUNT))
            {
                minLongReadCount = 2;
            }
        }

        double maxLongReadCount = max(reads.size() / WEAK_ASSEMBLY_UNPAIRED_READ_FACTOR, minLongReadCount);

        if(useIncreasedLimits)
            ++maxLongReadCount;

        maxLongReadCount = min(maxLongReadCount, WEAK_ASSEMBLY_UNPAIRED_MAX_READS);

        double upperFactor = WEAK_ASSEMBLY_UNPAIRED_LONG_EXT_FACTOR;
        double lowerFactor = useIncreasedLimits ? WEAK_ASSEMBLY_UNPAIRED_LONG_EXT_FACTOR_LOWER : upperFactor;

        for(int i = 0; i < extensionLengths.size() - 1; ++i)
        {
            int extLength = extensionLengths.get(i);
            int nextLength = extensionLengths.get(i + 1);

            int readCount = i + 1;

            if(readCount > maxLongReadCount)
                return false;

            double longExtFactor = readCount == maxLongReadCount ? upperFactor : lowerFactor;

            if(extLength > nextLength * longExtFactor)
                return true;
        }

        return false;
    }

    private static boolean readsCloseMatch(final SupportRead read1, final SupportRead read2)
    {
        if(read1.orientation() != read2.orientation())
            return false;

        int threePrime1 = read1.orientation().isForward() ? read1.untrimmedEnd() : read1.untrimmedStart();
        int threePrime2 = read2.orientation().isForward() ? read2.untrimmedEnd() : read2.untrimmedStart();

        return abs(threePrime1 - threePrime2) <= SBX_PRIME_POSITION_RANGE_THRESHOLD;
    }
}
