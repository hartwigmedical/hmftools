package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_JUNCTION_OVERLAP_MIN;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_SCORE_MIN_BASELINE;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_SCORE_MIN_RATIO;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_LOCATION_MATCH_DISTANCE;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.bwa.BwaMemAlignParams;
import com.hartwig.hmftools.common.bwa.BwaMemAligner;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMSequenceDictionary;

public class SagaMatcher
{
    private final SagaResource mSagaResource;
    private final BwaMemAligner mAligner;
    private final SAMSequenceDictionary mSagaDict;

    public SagaMatcher(final SagaResource sagaResource)
    {
        mSagaResource = sagaResource;
        mAligner = new BwaMemAligner(sagaResource.bwaIndexImagePath(), makeAlignerParams());
        mSagaDict = sagaResource.samDict();
    }

    public record MatchByLocation(
            SagaResource.Variant variant,
            int distance
    )
    {
    }

    @Nullable
    public MatchByLocation matchByLocation(final String chromosome, int position)
    {
        Map<String, List<SagaResource.IndexedBreakend>> breakends = mSagaResource.searchableBreakends();
        List<SagaResource.IndexedBreakend> chrBreakends = breakends.get(chromosome);
        return matchByLocationOnChromosome(position, chrBreakends);
    }

    @Nullable
    private MatchByLocation matchByLocationOnChromosome(int position, List<SagaResource.IndexedBreakend> chrBreakends)
    {
        String bestVariant = null;
        int bestDistance = SAGA_LOCATION_MATCH_DISTANCE + 1;
        if(chrBreakends != null)
        {
            // TODO: can binary search for better performance
            for(SagaResource.IndexedBreakend breakend : chrBreakends)
            {
                int distance = abs(breakend.position().Position - position);
                if(distance < bestDistance)
                {
                    bestVariant = breakend.variantId();
                    bestDistance = distance;
                }
            }
        }
        if(bestVariant == null)
        {
            return null;

        }
        else
        {
            SagaResource.Variant variant = mSagaResource.getVariantById(bestVariant);
            return new MatchByLocation(variant, bestDistance);
        }
    }

    public record MatchBySequence(
            SagaResource.Variant variant,
            Cigar cigar,
            int alignScore
    )
    {
    }

    // junctionOffsets are the indices just after the junction in sequence.
    // E.g.
    // sequence = RRRJJJRRR
    // junctionOffset[0] = 3
    // junctionOffset[1] = 6
    @Nullable
    public MatchBySequence matchBySequence(final byte[] sequence, List<Integer> junctionOffsets)
    {
        List<BwaMemAlignment> alignments = mAligner.alignSequence(sequence);
        return matchFromAlignments(sequence, alignments, junctionOffsets);
    }

    @Nullable
    private MatchBySequence matchFromAlignments(final byte[] sequence, List<BwaMemAlignment> alignments,
            List<Integer> junctionOffsets)
    {
        Stream<AlignmentCandidate> candidates = alignments.stream()
                .map(alignment -> evaluateAlignmentCandidate(sequence, alignment, junctionOffsets))
                .filter(Objects::nonNull)
                .sorted();
        return candidates.findFirst()
                .map(candidate -> new MatchBySequence(candidate.variant(), candidate.cigar(), candidate.alignScore()))
                .orElse(null);
    }

    private record AlignmentCandidate(
            SagaResource.Variant variant,
            Cigar cigar,
            int alignScore,
            List<List<Integer>> esveeJunctionOverlaps,
            List<List<Integer>> sagaJunctionOverlaps)
            implements Comparable<AlignmentCandidate>
    {
        @Override
        public int compareTo(@NotNull final AlignmentCandidate other)
        {
            // Higher alignment score is better.
            return Integer.compare(-alignScore, -other.alignScore);
        }
    }

    private AlignmentCandidate evaluateAlignmentCandidate(final byte[] sequence, final BwaMemAlignment alignment,
            List<Integer> seqJunctionOffsets)
    {
        if(alignment.getRefId() < 0 || SamRecordUtils.isFlagSet(alignment.getSamFlag(), SAMFlag.READ_UNMAPPED))
        {
            // Not a real alignment, means the query is not aligned.
            return null;
        }

        Cigar cigar = cigarFromStr(alignment.getCigar());
        int leftClip = leftClipLength(cigar);
        int rightClip = rightClipLength(cigar);

        int sagaAlignLength = alignment.getRefEnd() - alignment.getRefStart();
        int seqAlignLength = sequence.length - (leftClip + rightClip);
        int minAlignScore = calcMinAlignScore(sagaAlignLength, seqAlignLength);
        // Low alignment score means low sequence similarity, so we think this is not a good match.
        if(alignment.getAlignerScore() < minAlignScore)
        {
            return null;
        }

        boolean isForward = !SamRecordUtils.isFlagSet(alignment.getSamFlag(), SAMFlag.READ_REVERSE_STRAND);
        int startClip = isForward ? leftClip : rightClip;
        int endClip = isForward ? rightClip : leftClip;
        int seqStart = startClip;
        int seqEnd = sequence.length - endClip;

        // The novel sequence created by the variant junction needs to be included in the alignment.
        // Otherwise we could align onto the surrounding ref bases which could be similar, but the variant is not the same.
        List<List<Integer>> seqJunctionOverlaps =
                seqJunctionOffsets.stream().map(offset -> calcJunctionOverlap(seqStart, seqEnd, offset)).toList();
        if(!isJunctionOverlapOk(seqJunctionOverlaps))
        {
            return null;
        }

        String refContig = mSagaDict.getSequence(alignment.getRefId()).getSequenceName();
        SagaResource.AssemblyMetadata sagaAssembly = mSagaResource.getAssemblyByFastaLabel(refContig);

        List<List<Integer>> sagaJunctionOverlaps = sagaAssembly.junctionOffsets()
                .stream()
                .map(offset -> calcJunctionOverlap(alignment.getRefStart(), alignment.getRefEnd(), offset))
                .toList();
        if(!isJunctionOverlapOk(sagaJunctionOverlaps))
        {
            return null;
        }

        return new AlignmentCandidate(sagaAssembly.variant(), cigar, alignment.getAlignerScore(), seqJunctionOverlaps, sagaJunctionOverlaps);
    }

    private static int calcMinAlignScore(int alignLength1, int alignLength2)
    {
        return (int) ceil(min(alignLength1, alignLength2) * SAGA_ALIGN_SCORE_MIN_RATIO);
    }

    private static List<Integer> calcJunctionOverlap(int alignStart, int alignEnd, int junctionOffset)
    {
        // Example to help maths.
        // 012 345 678
        // RRR JJJ RRR
        // junctionOffsets[0] = 3
        // junctionOffsets[1] = 6
        return List.of(junctionOffset - alignStart, alignEnd - junctionOffset);
    }

    private static boolean isJunctionOverlapOk(List<List<Integer>> overlaps)
    {
        // Require the alignment to cover the junction +/- some tolerance bases.
        // At least one junction must be covered in this manner. We don't require all the junctions to be covered because we could be
        // matching just one side of the variant (e.g. for junction assembly, or an SGL).
        return overlaps.stream()
                .anyMatch(junctionOverlaps -> junctionOverlaps.stream().allMatch(overlap -> overlap >= SAGA_ALIGN_JUNCTION_OVERLAP_MIN));
    }

    private static BwaMemAligner.Params makeAlignerParams()
    {
        return new BwaMemAligner.Params(
                // TODO? may have to relax some params? but measure performance hit
                BwaMemAlignParams.DEFAULT,
                true,
                SAGA_ALIGN_SCORE_MIN_BASELINE,
                // Single threaded because ESVEE is already multithreaded appropriately.
                1,
                // Don't care about batching because we align 1 sequence at a time.
                null
        );
    }
}
