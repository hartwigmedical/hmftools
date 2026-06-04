package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_JUNCTION_INDEL_DISTANCE;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_JUNCTION_INDEL_MAX_LENGTH;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_JUNCTION_OVERLAP_MIN;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_LENGTH_MIN_RATIO;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_SCORE_MIN_BASELINE;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_SCORE_MIN_RATIO;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.bwa.BwaMemAlignParams;
import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.bwa.IBwaMemAligner;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;

// Matches assemblies to SAGA variants by base sequence similarity.
public class SagaSequenceMatcher
{
    private final Config mConfig;
    private final IBwaMemAligner mAligner;
    private final Map<Integer, SagaAssembly> mAssembliesByContigId;

    // Do not instantiate directly - use SagaMatcherFactory.
    public SagaSequenceMatcher(final Config config, final IBwaMemAligner aligner, final Map<Integer, SagaAssembly> assembliesByContigId)
    {
        mConfig = config;
        mAligner = aligner;
        mAssembliesByContigId = assembliesByContigId;
    }

    public record Config(
            double alignLengthRatioMin,
            int alignScoreMin,
            double alignScoreRatioMin,
            int junctionOverlapMin
    )
    {
        static Config DEFAULT = new Config(
                SAGA_ALIGN_LENGTH_MIN_RATIO,
                SAGA_ALIGN_SCORE_MIN_BASELINE,
                SAGA_ALIGN_SCORE_MIN_RATIO,
                SAGA_ALIGN_JUNCTION_OVERLAP_MIN
        );
    }

    // junctionOffsets are the indices just after the junction in sequence.
    // E.g.
    // sequence = RRRJJJRRR
    // junctionOffset[0] = 3
    // junctionOffset[1] = 6
    @Nullable
    public SagaMatchBySequence matchBySequence(final byte[] sequence, List<Integer> junctionOffsets)
    {
        List<BwaMemAlignment> alignments = mAligner.alignSequence(sequence);
        return matchFromAlignments(sequence, alignments, junctionOffsets);
    }

    @Nullable
    public SagaMatchBySequence matchBySequence(final String sequence, List<Integer> junctionOffsets)
    {
        return matchBySequence(sequence.getBytes(), junctionOffsets);
    }

    @Nullable
    private SagaMatchBySequence matchFromAlignments(final byte[] sequence, List<BwaMemAlignment> alignments,
            List<Integer> junctionOffsets)
    {
        Stream<AlignmentCandidate> candidates = alignments.stream()
                .map(alignment -> evaluateAlignmentCandidate(sequence, alignment, junctionOffsets))
                .filter(Objects::nonNull)
                .sorted(new AlignmentCandidateComparator(sequence.length));
        return candidates.findFirst()
                .map(candidate -> new SagaMatchBySequence(candidate.variant(), candidate.cigar(), candidate.alignScore()))
                .orElse(null);
    }

    private record AlignmentCandidate(
            SagaAssembly assembly,
            Cigar cigar,
            int alignScore,
            List<List<Integer>> esveeJunctionOverlaps,
            List<List<Integer>> sagaJunctionOverlaps)
    {
        public SagaVariant variant()
        {
            return assembly.variant();
        }
    }

    private record AlignmentCandidateComparator(int sequenceLength) implements Comparator<AlignmentCandidate>
    {
        @Override
        public int compare(final AlignmentCandidate o1, final AlignmentCandidate o2)
        {
            // Higher alignment score is better.
            int res = Integer.compare(-o1.alignScore, -o2.alignScore);
            if(res != 0)
            {
                return res;
            }

            // Then prefer the variant whose length is closer to the ESVEE assembly length (rarely occurs but can help tie-break different alleles).
            int distance1 = abs(o1.assembly().assemblyLength() - sequenceLength);
            int distance2 = abs(o2.assembly().assemblyLength() - sequenceLength);
            return Integer.compare(distance1, distance2);
        }
    }

    private AlignmentCandidate evaluateAlignmentCandidate(
            final byte[] sequence, final BwaMemAlignment alignment, final List<Integer> seqJunctionOffsets)
    {
        if(alignment.getRefId() < 0 || SamRecordUtils.isFlagSet(alignment.getSamFlag(), SAMFlag.READ_UNMAPPED))
        {
            // Not a real alignment, means the query is not aligned.
            return null;
        }

        boolean isForward = !SamRecordUtils.isFlagSet(alignment.getSamFlag(), SAMFlag.READ_REVERSE_STRAND);

        // ESVEE assemblies are implied as forward strand and so are SAGA assemblies, so there should never be a reverse strand match.
        if(!isForward)
        {
            return null;
        }

        SagaAssembly sagaAssembly = getAssemblyByContigId(alignment.getRefId());

        Cigar cigar = cigarFromStr(alignment.getCigar());
        int leftClip = leftClipLength(cigar);
        int rightClip = rightClipLength(cigar);

        int sagaAlignLength = alignment.getRefEnd() - alignment.getRefStart();
        int seqAlignLength = sequence.length - (leftClip + rightClip);

        // Only match if the majority of the sequences align, where possible.
        // I.e. don't allow a small subsequence match where much more sequence was possible to match.
        if(!isAlignLengthOk(seqAlignLength, sagaAlignLength, sequence.length, sagaAssembly.assemblyLength()))
        {
            return null;
        }

        int minAlignScore = calcMinAlignScore(sagaAlignLength, seqAlignLength);
        // Low alignment score means low sequence similarity, so we think this is not a good match.
        if(alignment.getAlignerScore() < minAlignScore)
        {
            return null;
        }

        int startClip = isForward ? leftClip : rightClip;
        int endClip = isForward ? rightClip : leftClip;
        int seqStart = startClip;
        int seqEnd = sequence.length - endClip;

        // The novel sequence created by the variant junction needs to be included in the alignment.
        // Otherwise, we could align onto the surrounding ref bases which could be similar, but the variant is different.
        List<List<Integer>> seqJunctionOverlaps =
                seqJunctionOffsets.stream().map(offset -> calcJunctionOverlap(seqStart, seqEnd, offset)).toList();
        if(!isJunctionOverlapOk(seqJunctionOverlaps))
        {
            return null;
        }

        List<List<Integer>> sagaJunctionOverlaps = sagaAssembly.junctionOffsets()
                .stream()
                .map(offset -> calcJunctionOverlap(alignment.getRefStart(), alignment.getRefEnd(), offset))
                .toList();
        if(!isJunctionOverlapOk(sagaJunctionOverlaps))
        {
            return null;
        }

        // Ensure there are no large indels near the junction.
        // This would mean the junction sequence is different, despite the rest of the sequence aligning.
        // Potential FIXME: limit to only the matched junctions.
        if(!checkIndelsNotNearJunctions(cigar, alignment.getRefStart(), seqJunctionOffsets, sagaAssembly.junctionOffsets()))
        {
            return null;
        }

        return new AlignmentCandidate(sagaAssembly, cigar, alignment.getAlignerScore(), seqJunctionOverlaps, sagaJunctionOverlaps);
    }

    private boolean isAlignLengthOk(int seqAlignLength, int sagaAlignLength, int seqAssemblyLength, int sagaAssemblyLength)
    {
        int minAlignLength = calcMinAlignLength(seqAssemblyLength, sagaAssemblyLength);
        return min(seqAlignLength, sagaAlignLength) >= minAlignLength;
    }

    private int calcMinAlignLength(int assemblyLength1, int assemblyLength2)
    {
        return (int) ceil(min(assemblyLength1, assemblyLength2) * mConfig.alignLengthRatioMin);
    }

    private int calcMinAlignScore(int alignLength1, int alignLength2)
    {
        return (int) ceil(min(alignLength1, alignLength2) * mConfig.alignScoreRatioMin);
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

    private boolean isJunctionOverlapOk(final List<List<Integer>> overlaps)
    {
        // Require the alignment to cover the junction +/- some tolerance bases.
        // At least one junction must be covered in this manner. We don't require all the junctions to be covered because we could be
        // matching just one side of the variant (e.g. for junction assembly, or an SGL).
        int minOverlap = mConfig.junctionOverlapMin;
        return overlaps.stream()
                .anyMatch(junctionOverlaps -> junctionOverlaps.stream().allMatch(overlap -> overlap >= minOverlap));
    }

    private record ElementWithPositions(
            CigarElement element,
            int readStart,
            int readEnd,    // Exclusive
            int refStart,
            int refEnd      // Exclusive
    )
    {
        public CigarOperator operator()
        {
            return element.getOperator();
        }
    }

    private static List<ElementWithPositions> getElementPositions(final Cigar cigar, int refStart)
    {
        List<ElementWithPositions> result = new ArrayList<>();
        int readIndex = 0;
        int refIndex = refStart;
        for(CigarElement element : cigar.getCigarElements())
        {
            CigarOperator operator = element.getOperator();
            int nextReadIndex = readIndex;
            if(operator.consumesReadBases() || operator.isClipping())
            {
                nextReadIndex += element.getLength();
            }
            int nextRefIndex = refIndex;
            if(operator.consumesReferenceBases())
            {
                nextRefIndex += element.getLength();
            }
            result.add(new ElementWithPositions(element, readIndex, nextReadIndex, refIndex, nextRefIndex));

            readIndex = nextReadIndex;
            refIndex = nextRefIndex;
        }
        return result;
    }

    private static boolean checkIndelsNotNearJunctions(
            final Cigar cigar, int sagaStart, final List<Integer> seqJunctionOffsets, final List<Integer> sagaJunctionOffsets)
    {
        return getElementPositions(cigar, sagaStart).stream()
                .filter(e -> e.operator().isIndel())
                .allMatch(element -> checkIndelOnJunctions(element.element(), element.readStart(), element.readEnd(), seqJunctionOffsets)
                        && checkIndelOnJunctions(element.element(), element.refStart(), element.refEnd(), sagaJunctionOffsets));
    }

    private static boolean checkIndelOnJunctions(final CigarElement element, int start, int end, final List<Integer> junctions)
    {
        if(element.getLength() <= SAGA_ALIGN_JUNCTION_INDEL_MAX_LENGTH)
        {
            // Allow very small indels.
            return true;
        }
        // Ensure there are no large indels near the junctions.
        return junctions.stream().noneMatch(j -> isJunctionNearIndel(j, start, end));
    }

    private static boolean isJunctionNearIndel(int junctionOffset, int indelStart, int indelEnd)
    {
        boolean insideIndel = junctionOffset >= indelStart && junctionOffset < indelEnd;
        boolean nearIndelStart = abs(junctionOffset - indelStart) < SAGA_ALIGN_JUNCTION_INDEL_DISTANCE;
        boolean nearIndelEnd = abs(junctionOffset - indelEnd) < SAGA_ALIGN_JUNCTION_INDEL_DISTANCE;
        return insideIndel || nearIndelStart || nearIndelEnd;
    }

    public static BwaMemAligner.Params createAlignerParams(final Config matchConfig)
    {
        // single-threaded in since Esvee makes these calls within threads already, and don't batch since phased assemblies are aligned
        // one at a time, also in a threaded context
        return new BwaMemAligner.Params(
                // Note: default params may be too strict. If a case where BWA drops an alignment of a valid match, should relax these.
                BwaMemAlignParams.DEFAULT,
                true,
                matchConfig.alignScoreMin(),
                1,
                null
        );
    }

    private SagaAssembly getAssemblyByContigId(int contigId)
    {
        SagaAssembly assembly = mAssembliesByContigId.get(contigId);
        if(assembly == null)
        {
            throw new IllegalArgumentException("No SAGA variant with contig ID: " + contigId);
        }
        else
        {
            return assembly;
        }
    }
}
