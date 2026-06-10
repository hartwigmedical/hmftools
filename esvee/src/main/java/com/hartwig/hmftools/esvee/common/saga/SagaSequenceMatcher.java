package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;

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
    private final SagaSequenceMatcherConfig mConfig;
    private final IBwaMemAligner mAligner;
    private final Map<Integer, SagaAssembly> mAssembliesByContigId;

    // Do not instantiate directly - use SagaMatcherFactory.
    public SagaSequenceMatcher(final SagaSequenceMatcherConfig config, final IBwaMemAligner aligner,
            final Map<Integer, SagaAssembly> assembliesByContigId)
    {
        mConfig = config;
        mAligner = aligner;
        mAssembliesByContigId = assembliesByContigId;
    }

    private record MatchArguments(
            byte[] query,
            List<Integer> junctionOffsets,
            boolean lowerJunctionOverlap
    )
    {
    }

    // junctionOffsets are the indices just after the junction in sequence.
    // E.g.
    // sequence = RRRJJJRRR
    // junctionOffset[0] = 3
    // junctionOffset[1] = 6
    @Nullable
    public SagaMatchBySequence matchBySequence(final byte[] sequence, final List<Integer> junctionOffsets, boolean lowerJunctionOverlap)
    {
        MatchArguments args = new MatchArguments(sequence, junctionOffsets, lowerJunctionOverlap);
        List<BwaMemAlignment> alignments = mAligner.alignSequence(sequence);
        return matchFromAlignments(args, alignments);
    }

    @Nullable
    public SagaMatchBySequence matchBySequence(final String sequence, final List<Integer> junctionOffsets, boolean lowerJunctionOverlap)
    {
        return matchBySequence(sequence.getBytes(), junctionOffsets, lowerJunctionOverlap);
    }

    @Nullable
    private SagaMatchBySequence matchFromAlignments(final MatchArguments args, final List<BwaMemAlignment> rawAlignments)
    {
        List<SagaSequenceMatchCandidate> candidates = rawAlignments.stream()
                .map(alignment -> wrapAlignment(alignment, args.query().length))
                .filter(Objects::nonNull)
                .map(alignment -> evaluateAlignment(alignment, args))
                .sorted(new MatchCandidateComparator(args))
                .toList();
        return candidates.stream()
                .filter(SagaSequenceMatchCandidate::isAccepted)
                .findFirst()
                .map(candidate -> new SagaMatchBySequence(candidate.variant(), candidate.cigar(), candidate.alignScore()))
                .orElse(null);
    }

    private record MatchCandidateComparator(MatchArguments args) implements Comparator<SagaSequenceMatchCandidate>
    {
        @Override
        public int compare(final SagaSequenceMatchCandidate o1, final SagaSequenceMatchCandidate o2)
        {
            if(o1.isAccepted() != o2.isAccepted())
            {
                return o1.isAccepted() ? -1 : 1;
            }

            // Higher alignment score is better.
            int res = Integer.compare(-o1.alignScore(), -o2.alignScore());
            if(res != 0)
            {
                return res;
            }

            // Then prefer the variant whose length is closer to the ESVEE assembly length (rarely occurs but can help tie-break different alleles).
            int distance1 = abs(o1.sagaAssembly().assemblyLength() - args.query().length);
            int distance2 = abs(o2.sagaAssembly().assemblyLength() - args.query().length);
            return Integer.compare(distance1, distance2);
        }
    }

    @Nullable
    private SagaAlignment wrapAlignment(final BwaMemAlignment alignment, int queryLength)
    {
        if(alignment.getRefId() < 0 || SamRecordUtils.isFlagSet(alignment.getSamFlag(), SAMFlag.READ_UNMAPPED))
        {
            // Not a real alignment, means the query is not aligned.
            return null;
        }

        SagaAssembly sagaAssembly = getAssemblyByContigId(alignment.getRefId());
        Cigar cigar = cigarFromStr(alignment.getCigar());
        return new SagaAlignment(alignment, cigar, queryLength, sagaAssembly);
    }

    private SagaSequenceMatchCandidate evaluateAlignment(final SagaAlignment alignment, final MatchArguments args)
    {
        List<String> filters = new ArrayList<>();

        // ESVEE assemblies are implied as forward strand and so are SAGA assemblies, so there should never be a reverse strand match.
        if(!alignment.isForward())
        {
            filters.add("reverse_strand");
        }

        // Only match if the majority of the sequences align, where possible.
        // I.e. don't allow a small subsequence match where much more sequence was possible to match.
        if(!isAlignLengthOk(alignment))
        {
            filters.add("aligned_length");
        }

        // Low alignment score means low sequence similarity, so we think this is not a good match.
        if(!isAlignScoreOk(alignment))
        {
            filters.add("align_score");
        }

        // The novel sequence created by the variant junction needs to be included in the alignment.
        // Otherwise, we could align onto the surrounding ref bases which could be similar, but the variant is different.
        List<List<Integer>> seqJunctionOverlaps = args.junctionOffsets().stream()
                .map(offset -> calcJunctionOverlap(alignment.queryStart(), alignment.queryEnd(), offset))
                .toList();
        if(!isJunctionOverlapOk(seqJunctionOverlaps, args.lowerJunctionOverlap()))
        {
            filters.add("esvee_junction_overlap");
        }

        List<List<Integer>> sagaJunctionOverlaps = alignment.sagaAssembly().junctionOffsets()
                .stream()
                .map(offset -> calcJunctionOverlap(alignment.sagaStart(), alignment.sagaEnd(), offset))
                .toList();
        if(!isJunctionOverlapOk(sagaJunctionOverlaps, args.lowerJunctionOverlap()))
        {
            filters.add("saga_junction_overlap");
        }

        // Ensure there are no large indels near the junction.
        // This would mean the junction sequence is different, despite the rest of the sequence aligning.
        // Potential FIXME: limit to only the matched junctions.
        if(!checkIndelsNotNearJunctions(alignment.cigar(), alignment.sagaStart(), args.junctionOffsets(), alignment.sagaAssembly()
                .junctionOffsets()))
        {
            filters.add("indel_near_junction");
        }

        return new SagaSequenceMatchCandidate(alignment, seqJunctionOverlaps, sagaJunctionOverlaps, filters);
    }

    private int calcMaxUnalignedBases(int assemblyLength1, int assemblyLength2)
    {
        return (int) ceil(min(assemblyLength1, assemblyLength2) * (1 - mConfig.alignLengthRatioMin()));
    }

    private boolean isAlignLengthOk(final SagaAlignment alignment)
    {
        // Checking this way instead of just looking at the alignment lengths, because for a single-sided assembly there may not be enough
        // bases to align through to 80% of the other side, even though this is a legit match.
        // It's better to look at how many bases could've been aligned on either side of the aligned subsequence and see if this is too many.
        int unaligned = alignment.leftUnaligned() + alignment.rightUnaligned();
        int maxUnaligned = calcMaxUnalignedBases(alignment.queryLength(), alignment.sagaLength());
        return unaligned <= maxUnaligned;
    }

    private int calcMinAlignScore(int alignLength1, int alignLength2)
    {
        return (int) ceil(min(alignLength1, alignLength2) * mConfig.alignScoreRatioMin());
    }

    private boolean isAlignScoreOk(final SagaAlignment alignment)
    {
        int minAlignScore = calcMinAlignScore(alignment.queryAlignLength(), alignment.sagaAlignLength());
        return alignment.alignScore() >= minAlignScore;
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

    private boolean isJunctionOverlapOk(final List<List<Integer>> overlaps, boolean lowerJunctionOverlap)
    {
        // Require the alignment to cover the junction +/- some tolerance bases.
        // At least one junction must be covered in this manner. We don't require all the junctions to be covered because we could be
        // matching just one side of the variant (e.g. for junction assembly, or an SGL).
        int minOverlap = calcJunctionOverlapMin(lowerJunctionOverlap);
        return overlaps.stream()
                .anyMatch(junctionOverlaps -> junctionOverlaps.stream().allMatch(overlap -> overlap >= minOverlap));
    }

    private int calcJunctionOverlapMin(boolean lowerJunctionOverlap)
    {
        // For LINE assemblies, the extension is allowed to be shorter.
        return lowerJunctionOverlap ? mConfig.junctionOverlapMinLower() : mConfig.junctionOverlapMin();
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

    private boolean checkIndelsNotNearJunctions(
            final Cigar cigar, int sagaStart, final List<Integer> seqJunctionOffsets, final List<Integer> sagaJunctionOffsets)
    {
        return getElementPositions(cigar, sagaStart).stream()
                .filter(e -> e.operator().isIndel())
                .allMatch(element -> checkIndelOnJunctions(element.element(), element.readStart(), element.readEnd(), seqJunctionOffsets)
                        && checkIndelOnJunctions(element.element(), element.refStart(), element.refEnd(), sagaJunctionOffsets));
    }

    private boolean checkIndelOnJunctions(final CigarElement element, int start, int end, final List<Integer> junctions)
    {
        if(element.getLength() <= mConfig.junctionIndelLengthMax())
        {
            // Allow very small indels.
            return true;
        }
        // Ensure there are no large indels near the junctions.
        return junctions.stream().noneMatch(j -> isJunctionNearIndel(j, start, end));
    }

    private boolean isJunctionNearIndel(int junctionOffset, int indelStart, int indelEnd)
    {
        boolean insideIndel = junctionOffset >= indelStart && junctionOffset < indelEnd;
        boolean nearIndelStart = abs(junctionOffset - indelStart) < mConfig.junctionIndelDistance();
        boolean nearIndelEnd = abs(junctionOffset - indelEnd) < mConfig.junctionIndelDistance();
        return insideIndel || nearIndelStart || nearIndelEnd;
    }

    public static BwaMemAligner.Params createAlignerParams(final SagaSequenceMatcherConfig matchConfig)
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
