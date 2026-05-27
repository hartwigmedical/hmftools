package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.Objects.requireNonNull;
import static java.util.function.UnaryOperator.identity;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_JUNCTION_OVERLAP_MIN;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_LENGTH_MIN_RATIO;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_SCORE_MIN_BASELINE;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_SCORE_MIN_RATIO;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_LOCATION_MATCH_DISTANCE;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.bwa.BwaMemAlignParams;
import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.bwa.IBwaMemAligner;
import com.hartwig.hmftools.common.region.BasePosition;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMSequenceRecord;

public class SagaMatcher
{
    private final Map<String, SagaResource.AssemblyMetadata> mAssembliesByVariantId;
    // Map from contig to list of breakends sorted by position ascending.
    private final Map<String, List<IndexedBreakend>> mSearchableBreakends;
    private final IBwaMemAligner mAligner;
    // So we can quickly look up a variant from the alignment contig ID.
    private final Map<Integer, String> mContigIdToVariantId;
    private final MatchConfig mMatchConfig;

    private SagaMatcher(final List<SagaResource.AssemblyMetadata> assemblies, final IBwaMemAligner aligner,
            final Map<Integer, String> contigIdToName, final MatchConfig matchConfig)
    {
        mAssembliesByVariantId = assemblies.stream()
                .collect(Collectors.toMap(SagaResource.AssemblyMetadata::variantId, identity()));

        mSearchableBreakends = assemblies.stream()
                .flatMap(assembly ->
                        assembly.variant().breakends().map(breakend ->
                                new IndexedBreakend(breakend.position(), assembly.variantId())))
                .collect(Collectors.groupingBy(IndexedBreakend::chromosome));
        // Sort by position so they can be binary-searched.
        mSearchableBreakends.forEach((chr, breakends) -> breakends.sort(Comparator.comparing(IndexedBreakend::position)));

        // Precompute a map of the contig ID integer to the variant ID, so we can look it up from a BwaMemAlignment.
        Map<String, String> contigNameToVariantId = assemblies.stream()
                .collect(Collectors.toMap(SagaResource.AssemblyMetadata::fastaLabel, SagaResource.AssemblyMetadata::variantId));
        mContigIdToVariantId = contigIdToName.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, entry -> requireNonNull(contigNameToVariantId.get(entry.getValue()))));

        mAligner = aligner;

        mMatchConfig = matchConfig;
    }

    private SagaMatcher(final SagaResource sagaResource, final MatchConfig matchConfig)
    {
        this(
                sagaResource.assemblies(),
                new BwaMemAligner(sagaResource.bwaIndexImagePath(), createAlignerParams(matchConfig)),
                sagaResource.samDict().getSequences()
                        .stream()
                        .collect(Collectors.toMap(SAMSequenceRecord::getSequenceIndex, SAMSequenceRecord::getSequenceName)),
                matchConfig);
    }

    public SagaMatcher(final SagaResource sagaResource)
    {
        this(sagaResource, MatchConfig.DEFAULT);
    }

    public record MatchConfig(
            int locationDistanceMax,
            double alignLengthRatioMin,
            int alignScoreMin,
            double alignScoreRatioMin,
            int junctionOverlapMin
    )
    {
        static MatchConfig DEFAULT = new MatchConfig(
                SAGA_LOCATION_MATCH_DISTANCE,
                SAGA_ALIGN_LENGTH_MIN_RATIO,
                SAGA_ALIGN_SCORE_MIN_BASELINE,
                SAGA_ALIGN_SCORE_MIN_RATIO,
                SAGA_ALIGN_JUNCTION_OVERLAP_MIN
        );
    }

    private record IndexedBreakend(
            BasePosition basePosition,
            String variantId
    )
    {
        public String chromosome()
        {
            return basePosition.Chromosome;
        }

        public int position()
        {
            return basePosition.Position;
        }
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
        Map<String, List<IndexedBreakend>> breakends = mSearchableBreakends;
        List<IndexedBreakend> chrBreakends = breakends.get(chromosome);
        return matchByLocationOnChromosome(position, chrBreakends);
    }

    @Nullable
    private MatchByLocation matchByLocationOnChromosome(int position, List<IndexedBreakend> chrBreakends)
    {
        IndexedBreakend bestBreakend = null;
        int bestDistance = mMatchConfig.locationDistanceMax + 1;
        if(chrBreakends != null)
        {
            IndexedBreakend target =
                    new IndexedBreakend(new BasePosition(chrBreakends.get(0).chromosome(), position), "");
            int index = Collections.binarySearch(chrBreakends, target, Comparator.comparingInt(IndexedBreakend::position));
            if(index >= 0)
            {
                // Exact position found.
                bestBreakend = chrBreakends.get(index);
                bestDistance = 0;
            }
            else
            {
                // Exact position not found, so need to look at adjacent indices.
                index = -(index + 1);
                int startIndex = max(index - 1, 0);
                int endIndex = min(index + 1, chrBreakends.size() - 1);
                for(int i = startIndex; i <= endIndex; i++)
                {
                    IndexedBreakend breakend = chrBreakends.get(i);
                    int distance = abs(breakend.position() - position);
                    if(distance < bestDistance)
                    {
                        bestBreakend = breakend;
                        bestDistance = distance;
                    }
                }
            }
        }
        if(bestBreakend == null)
        {
            return null;

        }
        else
        {
            SagaResource.Variant variant = getAssemblyByVariantId(bestBreakend.variantId()).variant();
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
    public MatchBySequence matchBySequence(final String sequence, List<Integer> junctionOffsets)
    {
        return matchBySequence(sequence.getBytes(), junctionOffsets);
    }

    @Nullable
    private MatchBySequence matchFromAlignments(final byte[] sequence, List<BwaMemAlignment> alignments,
            List<Integer> junctionOffsets)
    {
        Stream<AlignmentCandidate> candidates = alignments.stream()
                .map(alignment -> evaluateAlignmentCandidate(sequence, alignment, junctionOffsets))
                .filter(Objects::nonNull)
                .sorted(new AlignmentCandidateComparator(sequence.length));
        return candidates.findFirst()
                .map(candidate -> new MatchBySequence(candidate.variant(), candidate.cigar(), candidate.alignScore()))
                .orElse(null);
    }

    private record AlignmentCandidate(
            SagaResource.AssemblyMetadata assembly,
            Cigar cigar,
            int alignScore,
            List<List<Integer>> esveeJunctionOverlaps,
            List<List<Integer>> sagaJunctionOverlaps)
    {
        public SagaResource.Variant variant()
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

    private AlignmentCandidate evaluateAlignmentCandidate(final byte[] sequence, final BwaMemAlignment alignment,
            List<Integer> seqJunctionOffsets)
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

        SagaResource.AssemblyMetadata sagaAssembly = getAssemblyByContigId(alignment.getRefId());

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

        return new AlignmentCandidate(sagaAssembly, cigar, alignment.getAlignerScore(), seqJunctionOverlaps, sagaJunctionOverlaps);
    }

    private boolean isAlignLengthOk(int seqAlignLength, int sagaAlignLength, int seqAssemblyLength, int sagaAssemblyLength)
    {
        int minAlignLength = calcMinAlignLength(seqAssemblyLength, sagaAssemblyLength);
        return min(seqAlignLength, sagaAlignLength) >= minAlignLength;
    }

    private int calcMinAlignLength(int assemblyLength1, int assemblyLength2)
    {
        return (int) ceil(min(assemblyLength1, assemblyLength2) * mMatchConfig.alignLengthRatioMin);
    }

    private int calcMinAlignScore(int alignLength1, int alignLength2)
    {
        return (int) ceil(min(alignLength1, alignLength2) * mMatchConfig.alignScoreRatioMin);
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
        int minOverlap = mMatchConfig.junctionOverlapMin;
        return overlaps.stream()
                .anyMatch(junctionOverlaps -> junctionOverlaps.stream().allMatch(overlap -> overlap >= minOverlap));
    }

    private static BwaMemAligner.Params createAlignerParams(final MatchConfig matchConfig)
    {
        return new BwaMemAligner.Params(
                // TODO? may have to relax some params? but measure performance hit
                BwaMemAlignParams.DEFAULT,
                true,
                matchConfig.alignScoreMin(),
                // Single threaded because ESVEE is already multithreaded appropriately.
                1,
                // Don't care about batching because we align 1 sequence at a time.
                null
        );
    }

    private SagaResource.AssemblyMetadata getAssemblyByVariantId(final String variantId)
    {
        SagaResource.AssemblyMetadata assembly = mAssembliesByVariantId.get(variantId);
        if(assembly == null)
        {
            throw new IllegalArgumentException("No SAGA variant with ID: " + variantId);
        }
        else
        {
            return assembly;
        }
    }

    private SagaResource.AssemblyMetadata getAssemblyByContigId(int contigId)
    {
        String variantId = mContigIdToVariantId.get(contigId);
        if(variantId == null)
        {
            throw new IllegalArgumentException("No SAGA variant with contig ID: " + contigId);
        }
        else
        {
            return getAssemblyByVariantId(variantId);
        }
    }
}
