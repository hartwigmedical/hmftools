package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.min;
import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
import com.hartwig.hmftools.common.bwa.BwaMemAlignParams;
import com.hartwig.hmftools.common.bwa.BwaMemAlignerConfig;
import com.hartwig.hmftools.common.bwa.IBwaMemAligner;
import com.hartwig.hmftools.common.utils.IntPair;

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
            List<SagaJunctionInfo> junctions,
            boolean lowerJunctionOverlap,
            boolean allowReverseStrand
    )
    {
    }

    @Nullable
    public SagaSequenceMatch matchBySequence(final byte[] sequence, final List<SagaJunctionInfo> junctions, boolean lowerJunctionOverlap,
            boolean allowReverseStrand)
    {
        MatchArguments args = new MatchArguments(sequence, junctions, lowerJunctionOverlap, allowReverseStrand);
        return matchBySequenceImpl(args);
    }

    @Nullable
    private SagaSequenceMatch matchBySequenceImpl(final MatchArguments args)
    {
        List<BwaMemAlignment> alignments = mAligner.alignSequence(args.query);
        return matchFromAlignments(args, alignments);
    }

    @Nullable
    private SagaSequenceMatch matchFromAlignments(final MatchArguments args, final List<BwaMemAlignment> rawAlignments)
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
                .map(candidate -> new SagaSequenceMatch(candidate.alignment()))
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
            int distance1 = abs(o1.sagaAssembly().length() - args.query().length);
            int distance2 = abs(o2.sagaAssembly().length() - args.query().length);
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
        SagaAlignment sagaAlignment = new SagaAlignment(alignment, cigar, queryLength, sagaAssembly);
        sagaAlignment.validate();
        return sagaAlignment;
    }

    private SagaSequenceMatchCandidate evaluateAlignment(final SagaAlignment alignment, final MatchArguments args)
    {
        Set<SagaSequenceMatchCandidateFilter> filters = new TreeSet<>();

        applyAlignmentFilters(alignment, args, filters);

        List<SagaJunctionMatchInfo> queryJunctionMatches = new ArrayList<>();
        List<SagaJunctionMatchInfo> sagaJunctionMatches = new ArrayList<>();
        applyJunctionFilters(args, alignment, queryJunctionMatches, sagaJunctionMatches, filters);

        return new SagaSequenceMatchCandidate(alignment, queryJunctionMatches, sagaJunctionMatches, filters);
    }

    private void applyAlignmentFilters(final SagaAlignment alignment, final MatchArguments args,
            Set<SagaSequenceMatchCandidateFilter> filters)
    {
        // For junction assemblies, the sequence may match to SAGA in reverse if it's part of an INV.
        // For phased assemblies, expect that it matches to the forward strand because ESVEE assemblies are always forward strand.
        if(!args.allowReverseStrand() && !alignment.isForward())
        {
            filters.add(SagaSequenceMatchCandidateFilter.REVERSE_STRAND);
        }

        // Only match if the majority of the sequences align, where possible.
        // I.e. don't allow a small subsequence match where much more sequence was possible to match.
        if(!checkAlignLength(alignment))
        {
            filters.add(SagaSequenceMatchCandidateFilter.ALIGNED_LENGTH);
        }

        // Low alignment score means low sequence similarity, so we think this is not a good match.
        if(!checkAlignScore(alignment))
        {
            filters.add(SagaSequenceMatchCandidateFilter.ALIGN_SCORE);
        }
    }

    private int calcMaxUnalignedBases(int assemblyLength1, int assemblyLength2)
    {
        return (int) ceil(min(assemblyLength1, assemblyLength2) * (1 - mConfig.alignLengthRatioMin()));
    }

    private boolean checkAlignLength(final SagaAlignment alignment)
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

    private boolean checkAlignScore(final SagaAlignment alignment)
    {
        int minAlignScore = calcMinAlignScore(alignment.queryAlignLength(), alignment.sagaAlignLength());
        return alignment.alignScore() >= minAlignScore;
    }

    private SagaJunctionMatchInfo evaluateJunctionAlignment(final SagaAlignment alignment, final List<CigarElemWithPos> cigarElements,
            final SagaJunctionInfo junctionInfo, boolean isQuery)
    {
        int alignStart = isQuery ? alignment.queryStart() : alignment.sagaStart();
        int alignEnd = isQuery ? alignment.queryEnd() : alignment.sagaEnd();
        IntPair junctionOverlap = calcJunctionOverlap(alignStart, alignEnd, junctionInfo.assemblyOffset());
        boolean indelNearby = junctionHasIndelNearby(junctionInfo, cigarElements, isQuery);
        return new SagaJunctionMatchInfo(junctionInfo, junctionOverlap.left, junctionOverlap.right, indelNearby);
    }

    private static IntPair calcJunctionOverlap(int alignStart, int alignEnd, int junctionOffset)
    {
        // Example to help maths.
        // 012 345 678
        // RRR JJJ RRR
        // junctionOffsets[0] = 3
        // junctionOffsets[1] = 6
        return new IntPair(junctionOffset - alignStart, alignEnd - junctionOffset);
    }

    private boolean checkJunctionOverlap(final SagaJunctionMatchInfo junctionMatchInfo, boolean lowerJunctionOverlap)
    {
        // Require the alignment to cover the junction +/- some tolerance bases.
        // Otherwise, we could align onto the surrounding ref bases which could be similar, but the variant is different.
        int minOverlap = calcJunctionOverlapMin(lowerJunctionOverlap);
        return junctionMatchInfo.alignOverlapLeft() >= minOverlap && junctionMatchInfo.alignOverlapRight() >= minOverlap;
    }

    private int calcJunctionOverlapMin(boolean lowerJunctionOverlap)
    {
        // For LINE assemblies, the extension is allowed to be shorter.
        return lowerJunctionOverlap ? mConfig.junctionOverlapMinLower() : mConfig.junctionOverlapMin();
    }

    private record CigarElemWithPos(
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

    private static List<CigarElemWithPos> getCigarElementPositions(final Cigar cigar, int refStart)
    {
        List<CigarElemWithPos> result = new ArrayList<>();
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
            result.add(new CigarElemWithPos(element, readIndex, nextReadIndex, refIndex, nextRefIndex));

            readIndex = nextReadIndex;
            refIndex = nextRefIndex;
        }
        return result;
    }

    private boolean junctionHasIndelNearby(final SagaJunctionInfo junctionInfo, final List<CigarElemWithPos> cigarElements, boolean isQuery)
    {
        return cigarElements.stream()
                .filter(e -> e.operator().isIndel() && !ignoreIndel(e.element()))
                .anyMatch(e -> isJunctionNearIndel(
                        junctionInfo.assemblyOffset(), isQuery ? e.readStart() : e.refStart(), isQuery ? e.readEnd() : e.refEnd()));
    }

    private boolean ignoreIndel(final CigarElement element)
    {
        // Allow very small indels.
        return element.getLength() <= mConfig.junctionIndelLengthMax();
    }

    private boolean isJunctionNearIndel(int junctionOffset, int indelStart, int indelEnd)
    {
        boolean insideIndel = junctionOffset >= indelStart && junctionOffset < indelEnd;
        boolean nearIndelStart = abs(junctionOffset - indelStart) < mConfig.junctionIndelDistance();
        boolean nearIndelEnd = abs(junctionOffset - indelEnd) < mConfig.junctionIndelDistance();
        return insideIndel || nearIndelStart || nearIndelEnd;
    }

    private List<SagaJunctionMatchInfo> filterJunctionMatches(final MatchArguments args,
            final List<SagaJunctionMatchInfo> junctionMatchInfos, boolean isQuery, Set<SagaSequenceMatchCandidateFilter> filters)
    {
        List<SagaJunctionMatchInfo> alignedJunctions = junctionMatchInfos.stream()
                .filter(j -> checkJunctionOverlap(j, args.lowerJunctionOverlap()))
                .toList();

        // Require at least 1 junction to be covered by the alignment.
        // Don't require all the junctions to be covered because we could be matching just one side of the variant (e.g. for junction assembly, or an SGL).
        if(alignedJunctions.isEmpty())
        {
            filters.add(isQuery
                    ? SagaSequenceMatchCandidateFilter.QUERY_JUNCTION_OVERLAP
                    : SagaSequenceMatchCandidateFilter.SAGA_JUNCTION_OVERLAP);
        }

        // Ensure there are no large indels near the junctions.
        // This would mean the junction sequence is different, despite the rest of the sequence aligning.
        // Potential FIXME: limit to only the matched junctions.
        if(junctionMatchInfos.stream().anyMatch(SagaJunctionMatchInfo::indelNearby))
        {
            filters.add(isQuery
                    ? SagaSequenceMatchCandidateFilter.QUERY_JUNCTION_INDEL
                    : SagaSequenceMatchCandidateFilter.SAGA_JUNCTION_INDEL);
            return emptyList();
        }
        else
        {
            return alignedJunctions;
        }
    }

    private void applyJunctionFilters(final MatchArguments args, final SagaAlignment alignment,
            List<SagaJunctionMatchInfo> queryJunctionMatches, List<SagaJunctionMatchInfo> sagaJunctionMatches,
            Set<SagaSequenceMatchCandidateFilter> filters)
    {
        List<CigarElemWithPos> cigarElements = getCigarElementPositions(alignment.cigar(), alignment.sagaStart());

        List<SagaJunctionMatchInfo> queryJunctionMatchInfos = args.junctions().stream()
                .map(j -> evaluateJunctionAlignment(alignment, cigarElements, j, true)).toList();
        List<SagaJunctionMatchInfo> sagaJunctionMatchInfos = alignment.sagaAssembly().junctions().stream()
                .map(j -> evaluateJunctionAlignment(alignment, cigarElements, j, false)).toList();

        queryJunctionMatches.addAll(filterJunctionMatches(args, queryJunctionMatchInfos, true, filters));
        sagaJunctionMatches.addAll(filterJunctionMatches(args, sagaJunctionMatchInfos, false, filters));
    }

    public static BwaMemAlignerConfig createAlignerConfig(final String indexImagePath, final SagaSequenceMatcherConfig matchConfig)
    {
        // single-threaded in since Esvee makes these calls within threads already, and don't batch since phased assemblies are aligned
        // one at a time, also in a threaded context
        return new BwaMemAlignerConfig(
                indexImagePath,
                createAlignParams(matchConfig),
                true,
                1,
                null
        );
    }

    private static BwaMemAlignParams createAlignParams(final SagaSequenceMatcherConfig matchConfig)
    {
        // Note: default params may be too strict. If a case where BWA drops an alignment of a valid match, should relax these.
        return BwaMemAlignParams.DEFAULT
                .withMinAlignScore(matchConfig.alignScoreMin());
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
