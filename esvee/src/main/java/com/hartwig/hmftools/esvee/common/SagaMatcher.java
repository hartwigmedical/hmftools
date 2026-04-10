package com.hartwig.hmftools.esvee.common;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.min;

import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_ALIGN_SCORE_MIN_RATIO;
import static com.hartwig.hmftools.esvee.common.SvConstants.SAGA_LOCATION_MATCH_DISTANCE;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.Stream;

import com.hartwig.hmftools.esvee.assembly.alignment.BwaAligner;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class SagaMatcher
{
    private final SagaResource mSagaResource;
    private final BwaAligner mAligner;

    public SagaMatcher(final SagaResource sagaResource)
    {
        mSagaResource = sagaResource;
        // TODO: set params? min alignment score? output multiple alignments?
        mAligner = new BwaAligner(sagaResource.bwaIndexImagePath());
    }

    public SagaResource.Variant matchByLocation(final String chromosome, int position)
    {
        Map<String, List<SagaResource.IndexedBreakend>> breakends = mSagaResource.searchableBreakends();
        List<SagaResource.IndexedBreakend> chrBreakends = breakends.get(chromosome);
        return matchByLocationOnChromosome(position, chrBreakends);
    }

    private SagaResource.Variant matchByLocationOnChromosome(int position, List<SagaResource.IndexedBreakend> chrBreakends)
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
            return mSagaResource.getVariantById(bestVariant);
        }
    }

    public SagaResource.Variant matchBySequence(final byte[] sequence, int junctionOffset1, @Nullable Integer junctionOffset2)
    {
        if(!(junctionOffset1 >= 0 && junctionOffset1 <= sequence.length))
        {
            throw new IllegalArgumentException("Invalid junctionOffset1: " + junctionOffset1);
        }
        if(junctionOffset2 != null && !(junctionOffset2 >= 0 && junctionOffset2 <= sequence.length))
        {
            throw new IllegalArgumentException("Invalid junctionOffset2: " + junctionOffset2);
        }

        List<BwaMemAlignment> alignments = mAligner.alignSequence(sequence);
        return matchFromAlignments(alignments, junctionOffset1, junctionOffset2).orElse(null);
    }

    private Optional<SagaResource.Variant> matchFromAlignments(List<BwaMemAlignment> alignments, int junctionOffset1,
            @Nullable Integer junctionOffset2)
    {
        Stream<AlignmentCandidate> candidates = alignments.stream()
                .map(alignment -> evaluateAlignmentCandidate(alignment, junctionOffset1, junctionOffset2))
                .filter(Objects::nonNull)
                .sorted();
        return candidates.findFirst().map(candidate -> mSagaResource.getVariantByFastaLabel(candidate.variantFastaLabel()));
    }

    private record AlignmentCandidate(
            String variantFastaLabel,
            int alignScore,
            // [2][2][2]
            int[][][] junctionOverlaps)
            implements Comparable<AlignmentCandidate>
    {
        @Override
        public int compareTo(@NotNull final AlignmentCandidate other)
        {
            // Higher alignment score is better.
            return Integer.compare(-alignScore, -other.alignScore);
        }
    }

    private static AlignmentCandidate evaluateAlignmentCandidate(final BwaMemAlignment alignment, int junctionOffset1,
            @Nullable Integer junctionOffset2)
    {
        int alignLength1 = alignment.getRefEnd() - alignment.getRefStart();
        int alignLength2;   // TODO

        int minAlignScore = calcMinAlignScore(alignLength1, alignLength2);
        if(alignment.getAlignerScore() < minAlignScore)
        {
            return null;
        }

        // TODO: check junction overlap
    }

    private static int calcMinAlignScore(int alignLength1, int alignLength2)
    {
        return (int) ceil(min(alignLength1, alignLength2) * SAGA_ALIGN_SCORE_MIN_RATIO);
    }
}
