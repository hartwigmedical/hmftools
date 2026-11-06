package com.hartwig.hmftools.viridian.profiling.align;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRUS_GENOME_ORIGIN_CLIP_TOLERANCE;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import htsjdk.samtools.Cigar;

// One alignment of a read to a virus genome, with just the data required for our analysis, in a nice format.
// This class is quite optimised for memory, otherwise it can explode since we keep multiple alignments per read.
public record ViralReadAlignment(
        ReadId readId,
        VirusGenome genome,
        int alignmentStart,
        Cigar cigar,
        int alignerScore,
        // NM + clipped bases: how poorly the virus genome explains the read
        int divergence,
        int[] mismatchPositions)
{
    public ViralReadAlignment
    {
        if(alignmentStart < 1)
        {
            throw new IllegalArgumentException("Invalid alignmentStart: " + alignmentStart);
        }
        if(cigar.isEmpty())
        {
            throw new IllegalArgumentException("Invalid cigar");
        }
        if(alignerScore < 0)
        {
            throw new IllegalArgumentException("Invalid alignerScore: " + alignerScore);
        }
        if(divergence < leftClipLength(cigar) + rightClipLength(cigar))
        {
            throw new IllegalArgumentException("Invalid divergence: " + divergence);
        }
        int alignmentEnd = alignmentStart + cigar.getReferenceLength() - 1;
        if(alignmentEnd > genome.length())
        {
            throw new IllegalArgumentException(String.format(
                    "Alignment end %d past virus genome %s length %d", alignmentEnd, genome.contigName(), genome.length()));
        }
    }

    public int alignmentEnd()
    {
        return alignmentStart + cigar.getReferenceLength() - 1;
    }

    public int leftClip()
    {
        return leftClipLength(cigar);
    }

    public int rightClip()
    {
        return rightClipLength(cigar);
    }

    public List<MatchedInterval> matchedIntervals()
    {
        return MatchedInterval.from(alignmentStart, cigar, mismatchPositions);
    }

    // The clipped bases project past a virus genome end, so the read straddles the circular genome's linearization
    // origin: those bases wrap to the other end and cannot align linearly. An artifact, not real divergence.
    public boolean clipsOverGenomeEnd()
    {
        int genomeLength = genome.length();
        return alignmentStart - leftClip() < 1 - VIRUS_GENOME_ORIGIN_CLIP_TOLERANCE
                || alignmentEnd() + rightClip() > genomeLength + VIRUS_GENOME_ORIGIN_CLIP_TOLERANCE;
    }

    @Override
    public boolean equals(Object other)
    {
        return other instanceof ViralReadAlignment that
                && alignmentStart == that.alignmentStart
                && alignerScore == that.alignerScore
                && divergence == that.divergence
                && readId.equals(that.readId)
                && genome.equals(that.genome)
                && cigar.equals(that.cigar)
                && Arrays.equals(mismatchPositions, that.mismatchPositions);
    }

    @Override
    public int hashCode()
    {
        return 31 * Objects.hash(readId, genome, alignmentStart, cigar, alignerScore, divergence)
                + Arrays.hashCode(mismatchPositions);
    }

    // Sort primarily by divergence because we are mostly interested in per-base difference.
    // Aligner score is not the best fit.
    public static final Comparator<ViralReadAlignment> BEST_FIRST_TOTAL_ORDER = Comparator
            .comparingInt(ViralReadAlignment::divergence)
            .thenComparing(ViralReadAlignment::alignerScore, Comparator.reverseOrder())
            .thenComparingInt(ViralReadAlignment::alignmentStart)
            .thenComparing(alignment -> alignment.cigar().toString());
}
