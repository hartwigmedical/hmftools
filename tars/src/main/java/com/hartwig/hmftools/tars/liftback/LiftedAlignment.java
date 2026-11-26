package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public record LiftedAlignment(
        BasePosition liftedPosition,
        String liftedCigar,
        int numMismatches,
        boolean fromTxContig,
        boolean forwardStrand,
        // +1 forward / -1 reverse for tx-contig alignments; 0 when no transcript strand is known.
        @Nullable Orientation transcriptStrand,
        List<Integer> mergedSupplementaryIndices,
        List<ChrBaseRegion> mergedSupplementaryIntrons,
        int mergedSupplementaryMapQuality,
        List<BaseRegion> alignedBlocks,
        // set when an overhang collapse made this alt a duplicate placement; excluded from the XA tag.
        boolean dropped,
        @Nullable Integer genomicScore)
{
    public LiftedAlignment
    {
        mergedSupplementaryIndices = List.copyOf(mergedSupplementaryIndices);
        mergedSupplementaryIntrons = List.copyOf(mergedSupplementaryIntrons);
    }

    public LiftedAlignment(
            final String liftedChromosome, final int liftedPos, final String liftedCigar, final int numMismatches,
            final boolean fromTxContig, final boolean forwardStrand, @Nullable final Orientation transcriptStrand)
    {
        this(
                new BasePosition(liftedChromosome, liftedPos), liftedCigar, numMismatches, fromTxContig, forwardStrand,
                transcriptStrand, Collections.emptyList(), Collections.emptyList(), -1,
                alignedBlocks(liftedPos, liftedCigar), false, null);
    }

    // GenomicScore is left uncomputed: it described the previous cigar
    public LiftedAlignment withLiftedCigar(final int liftedPos, final String liftedCigar)
    {
        return new LiftedAlignment(
                new BasePosition(liftedPosition.Chromosome, liftedPos), liftedCigar, numMismatches,
                fromTxContig, forwardStrand, transcriptStrand,
                mergedSupplementaryIndices, mergedSupplementaryIntrons, mergedSupplementaryMapQuality,
                alignedBlocks(liftedPos, liftedCigar), dropped, null);
    }

    public LiftedAlignment withSupplementaryMerge(
            final int liftedPos, final String liftedCigar, final List<Integer> mergedSupplementaryIndices,
            final List<ChrBaseRegion> mergedSupplementaryIntrons, final int mergedSupplementaryMapQuality,
            @Nullable final Orientation spliceStrand)
    {
        Orientation mergedTranscriptStrand = transcriptStrand != null ? transcriptStrand : spliceStrand;
        return new LiftedAlignment(
                new BasePosition(liftedPosition.Chromosome, liftedPos), liftedCigar, numMismatches,
                fromTxContig, forwardStrand, mergedTranscriptStrand,
                mergedSupplementaryIndices, mergedSupplementaryIntrons, mergedSupplementaryMapQuality,
                alignedBlocks(liftedPos, liftedCigar), false, null);
    }

    public LiftedAlignment withTranscriptStrand(final Orientation transcriptStrand)
    {
        return new LiftedAlignment(
                liftedPosition, liftedCigar, numMismatches, fromTxContig, forwardStrand, transcriptStrand,
                mergedSupplementaryIndices, mergedSupplementaryIntrons, mergedSupplementaryMapQuality,
                alignedBlocks, dropped, genomicScore);
    }

    public LiftedAlignment withDropped()
    {
        return new LiftedAlignment(
                liftedPosition, liftedCigar, numMismatches, fromTxContig, forwardStrand, transcriptStrand,
                mergedSupplementaryIndices, mergedSupplementaryIntrons, mergedSupplementaryMapQuality,
                alignedBlocks, true, genomicScore);
    }

    public LiftedAlignment withGenomicScore(final int genomicScore)
    {
        return new LiftedAlignment(
                liftedPosition, liftedCigar, numMismatches, fromTxContig, forwardStrand, transcriptStrand,
                mergedSupplementaryIndices, mergedSupplementaryIntrons, mergedSupplementaryMapQuality,
                alignedBlocks, dropped, genomicScore);
    }

    public boolean hasSupplementaryMerge()
    {
        return !mergedSupplementaryIndices.isEmpty();
    }

    // bwa XA entry "chrom,<sign>pos,cigar,NM;", the position sign carrying the strand.
    public String toXaEntry()
    {
        return liftedPosition.Chromosome + ',' + (forwardStrand ? '+' : '-') + liftedPosition.Position + ','
                + liftedCigar + ',' + numMismatches + ';';
    }

    public AlignmentKey key()
    {
        return AlignmentKey.from(this);
    }

    public int alignedEnd()
    {
        return liftedPosition.Position + CigarUtils.calcCigarAlignedLength(liftedCigar) - 1;
    }

    // Minimum genomic gap between aligned reference blocks. N operations split blocks because they are RNA introns,
    // not aligned sequence; treating the outer CIGAR span as one interval would make a mate inside an intron distance 0.
    public int alignedBlockDistance(final LiftedAlignment other)
    {
        if(!liftedPosition.Chromosome.equals(other.liftedPosition.Chromosome))
        {
            return Integer.MAX_VALUE;
        }

        int closest = Integer.MAX_VALUE;
        for(BaseRegion block : alignedBlocks)
        {
            for(BaseRegion otherBlock : other.alignedBlocks)
            {
                int distance;
                if(block.start() > otherBlock.end())
                {
                    distance = block.start() - otherBlock.end();
                }
                else if(otherBlock.start() > block.end())
                {
                    distance = otherBlock.start() - block.end();
                }
                else
                {
                    return 0;
                }
                closest = Math.min(closest, distance);
            }
        }
        return closest;
    }

    private static List<BaseRegion> alignedBlocks(final int liftedPos, final String liftedCigar)
    {
        List<BaseRegion> blocks = new ArrayList<>();
        int referencePosition = liftedPos;
        int blockStart = -1;
        int blockEnd = -1;

        for(CigarElement element : CigarUtils.cigarFromStr(liftedCigar).getCigarElements())
        {
            CigarOperator operator = element.getOperator();
            if(operator == CigarOperator.N)
            {
                if(blockStart >= 0)
                {
                    blocks.add(new BaseRegion(blockStart, blockEnd));
                    blockStart = -1;
                }
                referencePosition += element.getLength();
            }
            else if(operator.consumesReferenceBases())
            {
                if(blockStart < 0)
                {
                    blockStart = referencePosition;
                }
                referencePosition += element.getLength();
                blockEnd = referencePosition - 1;
            }
        }
        if(blockStart >= 0)
        {
            blocks.add(new BaseRegion(blockStart, blockEnd));
        }
        return blocks;
    }

    // Same-locus test shared by the locus count and the XA build: the two must agree. False against a null other, so an
    // absent primary collapses nothing.
    public boolean overlaps(final LiftedAlignment other)
    {
        if(!liftedPosition.Chromosome.equals(other.liftedPosition.Chromosome))
        {
            return false;
        }
        return positionsOverlap(liftedPosition.Position, alignedEnd(), other.liftedPosition.Position, other.alignedEnd());
    }

    public boolean cigarHasN()
    {
        return liftedCigar.indexOf('N') >= 0;
    }

    public boolean cigarHasSoftClip()
    {
        return liftedCigar.indexOf('S') >= 0;
    }
}
