package com.hartwig.hmftools.sage.common;

import static java.lang.Math.E;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.Arrays;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class RefSequence
{
    public final int Start;
    public final int End;
    public final byte[] Bases;

    private static final int BUFFER = 1000;

    public RefSequence(final ChrBaseRegion region, final RefGenomeInterface refGenome)
    {
        final int sequenceEnd = refGenome.getChromosomeLength(region.Chromosome);

        Start = max(1, region.start() - BUFFER);
        End = min(sequenceEnd, region.end() + BUFFER);
        Bases = refGenome.getBases(region.Chromosome, Start, End);
    }

    public RefSequence(final int startPosition, final byte[] refBases)
    {
        Start = startPosition;
        End = startPosition + refBases.length - 1;
        Bases = refBases;
    }

    public int index(int position)
    {
        return position - Start;
    }

    public int length() { return Bases.length; }

    public boolean containsPosition(int position)
    {
        int index = index(position);
        return index >= 0 && index < Bases.length;
    }

    public byte base(int position) { return Bases[index(position)]; }

    public byte[] baseRange(int posStart, int posEnd)
    {
        int indexStart = index(posStart);
        int indexEnd = index(posEnd);

        if(indexStart < 0 || indexEnd < indexStart || indexEnd >= Bases.length)
            return null;

        return Arrays.subsetArray(Bases, indexStart, indexEnd);
    }

    public String indexBases(int start, int end)
    {
        if(start < 0 || end >= Bases.length || start > end)
            return null;

        return new String(Bases, start, end - start + 1);
    }

    public String positionBases(int posStart, int posEnd)
    {
        int start = index(posStart);
        int length = posEnd - posStart;
        int end = start + length;

        if(start < 0 || end >= Bases.length || start > end)
            return null;

        return new String(Bases, start, length + 1);
    }

    public byte[] trinucleotideContext(int position)
    {
        int index = index(position);
        return new byte[] { Bases[index - 1], Bases[index], Bases[index + 1] };
    }

    public String toString() { return format("%d-%d len(%d)", Start, End, length()); }

    @VisibleForTesting
    public RefSequence(final ChrBaseRegion region, final ReferenceSequenceFile refGenome)
    {
        final int sequenceEnd = refGenome.getSequenceDictionary().getSequence(region.Chromosome).getSequenceLength();
        Start = max(1, region.start() - BUFFER);
        End = min(sequenceEnd, region.end() + BUFFER);
        Bases = refGenome.getSubsequenceAt(region.Chromosome, Start, End).getBases();
    }

    @VisibleForTesting
    public RefSequence(final ReferenceSequence sequence)
    {
        Start = sequence.getContigIndex() + 1;
        End = Start + sequence.getBases().length - 1;
        Bases = sequence.getBases();
    }
}
