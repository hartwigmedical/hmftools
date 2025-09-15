package com.hartwig.hmftools.redux.bqr;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class RefSequence
{
    public final int Start;
    public final int End;
    public final byte[] Bases;
    public final boolean IsValid;

    private static final int BUFFER = 5;

    public RefSequence(final String chromosome, final int posStart, final int posEnd, final RefGenomeInterface refGenome)
    {
        final int sequenceEnd = refGenome.getChromosomeLength(chromosome);

        Start = max(1, posStart - BUFFER);
        End = min(sequenceEnd, posEnd + BUFFER);
        Bases = refGenome.getBases(chromosome, Start, End);
        IsValid = true;
    }

    public int index(int position)
    {
        return position - Start;
    }
    public int position(int index) { return Start + index; }

    public boolean positionWithinBounds(int position)
    {
        return position >= Start + BUFFER && position <= End - BUFFER;
    }

    public boolean beforeStart(int position)
    {
        return position < Start + BUFFER;
    }

    public boolean afterEnd(int position)
    {
        return position > End - BUFFER;
    }

    public int length()
    {
        return Bases.length;
    }

    public byte base(int position)
    {
        return Bases[index(position)];
    }

    public byte[] trinucleotideContext(int position)
    {
        int index = index(position);
        return new byte[] { Bases[index - 1], Bases[index], Bases[index + 1] };
    }

    public void populateTrinucleotideContext(int position, final byte[] trinucleotideContext)
    {
        int index = index(position);
        trinucleotideContext[0] = Bases[index - 1];
        trinucleotideContext[1] = Bases[index];
        trinucleotideContext[2] = Bases[index + 1];
    }

    public String toString()
    {
        return format("%d-%d len(%d)", Start, End, length());
    }

}
