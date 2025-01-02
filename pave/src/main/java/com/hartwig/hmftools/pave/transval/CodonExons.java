package com.hartwig.hmftools.pave.transval;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * The position of a codon in a transcript together with the exon in which
 * it starts and the exon in which it ends, which is null if the codon
 * is contained in the first exon. (The shortest codon in the human genome
 * has length 4, so a codon cqn never require three exons.)
 */
public class CodonExons
{
    final public int codonStart;
    @NotNull final public ExonData FirstExon;
    @Nullable
    final public ExonData SecondExon;

    public CodonExons(final int codonStart, @NotNull final ExonData firstExon, @Nullable final ExonData secondExon)
    {
        Preconditions.checkArgument(codonStart >= firstExon.Start);
        Preconditions.checkArgument(codonStart <= firstExon.End);
//        Preconditions.checkArgument(secondExon == null || secondExon.Rank == firstExon.Rank + 1);
        Preconditions.checkArgument(secondExon != null || codonStart + 2 <= firstExon.End);

        this.codonStart = codonStart;
        this.FirstExon = firstExon;
        this.SecondExon = secondExon;
    }

    @NotNull public String retrieveCodon(final String chromosome, final RefGenomeInterface refGenome)
    {
        int numberOfBasesInSecondCodon = codonStart + 2 - FirstExon.End;
        if (numberOfBasesInSecondCodon < 1)
        {
            return refGenome.getBaseString(chromosome, codonStart, codonStart + 2);
        }
        String part1 = refGenome.getBaseString(chromosome, codonStart, FirstExon.End + 1);
        String part2 = refGenome.getBaseString(chromosome, SecondExon.Start, SecondExon.Start + numberOfBasesInSecondCodon);
        return part1 + part2;
    }
}
