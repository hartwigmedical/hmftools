package com.hartwig.hmftools.pave.transval;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * The position of a codon in a transcript together with the exon in which
 * it starts and the exon in which it ends, which is null if the codon
 * is contained in the first exon. (The shortest codon in the human genome
 * has length 4, so a codon cqn never require three exons.)
 */
public class CodonRegions
{
    final public int CodonStart;
    @NotNull final public ChrBaseRegion FirstExon;
    @Nullable
    final public ChrBaseRegion SecondExon;
    final boolean IsPositiveReadStrand;
    final int numberOfBaseInSecondExon;

    public CodonRegions(final int codonStart, @NotNull final ChrBaseRegion firstExon, @Nullable final ChrBaseRegion secondExon)
    {
        this(codonStart, firstExon, secondExon, true);
    }

    public CodonRegions(final int position, @NotNull final ChrBaseRegion firstExon, @Nullable final ChrBaseRegion secondExon,
            final boolean isPositiveReadStrand)
    {
        Preconditions.checkArgument(position >= firstExon.start());
        Preconditions.checkArgument(position <= firstExon.end());

        this.CodonStart = position;
        this.FirstExon = firstExon;
        this.SecondExon = secondExon;
        this.IsPositiveReadStrand = isPositiveReadStrand;

        numberOfBaseInSecondExon = IsPositiveReadStrand ? CodonStart + 2 - FirstExon.end() : FirstExon.start() + 2 - CodonStart;
    }

    public boolean codonIsInSingleExon()
    {
        return numberOfBaseInSecondExon  < 1;
    }

    @NotNull public String retrieveCodon(final RefGenomeInterface refGenome)
    {
        return IsPositiveReadStrand ? codonForPositiveRead(refGenome) : codonForNegativeRead(refGenome);
    }

    @NotNull private String codonForPositiveRead(final RefGenomeInterface refGenome)
    {
        if (numberOfBaseInSecondExon < 1)
        {
            return refGenome.getBaseString(FirstExon.Chromosome, CodonStart, CodonStart + 2);
        }
        String part1 = refGenome.getBaseString(FirstExon.Chromosome, CodonStart, FirstExon.end());
        String part2 = refGenome.getBaseString(FirstExon.Chromosome, SecondExon.start(), SecondExon.start() + numberOfBaseInSecondExon - 1);
        return part1 + part2;
    }

    @NotNull private String codonForNegativeRead(final RefGenomeInterface refGenome)
    {
        if (numberOfBaseInSecondExon < 1)
        {
            String positiveStrandBases = refGenome.getBaseString(FirstExon.Chromosome, CodonStart - 2, CodonStart);
            return Nucleotides.reverseComplementBases(positiveStrandBases);
        }
        String part1 = Nucleotides.reverseComplementBases(refGenome.getBaseString(FirstExon.Chromosome, FirstExon.start(), CodonStart));
        String part2 = Nucleotides.reverseComplementBases(refGenome.getBaseString(FirstExon.Chromosome, SecondExon.end() - numberOfBaseInSecondExon + 1, SecondExon.end()));
        return part1 + part2;
    }
}
