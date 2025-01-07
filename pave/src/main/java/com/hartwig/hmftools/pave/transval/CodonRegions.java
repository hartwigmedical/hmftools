package com.hartwig.hmftools.pave.transval;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * A codon may be contained in a single exon or may be split across two exons.
 * (The shortest exon in the human genome has length 4,
 * so a codon cqn never require three exons.)
 * This class encapsulates the mapping of the start of a codon within
 * an exon and its possible extension into the next exon.
 */
public class CodonRegions
{
    final public int CodonStart;
    final private boolean IsPositiveReadStrand;
    @NotNull
    final private ChrBaseRegion RegionInFirstExon;
    @Nullable
    final private ChrBaseRegion RegionInSecondExon;

    public CodonRegions(final int codonStart, @NotNull final ChrBaseRegion firstExon, @Nullable final ChrBaseRegion secondExon)
    {
        this(codonStart, firstExon, secondExon, true);
    }

    public CodonRegions(final int position, @NotNull final ChrBaseRegion firstExon, @Nullable final ChrBaseRegion secondExon,
            final boolean isPositiveReadStrand)
    {
        Preconditions.checkArgument(position >= firstExon.start());
        Preconditions.checkArgument(position <= firstExon.end());
        Preconditions.checkArgument(secondExon == null || secondExon.Chromosome.equals(firstExon.Chromosome));

        this.CodonStart = position;
        this.IsPositiveReadStrand = isPositiveReadStrand;

        if(IsPositiveReadStrand)
        {
            int numberOfBaseInSecondExon = Math.max(CodonStart + 2 - firstExon.end(), 0);
            if(numberOfBaseInSecondExon == 0)
            {
                RegionInFirstExon = new ChrBaseRegion(firstExon.chromosome(), CodonStart, CodonStart + 2);
                RegionInSecondExon = null;
            }
            else
            {
                assert secondExon != null;
                RegionInFirstExon = new ChrBaseRegion(firstExon.chromosome(), CodonStart, firstExon.end());
                RegionInSecondExon = new ChrBaseRegion(firstExon.chromosome(), secondExon.start(), secondExon.start() + numberOfBaseInSecondExon - 1);
            }
        }
        else
        {
            int numberOfBaseInSecondExon = Math.max(firstExon.start() + 2 - CodonStart, 0);
            if(numberOfBaseInSecondExon == 0)
            {
                RegionInFirstExon = new ChrBaseRegion(firstExon.Chromosome, CodonStart - 2, CodonStart);
                RegionInSecondExon = null;
            }
            else
            {
                assert secondExon != null;
                RegionInFirstExon = new ChrBaseRegion(firstExon.Chromosome, firstExon.start(), CodonStart);
                RegionInSecondExon = new ChrBaseRegion(firstExon.Chromosome, secondExon.end() - numberOfBaseInSecondExon + 1, secondExon.end());
            }
        }
    }

    private static String bases(ChrBaseRegion chrBaseRegion, RefGenomeInterface refGenome)
    {
        return refGenome.getBaseString(chrBaseRegion.Chromosome, chrBaseRegion.start(), chrBaseRegion.end());
    }

    public int translateCodonPosition(final int position)
    {
        Preconditions.checkArgument(position >= 0);
        Preconditions.checkArgument(position <= 2);
        if (IsPositiveReadStrand)
        {
            if (RegionInSecondExon == null )
            {
                return CodonStart + position;
            }
            switch (position) {
                case 0: return RegionInFirstExon.start();
                case 1: return RegionInFirstExon.baseLength() == 2 ? RegionInFirstExon.end() : RegionInSecondExon.start();
                default: return RegionInSecondExon.end();
            }
        }
        else
        {
            if (RegionInSecondExon == null)
            {
                return CodonStart - position;
            }
            switch (position) {
                case 0: return RegionInFirstExon.end();
                case 1: return RegionInFirstExon.baseLength() == 2 ? RegionInFirstExon.start() : RegionInSecondExon.end();
                default: return RegionInSecondExon.start();
            }
        }
    }

    public boolean codonIsInSingleExon()
    {
        return RegionInSecondExon == null;
    }

    @NotNull
    public String retrieveCodon(final RefGenomeInterface refGenome)
    {
        return IsPositiveReadStrand ? codonForPositiveRead(refGenome) : codonForNegativeRead(refGenome);
    }

    @NotNull
    private String codonForPositiveRead(final RefGenomeInterface refGenome)
    {
        String part1 = bases(RegionInFirstExon, refGenome);
        String part2 = RegionInSecondExon == null ? "" : bases(RegionInSecondExon, refGenome);
        return part1 + part2;
    }

    @NotNull
    private String codonForNegativeRead(final RefGenomeInterface refGenome)
    {
        String part1 = Nucleotides.reverseComplementBases(bases(RegionInFirstExon, refGenome));
        String part2 = RegionInSecondExon == null ? "" : Nucleotides.reverseComplementBases(bases(RegionInSecondExon, refGenome));
        return part1 + part2;
    }
}
