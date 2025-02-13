package com.hartwig.hmftools.pave.transval;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

public class ChangeContextBuilder
{
    final int ExonIndex;
    final int ChangeStart;
    final int ChangeEnd;
    final int PaddingInPreviousExon;
    final int PaddingInNextExon;
    final int AminoAcidNumberOfFirstAminoAcidStartingInExon;

    public ChangeContextBuilder(final int exonIndex,
            final int changeStart,
            final int changeEnd,
            final int paddingInPreviousExon,
            final int paddingInNextExon,
            final int numberOfCodonsStartingInPreviousExons)
    {
        ExonIndex = exonIndex;
        ChangeStart = changeStart;
        ChangeEnd = changeEnd;
        PaddingInPreviousExon = paddingInPreviousExon;
        PaddingInNextExon = paddingInNextExon;
        AminoAcidNumberOfFirstAminoAcidStartingInExon = numberOfCodonsStartingInPreviousExons + 1;
    }

    @NotNull
    ChangeContext build(String chromosome, RefGenomeInterface genome, List<ChrBaseRegion> codingRegions, boolean isPositiveStrand)
    {
        ChrBaseRegion exon = codingRegions.get(ExonIndex);
        String exonBases = genome.getBaseString(chromosome, exon.start(), exon.end());
        String prefix = genome.getBaseString(chromosome, exon.start() - 5, exon.start() - 1);
        String basesInPreviousExon = "";
        if(PaddingInPreviousExon > 0)
        {
            ChrBaseRegion previousExon = codingRegions.get(ExonIndex - 1);
            if(isPositiveStrand)
            {
                basesInPreviousExon = lastN(previousExon, PaddingInPreviousExon, chromosome, genome);
            }
            else
            {
                basesInPreviousExon = firstN(previousExon, PaddingInPreviousExon, chromosome, genome);
            }
        }
        String basesInNextExon = "";
        if(PaddingInNextExon > 0)
        {
            ChrBaseRegion nextExon = codingRegions.get(ExonIndex + 1);
            if(isPositiveStrand)
            {
                basesInNextExon = firstN(nextExon, PaddingInNextExon, chromosome, genome);
            }
            else
            {
                basesInNextExon = lastN(nextExon, PaddingInNextExon, chromosome, genome);
            }
        }

        PaddedExon containingExon;
        if(isPositiveStrand)
        {
            containingExon = new PaddedExon(basesInPreviousExon, basesInNextExon, exonBases, exon.start(), prefix);
        }
        else
        {
            containingExon = new PaddedExon(basesInNextExon, basesInPreviousExon, exonBases, exon.start(), prefix);
        }
        return new ChangeContext(containingExon, ChangeStart, ChangeEnd, isPositiveStrand, AminoAcidNumberOfFirstAminoAcidStartingInExon);
    }

    private String lastN(ChrBaseRegion exon, int n, String chromosome, RefGenomeInterface genome)
    {
        int end = exon.end();
        int start = end - n + 1;

        return genome.getBaseString(chromosome, start, end);
    }

    private String firstN(ChrBaseRegion exon, int n, String chromosome, RefGenomeInterface genome)
    {
        int start = exon.start();
        int end = start + n - 1;
        return genome.getBaseString(chromosome, start, end);
    }
}
