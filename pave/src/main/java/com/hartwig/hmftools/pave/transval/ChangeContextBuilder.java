package com.hartwig.hmftools.pave.transval;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.NotNull;

class ChangeContextData
{
    final int ExonIndex;
    final int ChangeStart;
    final int ChangeEnd;
    final int PaddingInPreviousExon;
    final int PaddingInNextExon;
    final int AminoAcidNumberOfFirstAminoAcidStartingInExon;

    ChangeContextData(final int exonIndex,
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
}

class ChangeContextBuilder
{
    @NotNull
    final ChangeContextData Data;
    final ChangeContextData CompanionData;

    public ChangeContextBuilder(@NotNull final ChangeContextData data, @NotNull final ChangeContextData companionData)
    {
        this.Data = data;
        this.CompanionData = companionData;
    }

    public ChangeContextBuilder(final int exonIndex,
            final int changeStart,
            final int changeEnd,
            final int paddingInPreviousExon,
            final int paddingInNextExon,
            final int numberOfCodonsStartingInPreviousExons)
    {
        this.Data = new ChangeContextData(exonIndex, changeStart, changeEnd, paddingInPreviousExon, paddingInNextExon, numberOfCodonsStartingInPreviousExons);
        this.CompanionData = null;
    }

    @NotNull
    ChangeContext build(String chromosome, RefGenomeInterface genome, List<ChrBaseRegion> codingRegions, boolean isPositiveStrand)
    {
        ChangeContext result = build(Data, chromosome, genome, codingRegions, isPositiveStrand);
        if (CompanionData != null)
        {
            result.setCompanionContext(build(CompanionData, chromosome, genome, codingRegions, isPositiveStrand));
        }
        return result;
    }

    private ChangeContext build(ChangeContextData ccd, String chromosome,RefGenomeInterface genome, List<ChrBaseRegion> codingRegions, boolean isPositiveStrand)
    {
        ChrBaseRegion exon = codingRegions.get(ccd.ExonIndex);
        String exonBases = genome.getBaseString(chromosome, exon.start(), exon.end());
        String prefix = genome.getBaseString(chromosome, exon.start() - 5, exon.start() - 1);
        String suffix = genome.getBaseString(chromosome, exon.end() + 1, exon.end() + 5);
        String basesInPreviousExon = "";
        if(ccd.PaddingInPreviousExon > 0)
        {
            ChrBaseRegion previousExon = codingRegions.get(ccd.ExonIndex - 1);
            if(isPositiveStrand)
            {
                basesInPreviousExon = lastN(previousExon, ccd.PaddingInPreviousExon, chromosome, genome);
            }
            else
            {
                basesInPreviousExon = firstN(previousExon, ccd.PaddingInPreviousExon, chromosome, genome);
            }
        }
        String basesInNextExon = "";
        if(ccd.PaddingInNextExon > 0)
        {
            ChrBaseRegion nextExon = codingRegions.get(ccd.ExonIndex + 1);
            if(isPositiveStrand)
            {
                basesInNextExon = firstN(nextExon, ccd.PaddingInNextExon, chromosome, genome);
            }
            else
            {
                basesInNextExon = lastN(nextExon, ccd.PaddingInNextExon, chromosome, genome);
            }
        }

        PaddedExon containingExon;
        if(isPositiveStrand)
        {
            containingExon = new PaddedExon(ccd.ExonIndex, basesInPreviousExon, basesInNextExon, exonBases, exon.start(), prefix, suffix);
        }
        else
        {
            containingExon = new PaddedExon(ccd.ExonIndex, basesInNextExon, basesInPreviousExon, exonBases, exon.start(), prefix, suffix);
        }
        return new ChangeContext(containingExon, ccd.ChangeStart, ccd.ChangeEnd, isPositiveStrand, ccd.AminoAcidNumberOfFirstAminoAcidStartingInExon);
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
