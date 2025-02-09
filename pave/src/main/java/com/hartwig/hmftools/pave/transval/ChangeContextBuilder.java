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
    ChangeContext build(String chromosome, RefGenomeInterface genome, List<ChrBaseRegion> codingRegions)
    {
        ChrBaseRegion exon = codingRegions.get(ExonIndex);
        String exonBases = genome.getBaseString(chromosome, exon.start(), exon.end());
        String prefix = genome.getBaseString(chromosome,  exon.start() - 5, exon.start() - 1);
        String basesInPreviousExon = "";
        if(PaddingInPreviousExon > 0)
        {
            ChrBaseRegion previousExon = codingRegions.get(ExonIndex - 1);
            int start = previousExon.end() - PaddingInPreviousExon + 1;
            int stop = previousExon.end();
            basesInPreviousExon = genome.getBaseString(chromosome, start, stop);
        }
        String basesInNextExon = "";
        if(PaddingInNextExon > 0)
        {
            ChrBaseRegion nextExon = codingRegions.get(ExonIndex + 1);
            int start = nextExon.start();
            int stop = start + PaddingInNextExon - 1 ;
            basesInNextExon = genome.getBaseString(chromosome, start, stop );
        }

        PaddedExon containingExon = new PaddedExon(basesInPreviousExon, basesInNextExon, exonBases, exon.start(), prefix);
        return new ChangeContext(containingExon, ChangeStart, ChangeEnd, true, AminoAcidNumberOfFirstAminoAcidStartingInExon);
    }
}
