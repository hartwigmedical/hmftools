package com.hartwig.hmftools.pavereverse.base;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class ChangeContextBuilder
{
    public final ChangeContextData Data;
    public final ChangeContextData CompanionData;

    public ChangeContextBuilder(ChangeContextData data, ChangeContextData companionData)
    {
        Data = data;
        CompanionData = companionData;
    }

    public ChangeContextBuilder(int exonIndex,
            int changeStart,
            int changeEnd,
            int paddingInPreviousExon,
            int paddingInNextExon,
            int numberOfCodonsStartingInPreviousExons)
    {
        Data =
                new ChangeContextData(exonIndex, changeStart, changeEnd, paddingInPreviousExon, paddingInNextExon, numberOfCodonsStartingInPreviousExons);
        CompanionData = null;
    }

    public ChangeContext build(String chromosome, RefGenomeInterface genome, List<ChrBaseRegion> codingRegions, boolean isPositiveStrand)
    {
        ChangeContext result = build(Data, chromosome, genome, codingRegions, isPositiveStrand);
        if(CompanionData != null)
        {
            result.setCompanionContext(build(CompanionData, chromosome, genome, codingRegions, isPositiveStrand));
        }
        return result;
    }

    private ChangeContext build(ChangeContextData ccd, String chromosome, RefGenomeInterface genome, List<ChrBaseRegion> codingRegions,
            boolean isPositiveStrand)
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
            // If a transcript is incomplete, the sum of its coding region lengths
            // might not be a multiple of 3, in which case the PaddingInNextExon
            // value is non-zero, despite there not being a next exon.
            if(ccd.ExonIndex == codingRegions.size() - 1)
            {
                throw new IllegalArgumentException("Transcript is incomplete.");
            }
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
