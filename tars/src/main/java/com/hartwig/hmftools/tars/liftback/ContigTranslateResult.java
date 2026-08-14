package com.hartwig.hmftools.tars.liftback;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// genome-space result of translating a transcript-contig alignment
public record ContigTranslateResult(
        String chromosome,
        int genomicStart,
        Cigar genomicCigar)
{
    // the genomic CIGAR's N gaps, i.e. the re-inserted introns
    public List<BaseRegion> impliedIntrons()
    {
        List<BaseRegion> introns = new ArrayList<>();
        int pos = genomicStart;
        for(CigarElement element : genomicCigar.getCigarElements())
        {
            if(element.getOperator() == CigarOperator.N)
                introns.add(new BaseRegion(pos, pos + element.getLength() - 1));

            if(element.getOperator().consumesReferenceBases())
                pos += element.getLength();
        }
        return introns;
    }
}
