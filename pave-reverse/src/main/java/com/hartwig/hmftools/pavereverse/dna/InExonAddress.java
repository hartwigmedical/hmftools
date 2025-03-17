package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.pavereverse.gene.ExonForLocation;
import com.hartwig.hmftools.pavereverse.gene.GeneTranscript;

public class InExonAddress implements HgvsAddress
{
    public final int IndexOfBaseInCodingBases;

    public InExonAddress(final int indexOfBaseInCodingBases)
    {
        this.IndexOfBaseInCodingBases = indexOfBaseInCodingBases;
    }

    @Override
    public int toStrandLocation(GeneTranscript geneTranscript)
    {
        ExonForLocation exonForLocation = new ExonForLocation(IndexOfBaseInCodingBases - 1, geneTranscript.codingRegionLengths());
        ChrBaseRegion exonData = geneTranscript.CodingRegions.get(exonForLocation.ExonIndex);
        return exonData.start() + exonForLocation.LocationInExon; // start is location of base 1
    }
}
