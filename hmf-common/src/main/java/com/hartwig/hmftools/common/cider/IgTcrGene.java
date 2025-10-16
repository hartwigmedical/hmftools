package com.hartwig.hmftools.common.cider;

import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.region.BaseRegion;

import org.jetbrains.annotations.Nullable;

public record IgTcrGene(
        String geneName,
        String allele, // 01 etc
        IgTcrRegion region,
        IgTcrFunctionality functionality,
        String sequence,
        @Nullable String contigName,    // From the reference genome. Can be a chromosome or another contig type.
        @Nullable BaseRegion genePosition,
        @Nullable Strand geneStrand,
        @Nullable String anchorSequence,  // only valid for V / J gene
        @Nullable BaseRegion anchorPosition)
{
    public String getGeneAllele()
    {
        return geneName + "*" + allele;
    }

    public boolean isFunctional()
    {
        return functionality == IgTcrFunctionality.FUNCTIONAL;
    }
}
