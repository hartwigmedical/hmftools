package com.hartwig.hmftools.common.cider;

import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.jetbrains.annotations.Nullable;

public record IgTcrGene(
        String geneName,
        String allele, // 01 etc
        IgTcrRegion region,
        IgTcrFunctionality functionality,
        @Nullable ChrBaseRegion geneLocation,
        @Nullable Strand geneStrand,
        @Nullable String altAssemblyName, // not null if is not in primary assembly
        @Nullable String anchorSequence,  // only valid for V / J gene
        @Nullable ChrBaseRegion anchorLocation)
{
    public String geneAllele()
    {
        return geneName + "*" + allele;
    }

    public boolean isFunctional()
    {
        return functionality == IgTcrFunctionality.FUNCTIONAL;
    }

    public boolean inPrimaryAssembly()
    {
        return altAssemblyName == null;
    }
}
