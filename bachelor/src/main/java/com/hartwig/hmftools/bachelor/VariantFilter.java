package com.hartwig.hmftools.bachelor;

import com.hartwig.hmftools.common.variant.CodingEffect;

public class VariantFilter
{
    public final String Gene;
    public final String TranscriptId;
    public final String Chromosome;
    public final long Position;
    public final String Ref;
    public final String Alt;
    public final CodingEffect Effect;
    public final String HgvsProteinCodon;
    public final String DBSnpId;
    public final int MinCodon;
    public final String ClinvarSignificance;

    public VariantFilter(final String gene, final String trancriptId,
            final String chromosome, long position, final String ref, final String alt, final CodingEffect effect,
            final String hgvsProteinCodon, final String dbSnpId, final String clinvarSig, int minCodon)
    {
        Gene = gene;
        TranscriptId = trancriptId;
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Effect = effect;
        HgvsProteinCodon = hgvsProteinCodon;
        DBSnpId = dbSnpId;
        MinCodon = minCodon;
        ClinvarSignificance = clinvarSig;
    }

}
