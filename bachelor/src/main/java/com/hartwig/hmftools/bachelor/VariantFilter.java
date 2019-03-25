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

    public VariantFilter(final String gene, final String trancriptId,
            final String chromosome, long position, final String ref, final String alt,
            final CodingEffect effect, final String hgvsProteinCodon, final String dbSnpId, int minCodon)
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
    }

    /*
    public boolean matches(final BachelorGermlineVariant var)
    {
        if(!var.chromosome().equals(Chromosome) || var.position() != Position
                || !var.ref().equals(Ref) || !var.alts().equals(Alt))
        {
            return false;
        }

        for(final String effect : var.effectsList())
        {
            if (STOP_GAINED.isParentTypeOf(var.effects()) && Effects.contains("nonsense"))
                return true;

            if (Effects.contains(effect))
                return true;
        }

        return false;
    }
    */


}
