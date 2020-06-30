package com.hartwig.hmftools.bachelor.types;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class VariantFilter
{
    private static final Logger LOGGER = LogManager.getLogger(VariantFilter.class);

    public final String Gene;
    public final String TranscriptId;
    public final String Chromosome;
    public final long Position;
    public final String Ref;
    public final String Alt;
    public final CodingEffect Effect;
    public final String HgvsProteinCodon;
    public final int MinCodon;
    public final String ClinvarSignificance;
    public final String ClinvarSigInfo;

    public VariantFilter(final String gene, final String transcriptId,
            final String chromosome, long position, final String ref, final String alt, final CodingEffect effect,
            final String hgvsProteinCodon, final String clinvarSig, final String clinvarSigInfo, int minCodon)
    {
        Gene = gene;
        TranscriptId = transcriptId;
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Effect = effect;
        HgvsProteinCodon = hgvsProteinCodon;
        MinCodon = minCodon;
        ClinvarSignificance = clinvarSig;
        ClinvarSigInfo = clinvarSigInfo;
    }

    public boolean blacklistMatch(String gene, String chromosome, long position, String ref, String alt, int proteinPosition)
    {
        if (MinCodon >= 0)
        {
            if (proteinPosition >=0 && MinCodon <= proteinPosition)
            {
                LOGGER.debug("Gene({}) var({}:{}) ref({}) alt({}) matches filter on minCodon({})",
                        gene, chromosome, position, ref, alt, proteinPosition);

                return true;
            }
        }

        if (Chromosome.equals(chromosome) && Position == position && Ref.equals(ref) && Alt.equals(alt))
        {
            LOGGER.debug("Gene({}) var({}:{}) ref({}) alt({}) matches filter on position, ref & alt",
                    gene, chromosome, position, ref, alt);

            return true;
        }

        return false;
    }

    public boolean whitelistMatch(String gene, String chromosome, long position, String ref, String alt,
            CodingEffect codingEffect, String hgvsProtein)
    {
        if(codingEffect == MISSENSE)
        {
            if (!HgvsProteinCodon.isEmpty() && HgvsProteinCodon.equals(hgvsProtein))
            {
                LOGGER.debug("Gene({}) matches filter on hgvsProtein({})", gene, hgvsProtein);

                return true;
            }
        }
        else
        {
            if (Chromosome.equals(chromosome) && Position == position && Ref.equals(ref) && Alt.equals(alt))
            {
                LOGGER.debug("Gene({}) var({}:{}) ref({}) alt({}) matches filter on position, ref & alt",
                        gene, chromosome, position, ref, alt);

                return true;
            }
        }

        return false;
    }
}
