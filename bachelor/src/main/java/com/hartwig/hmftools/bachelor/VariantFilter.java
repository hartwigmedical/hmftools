package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

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

    private static final Logger LOGGER = LogManager.getLogger(VariantFilter.class);

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

    public boolean blacklistMatch(String gene, String chromosome, long position, String ref, String alt, int proteinPosition)
    {
        /*
        if (!HgvsProteinCodon.isEmpty())
        {
            if (HgvsProteinCodon.equals(hgvsProtein))
            {
                LOGGER.debug("gene({}) var({}:{}) ref({}) alt({}) matches filter on hgvsProtein({})",
                        gene, chromosome, position, ref, alt, hgvsProtein);

                return true;
            }
        }

        if (!DBSnpId.isEmpty())
        {
            if (varId.contains(DBSnpId))
            {
                LOGGER.debug("gene({}) var({}:{}) ref({}) alt({}) matches filter on DBSnpId({})",
                        gene, chromosome, position, ref, alt, varId);

                return true;
            }
        }
        */

        if (MinCodon >= 0)
        {
            if (proteinPosition >=0 && MinCodon <= proteinPosition)
            {
                LOGGER.debug("gene({}) var({}:{}) ref({}) alt({}) matches filter on minCodon({})",
                        gene, chromosome, position, ref, alt, proteinPosition);

                return true;
            }
        }

        if (Chromosome.equals(chromosome) && Position == position && Ref.equals(ref) && Alt.equals(alt))
        {
            LOGGER.debug("gene({}) var({}:{}) ref({}) alt({}) matches filter on position, ref & alt",
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
                LOGGER.debug("gene({}) matches filter on hgvsProtein({})", gene, hgvsProtein);

                return true;
            }
        }
        else
        {
            if (Chromosome.equals(chromosome) && Position == position && Ref.equals(ref) && Alt.equals(alt))
            {
                LOGGER.debug("gene({}) var({}:{}) ref({}) alt({}) matches filter on position, ref & alt",
                        gene, chromosome, position, ref, alt);

                return true;
            }
        }

        return false;
    }


    }
