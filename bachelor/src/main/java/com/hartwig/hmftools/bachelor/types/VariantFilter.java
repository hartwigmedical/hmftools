package com.hartwig.hmftools.bachelor.types;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BACH_LOGGER;
import static com.hartwig.hmftools.bachelor.types.PathogenicType.CLINVAR_BENIGN;
import static com.hartwig.hmftools.bachelor.types.PathogenicType.CLINVAR_CONFLICTING;
import static com.hartwig.hmftools.bachelor.types.PathogenicType.CLINVAR_LIKELY_BENIGN;
import static com.hartwig.hmftools.bachelor.types.PathogenicType.CLINVAR_LIKELY_PATHOGENIC;
import static com.hartwig.hmftools.bachelor.types.PathogenicType.CLINVAR_PATHOGENIC;
import static com.hartwig.hmftools.bachelor.types.PathogenicType.UNANNOTATED;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;

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
    public final int MinCodon;
    public final boolean Configured;
    public final String ClinvarSignificance;
    public final String ClinvarSigInfo;

    public VariantFilter(final String gene, final String transcriptId,
            final String chromosome, long position, final String ref, final String alt, final CodingEffect effect,
            final String hgvsProteinCodon, final String clinvarSig, final String clinvarSigInfo, int minCodon, boolean configured)
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
        Configured = configured;
    }

    public boolean blacklistMatch(String gene, String chromosome, long position, String ref, String alt, int proteinPosition)
    {
        return blacklistMinCodonMatch(gene, chromosome, position, ref, alt, proteinPosition)
                || blacklistSpecificVariantMatch(gene, chromosome, position, ref, alt);
    }

    public boolean blacklistMinCodonMatch(String gene, String chromosome, long position, String ref, String alt, int proteinPosition)
    {
        if (MinCodon >= 0 && proteinPosition >=0 && MinCodon <= proteinPosition)
        {
            BACH_LOGGER.debug("Gene({}) var({}:{}) ref({}) alt({}) matches filter on minCodon({})",
                    gene, chromosome, position, ref, alt, proteinPosition);

            return true;
        }

        return false;
    }

    public boolean blacklistSpecificVariantMatch(String gene, String chromosome, long position, String ref, String alt)
    {
        if (Chromosome.equals(chromosome) && Position == position && Ref.equals(ref) && Alt.equals(alt))
        {
            BACH_LOGGER.debug("Gene({}) var({}:{}) ref({}) alt({}) matches filter on position, ref & alt",
                    gene, chromosome, position, ref, alt);

            return true;
        }

        return false;
    }

    public boolean whitelistMatch(String gene, String chromosome, long position, String ref, String alt,
            CodingEffect codingEffect, String hgvsProtein)
    {
        if (!HgvsProteinCodon.isEmpty() && HgvsProteinCodon.equals(hgvsProtein))
        {
            BACH_LOGGER.debug("Gene({}) matches filter on hgvsProtein({})", gene, hgvsProtein);
            return true;
        }
        else if (Chromosome.equals(chromosome) && Position == position && Ref.equals(ref) && Alt.equals(alt))
        {
            BACH_LOGGER.debug("Gene({}) var({}:{}) ref({}) alt({}) matches filter on position, ref & alt",
                    gene, chromosome, position, ref, alt);
            return true;
        }

        return false;
    }

    private static final String CLINVAR_STR_PATHOGENIC = "Pathogenic";
    private static final String CLINVAR_STR_LIKELY_PATHOGENIC = "Likely_pathogenic";
    private static final String CLINVAR_STR_BENIGN = "Benign";
    private static final String CLINVAR_STR_LIKELY_BENIGN = "Likely_benign";
    private static final String CLINVAR_STR_CONFLICTING = "Conflicting";

    public static boolean hasPathogenic(final String clinvarStr) { return clinvarStr.contains(CLINVAR_STR_PATHOGENIC); }
    public static boolean hasLikelyPathogenic(final String clinvarStr) { return clinvarStr.contains(CLINVAR_STR_LIKELY_PATHOGENIC); }
    public static boolean hasBenign(final String clinvarStr) { return clinvarStr.contains(CLINVAR_STR_BENIGN); }
    public static boolean hasLikelyBenign(final String clinvarStr) { return clinvarStr.contains(CLINVAR_STR_LIKELY_BENIGN); }
    public static boolean isConflicting(final String clinvarStr)
    {
        return clinvarStr.contains(CLINVAR_STR_CONFLICTING);
    }

    public PathogenicType determinePathogenicType()
    {
        if(ClinvarSignificance.isEmpty())
            return UNANNOTATED;

        /*
         * CLINVAR_PATHOGENIC - At least 1 interpretation of 'PATHOGENIC' and none ‘BENIGN’ or ‘LIKELY_BENIGN’
         * CLINVAR_LIKELY_PATHOGENIC - No interpretation of PATHOGENIC, but at least 1 interpretation of 'LIKELY_PATHOGENIC' and none ‘BENIGN’ or ‘LIKELY_BENIGN’
         * CLINVAR_CONFLICTING - Variant has both likely 'BENIGN'/'LIKELY_BENIGN' and 'PATHOGENIC'/'LIKELY_PATHOGENIC' interpretations
         * CLINVAR_LIKELY_BENIGN - No interpretation of 'BENIGN' and at least 1 interpretation of 'LIKELY_BENIGN' and none ‘PATHOGENIC’ or ‘LIKELY_PATHOGENIC’
         * CLINVAR_BENIGN - At least 1 interpretation of 'BENIGN' and none ‘PATHOGENIC’ or ‘LIKELY_PATHOGENIC’
         */

        if(hasPathogenic(ClinvarSignificance))
            return CLINVAR_PATHOGENIC;

        if(hasLikelyPathogenic(ClinvarSignificance))
            return CLINVAR_LIKELY_PATHOGENIC;

        if(hasBenign(ClinvarSignificance))
            return CLINVAR_BENIGN;

        if(hasLikelyBenign(ClinvarSignificance))
            return CLINVAR_LIKELY_BENIGN;

        if(isConflicting(ClinvarSignificance))
        {
            if(hasPathogenic(ClinvarSigInfo) && !hasBenign(ClinvarSigInfo) && !hasLikelyBenign(ClinvarSigInfo))
                return CLINVAR_PATHOGENIC;

            if(hasLikelyPathogenic(ClinvarSigInfo) && !hasPathogenic(ClinvarSigInfo) && !hasBenign(ClinvarSigInfo) && !hasLikelyBenign(ClinvarSigInfo))
                return CLINVAR_LIKELY_PATHOGENIC;

            if(hasBenign(ClinvarSigInfo) && !hasPathogenic(ClinvarSigInfo) && !hasLikelyPathogenic(ClinvarSigInfo))
                return CLINVAR_BENIGN;

            if(hasLikelyBenign(ClinvarSigInfo) && !hasPathogenic(ClinvarSigInfo) && !hasLikelyPathogenic(ClinvarSigInfo) && !hasBenign(ClinvarSigInfo))
                return CLINVAR_LIKELY_BENIGN;

            if((hasPathogenic(ClinvarSigInfo) || hasLikelyPathogenic(ClinvarSigInfo)) && (hasBenign(ClinvarSigInfo) || hasLikelyBenign(ClinvarSigInfo)))
            {
                return CLINVAR_CONFLICTING;
            }
        }

        return UNANNOTATED;
    }
}
