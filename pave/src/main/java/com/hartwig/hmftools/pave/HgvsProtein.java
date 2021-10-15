package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.isInframe;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.isNonsenseOrFrameshift;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_TYPE_DEL;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_TYPE_INS;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_UNKNOWN;
import static com.hartwig.hmftools.pave.HgvsCoding.isDuplication;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

public final class HgvsProtein
{
    /* Rules and conventions
     */

    private static final String PROTEIN_ID = "p.";

    private static final Map<Character,String> AMINO_ACID_SINGLE_TO_TRI_MAP = Maps.newHashMap();

    static
    {
        for(Map.Entry<String,String> entry : AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.entrySet())
        {
            AMINO_ACID_SINGLE_TO_TRI_MAP.put(entry.getValue().charAt(0), entry.getKey());
            AMINO_ACID_SINGLE_TO_TRI_MAP.put(STOP_AMINO_ACID, "Ter");
        }
    }

    public static boolean reportProteinImpact(final VariantEffect effect)
    {
        return isInframe(effect) || isNonsenseOrFrameshift(effect) || effect == SYNONYMOUS || effect == MISSENSE;
    }

    public static String generate(final VariantData variant, final ProteinContext proteinContext, final VariantEffect effect)
    {
        if(!proteinContext.validRefCodon() || !proteinContext.validAltCodon())
            return HGVS_UNKNOWN;

        StringBuilder sb = new StringBuilder();
        sb.append(PROTEIN_ID);

        switch(effect)
        {
            case SYNONYMOUS:
                formSynonymous(proteinContext, sb);
                break;

            case MISSENSE:
                formMissense(proteinContext, sb);
                break;

            case STOP_LOST:
                formStopLost(variant, proteinContext, sb);
                break;

            case STOP_GAINED:
                formStopGained(proteinContext, sb);
                break;

            case START_LOST:
                formStartLost(proteinContext, sb);
                break;

            case INFRAME_DELETION:
            case PHASED_INFRAME_DELETION:
                formInframeDeletion(variant, proteinContext, sb);
                break;

            case INFRAME_INSERTION:
            case PHASED_INFRAME_INSERTION:
                formInframeInsertion(variant, proteinContext, sb);
                break;

            case FRAMESHIFT:
                formFrameshift(variant, proteinContext, sb);
                break;

            default:
                sb.append(HGVS_UNKNOWN);
                break;
        }

        return sb.toString();
    }

    private static void formSynonymous(final ProteinContext proteinContext, final StringBuilder sb)
    {
        // ref and alt codons will match
        // SnpEff uses index only, but instead use form:  p.Leu54= for 1 codon or p.Leu54_Arg55= for MNV
        sb.append(convertToTriLetters(proteinContext.RefAminoAcids.charAt(0)));
        sb.append(proteinContext.CodonIndex);

        if(proteinContext.RefAminoAcids.length() > 1)
        {
            sb.append('_');
            sb.append(convertToTriLetters(proteinContext.RefAminoAcids.charAt(1)));
            sb.append(proteinContext.CodonIndex + 1);
        }

        sb.append('=');
    }

    private static void formMissense(final ProteinContext proteinContext, final StringBuilder sb)
    {
        // from SNV or MNV, 1 codon p.Arg64Lys or 2 p.Ala100_Val101delinsArgTrp
        // expect to have to matching number of ref and alt codons, typically 1 or 2, possible 3 for a longer MNV (none in prod)
        // current codon index will be the most upstream position since based on most upstream codon (open or already complete)
        if(proteinContext.RefAminoAcids.length() == 1)
        {
            sb.append(convertToTriLetters(proteinContext.RefAminoAcids.charAt(0)));
            sb.append(proteinContext.CodonIndex);
            sb.append(convertToTriLetters(proteinContext.AltAminoAcids.charAt(0)));
        }
        else
        {
            sb.append(convertToTriLetters(proteinContext.RefAminoAcids.charAt(0)));
            sb.append(proteinContext.CodonIndex);

            sb.append('_');
            sb.append(convertToTriLetters(proteinContext.RefAminoAcids.charAt(1)));
            sb.append(proteinContext.CodonIndex + 1);
            sb.append(HGVS_TYPE_DEL);
            sb.append(HGVS_TYPE_INS);
            sb.append(convertToTriLetters(proteinContext.AltAminoAcids.charAt(0)));
            sb.append(convertToTriLetters(proteinContext.AltAminoAcids.charAt(1)));
        }
    }

    private static void formStopGained(final ProteinContext proteinContext, final StringBuilder sb)
    {
        // might just be better to add a '*' onto the end
        // SNV p.Cys495*
        // MNV a) first codon becomes stop, as per SNV b) first changes and second becomes stop: p.Cys495_Val496delinsArg*
        // any INS: a) insert of stop codon by itself
        // any INS: b) insert of other new AAs then stop codon
        // any INS: c) change of first codon, insert of other new AAs then stop codon
        // any DEL: c) change of first codon, insert of other new AAs then stop codon
        //
        // single codon
        // p.Cys495_Val496delins*
        // p.CysVal495*, looks like the changed (ie Ref) AAs then *
        // may use
        sb.append("*");
    }

    private static void formStopLost(final VariantData variant, final ProteinContext proteinContext, final StringBuilder sb)
    {
        // - stop_lost - p.Ter494ext*? - don't show any new AAs after the stop since considered irrelevant
        sb.append(convertToTriLetters(proteinContext.RefAminoAcids.charAt(0)));
        sb.append(proteinContext.CodonIndex);
        sb.append(convertToTriLetters(proteinContext.AltAminoAcids.charAt(0)));
        sb.append("ext*?");
    }

    private static void formStartLost(final ProteinContext proteinContext, final StringBuilder sb)
    {
        // no other details are required
        sb.append("p.Met1?");
    }

    private static void formInframeDeletion(final VariantData variant, final ProteinContext proteinContext, final StringBuilder sb)
    {
    }

    private static void formInframeInsertion(final VariantData variant, final ProteinContext proteinContext, final StringBuilder sb)
    {
    }

    private static void formFrameshift(final VariantData variant, final ProteinContext proteinContext, final StringBuilder sb)
    {
        sb.append("fs");
    }

/*
    int aaIndex = proteinContext.CodonIndex;
        String refAminoAcids = proteinContext.RefAminoAcids;
        String altAminoAcids = proteinContext.AltAminoAcids;

        // strip out any matching AA from the start if the ref has more than 1
        if(!refAminoAcids.isEmpty() && !altAminoAcids.isEmpty() && refAminoAcids.charAt(0) == altAminoAcids.charAt(0))
        {
            refAminoAcids = refAminoAcids.substring(1);
            altAminoAcids = altAminoAcids.substring(1);
            aaIndex++;
        }

        String refTriAAs = convertToTriLetters(refAminoAcids);
        String altTriAAs = !altAminoAcids.isEmpty() ? convertToTriLetters(altAminoAcids) : "";

        if(isInframe(effect))
        {
            formInframe(variant, proteinContext, sb);
        }
        else
        {
            sb.append(refTriAAs);
            sb.append(aaIndex);

            if(effect == SYNONYMOUS)
            {
                sb.append("=");

                // no Stop retained
            }
            else if(effect == MISSENSE)
            {
                sb.append(altAminoAcids);
            }
            else if(effect == STOP_GAINED)
            {
                sb.append("*");
            }
            else if(effect == START_LOST)
            {
                sb.append("?");
            }
            else if(effect == STOP_LOST)
            {
                sb.append(altTriAAs);
                sb.append("ext*?");
            }
            else if(effect == FRAMESHIFT)
            {
                sb.append("fs");
            }
            else
            {
                sb.append(HGVS_UNKNOWN);
            }
        }

        return sb.toString();
    }

    private static void formInframe(
            final VariantData variant, final ProteinContext proteinContext,
            final String refTriAAs, final String altTriAAs, int aaIndex, final StringBuilder sb)
    {
        if(variant.isDeletion())
        {
            sb.append(refTriAAs);
            sb.append(aaIndex);
            sb.append(altTriAAs);

            int delLength = refTriAAs.length() - altTriAAs.length();

            if(delLength == 1)
            {
                sb.append(HGVS_TYPE_DEL);
            }
            else
            {
                sb.append("_");
            }
        }
        else
        {
            // p.Glu290_Glu291insAla
            // p.Leu631_Ala632insGlyLeu
            // p.Ser941_Asn942insThr
            // if(isDuplication(variant))

            if()

            sb.append(refTriAAs);
            sb.append(aaIndex);
            sb.append(altTriAAs);


        }
    }

    private static void formInsertion(final VariantData variant, final ProteinContext proteinContext, final StringBuilder sb)
    {
        // Insertion: p.Lys2_Leu3insGlnSer
        // Insertion (non conservative): p.Cys28delinsTrpVal

    }

    private static void formDuplication(final VariantData variant, final ProteinContext proteinContext, final StringBuilder sb)
    {
        // Duplication (single AA): p.Gln8dup
        // Duplication (range): p.Gly4_Gln6dup

    }

    private static void formDeletion(
            final VariantData variant, final ProteinContext proteinContext, final StringBuilder sb)
    {
        // shows first and last AA that are preserved ??
        // Deletion (single AA): p.Lys2del
        // Deletion (range): p.Gly4_Gln6del
        // Deletion (non conservative): p.Cys28_Lys29delinsTrp
        if(variant.isDeletion())
        {
            int delLength = proteinContext.RefAminoAcids.length() - proteinContext.AltAminoAcids.length();

            if(delLength == 1)
            {
                sb.append(HGVS_TYPE_DEL);
            }
            else
            {
                sb.append("_");
            }

            return;
        }

    }
    */

    private static String convertToTriLetters(final char aminoAcid)
    {
        String aaTriCode = AMINO_ACID_SINGLE_TO_TRI_MAP.get(aminoAcid);
        return aaTriCode != null ? aaTriCode : HGVS_UNKNOWN;
    }

    private static String convertToTriLetters(final String aminoAcids)
    {
        StringBuilder sb = new StringBuilder();

        for(int i = 0; i < aminoAcids.length(); ++i)
        {
            sb.append(convertToTriLetters(aminoAcids.charAt(i)));
        }

        return sb.toString();
    }
}
