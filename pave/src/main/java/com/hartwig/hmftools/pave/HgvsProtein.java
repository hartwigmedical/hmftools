package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.START_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.isInframe;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_TYPE_DEL;
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
        }
    }

    public static String generate(final VariantData variant, final ProteinContext proteinContext, final VariantEffect effect)
    {
        if(!proteinContext.validRefCodon())
            return HGVS_UNKNOWN;

        StringBuilder sb = new StringBuilder();
        sb.append(PROTEIN_ID);

        if(isInframe(effect))
        {
            formInframe(variant, proteinContext, sb);
        }
        else
        {
            String refAminoAcids = convertToTriLetters(proteinContext.RefAminoAcids);
            String altAminoAcids = !proteinContext.AltAminoAcids.isEmpty() ? convertToTriLetters(proteinContext.AltAminoAcids) : "";
            int aaIndex = proteinContext.CodonIndex;

            sb.append(refAminoAcids);
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
                sb.append(altAminoAcids);
                sb.append("ext*?");
            }
            else if(effect == FRAMESHIFT)
            {
                sb.append("fs");
            }
            else if(isInframe(effect))
            {
                formInframe(variant, proteinContext, sb);
            }
            else
            {
                sb.append(HGVS_UNKNOWN);
            }
        }

        return sb.toString();
    }

    private static void formInframe(
            final VariantData variant, final ProteinContext proteinContext, final StringBuilder sb)
    {
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

        sb.append(refTriAAs);
        sb.append(aaIndex);

        if(variant.isDeletion())
        {
            int delLength = refAminoAcids.length() - altAminoAcids.length();

            if(delLength == 1)
            {
                sb.append(HGVS_TYPE_DEL);
            }
            else
            {
                sb.append("_");
            }
        }
        else if(variant.isBaseChange())
        {

        }
        else if(isDuplication(variant))
        {

        }

        /*
        if(variant.isDeletion())
        {
            formDeletion(variant, proteinContext, sb);
        }
        else if(isDuplication(variant))
        {
            formDuplication(variant, proteinContext, sb);
        }
        else
        {
            formInsertion(variant, proteinContext, sb);
        }
        */
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


    private static String convertToTriLetters(final String aminoAcids)
    {
        StringBuilder sb = new StringBuilder();

        for(int i = 0; i < aminoAcids.length(); ++i)
        {
            sb.append(AMINO_ACID_SINGLE_TO_TRI_MAP.get(aminoAcids.charAt(i)));
        }

        return sb.toString();
    }
}
