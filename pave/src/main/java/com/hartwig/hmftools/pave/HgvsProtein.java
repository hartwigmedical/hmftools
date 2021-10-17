package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.START_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.isInframe;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.isNonsenseOrFrameshift;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_TYPE_DEL;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_TYPE_DUP;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_TYPE_INS;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_UNKNOWN;
import static com.hartwig.hmftools.pave.HgvsCoding.isDuplication;
import static com.hartwig.hmftools.pave.ProteinUtils.trimAminoAcids;

import java.util.List;
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

    public static String generate(final VariantData variant, final ProteinContext proteinContext, final List<VariantEffect> effects)
    {
        if(effects.size() == 1)
            return generate(variant, proteinContext, effects.get(0));

        boolean hasStopLost = false;
        boolean hasStopGained = false;
        VariantEffect topEffect = null;

        for(VariantEffect effect : effects)
        {
            if(effect == STOP_GAINED)
                hasStopGained = true;
            else if(effect == STOP_LOST)
                hasStopLost = true;
            else if(topEffect == null)
                topEffect = effect;
        }

        String hgvs = generate(variant, proteinContext, topEffect);

        if(hasStopGained)
            hgvs += "*";
        else if(hasStopLost)
            hgvs += "ext*?";

        return hgvs;
    }

    public static String generate(final VariantData variant, final ProteinContext proteinContext, final VariantEffect effect)
    {
        if(!proteinContext.validRefCodon() || !proteinContext.validAltCodon())
            return HGVS_UNKNOWN;

        trimAminoAcids(proteinContext);

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

            case START_LOST:
                formStartLost(proteinContext, sb);
                break;

            case INFRAME_DELETION:
            case PHASED_INFRAME_DELETION:
                formInframeDeletion(proteinContext, sb);
                break;

            case INFRAME_INSERTION:
            case PHASED_INFRAME_INSERTION:
                formInframeInsertion(variant, proteinContext, sb);
                break;

            case FRAMESHIFT:
                formFrameshift(variant, proteinContext, sb);
                break;

            /*
            case STOP_LOST:
                formStopLost(variant, proteinContext, sb);
                break;

            case STOP_GAINED:
                formStopGained(proteinContext, sb);
                break;
            */

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

    private static void formInframeDeletion(final ProteinContext proteinContext, final StringBuilder sb)
    {
        // conservative means only whole codons are deleted, eg p.Gly4_Gln6del
        // conservative single AA: p.Lys2del
        // conservative multi: p.Gly4_Gln6del
        // non-con: p.Cys28_Lys29delinsTrp
        // only show the straddling bases (ie first and last) - but check that the alt and ref aren't including the same one at start or end
        int aaIndexStart = proteinContext.NetCodonIndexRange[SE_START];
        int aaIndexEnd = proteinContext.NetCodonIndexRange[SE_END];
        String refAminoAcids = proteinContext.NetRefAminoAcids;
        String altAminoAcids = proteinContext.NetAltAminoAcids;

        if(refAminoAcids.length() == 1)
        {
            sb.append(convertToTriLetters(refAminoAcids.charAt(0)));
            sb.append(aaIndexStart);
            sb.append(HGVS_TYPE_DEL);
        }
        else
        {
            sb.append(convertToTriLetters(refAminoAcids.charAt(0)));
            sb.append(aaIndexStart);

            sb.append('_');
            sb.append(convertToTriLetters(refAminoAcids.charAt(refAminoAcids.length() - 1)));
            sb.append(aaIndexEnd);
            sb.append(HGVS_TYPE_DEL);

            if(!altAminoAcids.isEmpty())
            {
                sb.append(HGVS_TYPE_INS);
                sb.append(convertToTriLetters(altAminoAcids.charAt(0)));
            }
        }
    }

    private static void formInframeInsertion(final VariantData variant, final ProteinContext proteinContext, final StringBuilder sb)
    {
        int aaIndexStart = proteinContext.CodonIndex;
        int aaIndexEnd = aaIndexStart + 1;
        String refAminoAcids = proteinContext.RefAminoAcids;
        String altAminoAcids = proteinContext.NetAltAminoAcids;

        if(isDuplication(variant))
        {
            // duplication (single AA) p.Gln8dup
            // duplication (range) p.Gly4_Gln6dup
            String insertedBases = variant.Alt.substring(1);
            int insertLength = insertedBases.length();
            int dupCount = insertLength / 3;

            sb.append(convertToTriLetters(refAminoAcids.charAt(0)));
            sb.append(aaIndexStart);

            if(dupCount > 1)
            {
                sb.append('_');
                sb.append(convertToTriLetters(refAminoAcids.charAt(refAminoAcids.length() - 1)));
                sb.append(aaIndexStart + 1);
            }

            sb.append(HGVS_TYPE_DUP);
        }
        else
        {
            // insertion p.Lys2_Leu3insGlnSer
            // insertion (conservative stop) p.Ser81_Val82ins*
            // insertion (non conservative) p.Cys28delinsTrpVal
            sb.append(convertToTriLetters(refAminoAcids.charAt(0)));
            sb.append(aaIndexStart);

            if(proteinContext.RefAminoAcids.charAt(0) == proteinContext.AltAminoAcids.charAt(0))
            {
                // conservative
                sb.append('_');
                sb.append(convertToTriLetters(refAminoAcids.charAt(refAminoAcids.length() - 1)));
                sb.append(aaIndexStart + 1);
                sb.append(HGVS_TYPE_INS);
                sb.append(convertToTriLetters(altAminoAcids));
            }
            else
            {
                sb.append(HGVS_TYPE_DEL);
                sb.append(HGVS_TYPE_INS);
                sb.append(convertToTriLetters(altAminoAcids));
            }
        }
    }

    private static void formFrameshift(final VariantData variant, final ProteinContext proteinContext, final StringBuilder sb)
    {
        if(variant.isDeletion())
        {
            formInframeDeletion(proteinContext, sb);
        }
        else
        {
            formInframeInsertion(variant, proteinContext, sb);
        }

        sb.append("fs");
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
