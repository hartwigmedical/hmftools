package com.hartwig.hmftools.pave;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.MISSENSE;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.START_LOST;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.SYNONYMOUS;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.isInframe;
import static com.hartwig.hmftools.common.variant.impact.VariantEffect.isNonsenseOrFrameshift;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_TYPE_DEL;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_TYPE_DUP;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_TYPE_INS;
import static com.hartwig.hmftools.pave.HgvsCoding.HGVS_UNKNOWN;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

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
    private static final String HGVS_FRAMESHIFT = "fs";
    private static final String HGVS_STOP_LOST = "ext*?";
    private static final String HGVS_STOP_GAINED = "*";
    private static final String HGVS_STOP_TRI_CODE = "Ter";
    private static final String HGVS_SYNONYMOUS = "=";

    private static final Map<Character,String> AMINO_ACID_SINGLE_TO_TRI_MAP = Maps.newHashMap();

    static
    {
        for(Map.Entry<String,String> entry : AminoAcids.TRI_LETTER_AMINO_ACID_TO_SINGLE_LETTER.entrySet())
        {
            AMINO_ACID_SINGLE_TO_TRI_MAP.put(entry.getValue().charAt(0), entry.getKey());
            AMINO_ACID_SINGLE_TO_TRI_MAP.put(STOP_AMINO_ACID, HGVS_STOP_TRI_CODE);
        }
    }

    public static boolean reportProteinImpact(final VariantEffect effect)
    {
        return isInframe(effect) || isNonsenseOrFrameshift(effect) || effect == SYNONYMOUS || effect == MISSENSE;
    }

    public static String generate(final VariantData variant, final ProteinContext proteinContext, final List<VariantEffect> effects)
    {
        boolean hasStopGained = false;
        VariantEffect topEffect = null;

        for(VariantEffect effect : effects)
        {
            if(effect == STOP_GAINED)
                hasStopGained = true;
            else if(topEffect == null)
                topEffect = effect;
        }

        if(topEffect == null && hasStopGained) // re-instate for the protein string if effect was removed
            topEffect = MISSENSE;

        String hgvs = generate(proteinContext, topEffect);

        if(hasStopGained)
        {
            if(hgvs.contains(HGVS_STOP_TRI_CODE))
            {
                // replace from Ter onwards
                int stopIndex = hgvs.indexOf(HGVS_STOP_TRI_CODE);
                hgvs = hgvs.substring(0, stopIndex) + HGVS_STOP_GAINED;
            }
            else if(topEffect != FRAMESHIFT)
            {
                hgvs += HGVS_STOP_GAINED;
            }
        }

        return hgvs;
    }

    public static String generate(final ProteinContext proteinContext, final VariantEffect effect)
    {
        if(!proteinContext.validRefCodon())
            return HGVS_UNKNOWN;

        if(!proteinContext.validAltCodon() && effect != FRAMESHIFT) // all the others should have been able to get complete codons
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

            case START_LOST:
                formStartLost(proteinContext, sb);
                break;

            case STOP_LOST:
                formStopLost(proteinContext, sb);
                break;

            case INFRAME_DELETION:
            case PHASED_INFRAME_DELETION:
                formInframeDeletion(proteinContext, sb);
                break;

            case INFRAME_INSERTION:
            case PHASED_INFRAME_INSERTION:
                formInframeInsertion(proteinContext, sb);
                break;

            case FRAMESHIFT:
                formFrameshift(proteinContext, sb);
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

        sb.append(HGVS_SYNONYMOUS);
    }

    private static void formMissense(final ProteinContext proteinContext, final StringBuilder sb)
    {
        // from SNV or MNV, 1 codon p.Arg64Lys or 2 p.Ala100_Val101delinsArgTrp
        // expect to have to matching number of ref and alt codons, typically 1 or 2, possible 3 for a longer MNV (none in prod)
        // current codon index will be the most upstream position since based on most upstream codon (open or already complete)
        int aaIndexStart = proteinContext.NetCodonIndexRange[SE_START];
        int aaIndexEnd = proteinContext.NetCodonIndexRange[SE_END];
        String refAminoAcids = proteinContext.NetRefAminoAcids;
        String altAminoAcids = proteinContext.NetAltAminoAcids;

        if(proteinContext.NetRefAminoAcids.length() == 1)
        {
            sb.append(convertToTriLetters(refAminoAcids.charAt(0)));
            sb.append(aaIndexStart);
            sb.append(convertToTriLetters(altAminoAcids.charAt(0)));
        }
        else
        {
            sb.append(convertToTriLetters(refAminoAcids.charAt(0)));
            sb.append(aaIndexStart);

            sb.append('_');
            sb.append(convertToTriLetters(refAminoAcids.charAt(1)));
            sb.append(aaIndexEnd);
            sb.append(HGVS_TYPE_DEL);
            sb.append(HGVS_TYPE_INS);
            sb.append(convertToTriLetters(altAminoAcids.charAt(0)));
            sb.append(convertToTriLetters(altAminoAcids.charAt(1)));
        }
    }

    private static void formInframeDeletion(final ProteinContext proteinContext, final StringBuilder sb)
    {
        if(proteinContext.IsPhased && proteinContext.NetRefAminoAcids.length() > 1 && proteinContext.NetAltAminoAcids.length() > 1)
        {
            formPhasedInframe(proteinContext, sb);
            return;
        }

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

    private static void formInframeInsertion(final ProteinContext proteinContext, final StringBuilder sb)
    {
        if(proteinContext.IsPhased && proteinContext.NetRefAminoAcids.length() > 1 && proteinContext.NetAltAminoAcids.length() > 1)
        {
            formPhasedInframe(proteinContext, sb);
            return;
        }

        if(proteinContext.IsDuplication)
        {
            // duplication (single AA) p.Gln8dup
            // duplication (range) p.Gly4_Gln6dup - variety of AAs
            // dup range, single AA repeated: c.1014_1019dupTGCTGC	p.Ala339_Ala340dup
            String altAminoAcids = proteinContext.NetAltAminoAcids;
            int aaIndexStart = proteinContext.NetCodonIndexRange[SE_START];
            int aaIndexEnd = proteinContext.NetCodonIndexRange[SE_END];

            if(aaIndexEnd > aaIndexStart)
            {
                sb.append(convertToTriLetters(altAminoAcids.charAt(0)));
                sb.append(aaIndexStart);
                sb.append('_');
                sb.append(convertToTriLetters(altAminoAcids.charAt(altAminoAcids.length() - 1)));
                sb.append(proteinContext.NetCodonIndexRange[SE_END]);
            }
            else
            {
                sb.append(convertToTriLetters(altAminoAcids.charAt(0)));
                sb.append(aaIndexStart);
            }

            sb.append(HGVS_TYPE_DUP);
        }
        else
        {
            // insertion conservative single: AA AC -> ARC = p.A2_C3insR from net: '' -> R
            // insertion conservative multiple: KL -> KGSL = p.K2_L3insGS from net: '' -> GS
            // insertion conservative multiple: DN -> DKDN = p.D295_N296insKD from net: '' -> DN
            // insertion non-conservative: DCN -> DWVN = p.C28delinsWV from net: C -> WV
            // insertion conservative stop: p.Ser81_Val82ins*
            String refAminoAcids = proteinContext.RefAminoAcids;

            if(proteinContext.NetRefAminoAcids.isEmpty())
            {
                // skip past any redundant ref codons upstream
                int i = 0;

                for(; i < min(refAminoAcids.length() - 1, proteinContext.AltAminoAcids.length()); ++i)
                {
                    if(refAminoAcids.charAt(i) != proteinContext.AltAminoAcids.charAt(i))
                        break;
                }

                --i;

                if(i < 0 || i >= refAminoAcids.length() - 1)
                {
                    sb.append(HGVS_UNKNOWN);
                    return;
                }

                int aaStartIndex = proteinContext.CodonIndex + i;
                int aaEndIndex = aaStartIndex + 1;

                // conservative: only an insert, no deleted AAs so quote the range
                sb.append(convertToTriLetters(refAminoAcids.charAt(i)));
                sb.append(aaStartIndex);
                sb.append('_');
                sb.append(convertToTriLetters(refAminoAcids.charAt(i + 1)));
                sb.append(aaEndIndex);
            }
            else
            {
                // non-conservative: an AA has been deleted
                sb.append(convertToTriLetters(proteinContext.NetRefAminoAcids.charAt(0)));
                sb.append(proteinContext.NetCodonIndexRange[SE_START]);
                sb.append(HGVS_TYPE_DEL);
            }

            sb.append(HGVS_TYPE_INS);
            sb.append(convertToTriLetters(proteinContext.NetAltAminoAcids));
        }
    }

    private static void formPhasedInframe(final ProteinContext proteinContext, final StringBuilder sb)
    {
        // multiple AAs deleted and changed/inserted
        int aaIndexStart = proteinContext.NetCodonIndexRange[SE_START];
        int aaIndexEnd = proteinContext.NetCodonIndexRange[SE_END];
        String refAminoAcids = proteinContext.NetRefAminoAcids;
        String altAminoAcids = proteinContext.NetAltAminoAcids;

        sb.append(convertToTriLetters(refAminoAcids.charAt(0)));
        sb.append(aaIndexStart);
        sb.append('_');
        sb.append(convertToTriLetters(refAminoAcids.charAt(refAminoAcids.length() - 1)));
        sb.append(aaIndexEnd);
        sb.append(HGVS_TYPE_DEL);
        sb.append(HGVS_TYPE_INS);
        sb.append(convertToTriLetters(altAminoAcids));
    }

    private static void formFrameshift(final ProteinContext proteinContext, final StringBuilder sb)
    {
        // report first changed AA (ie the ref) downstream and its index
        int i = 0;
        for(; i < proteinContext.RefAminoAcids.length(); ++i)
        {
            if(i >= proteinContext.AltAminoAcids.length())
                break;

            if(proteinContext.RefAminoAcids.charAt(i) != proteinContext.AltAminoAcids.charAt(i))
                break;
        }

        if(i == proteinContext.RefAminoAcids.length())
            --i;

        sb.append(convertToTriLetters(proteinContext.RefAminoAcids.charAt(i)));
        sb.append(proteinContext.CodonIndex + i);
        sb.append(HGVS_FRAMESHIFT);
    }

    private static void formStopLost(final ProteinContext proteinContext, final StringBuilder sb)
    {
        // - stop_lost - p.Ter494ext*? - don't show any new AAs after the stop since considered irrelevant
        sb.append(convertToTriLetters(proteinContext.RefAminoAcids.charAt(0)));
        sb.append(proteinContext.CodonIndex);

        if(proteinContext.RefAminoAcids.length() == proteinContext.AltAminoAcids.length())
        {
            int lastAaIndex = proteinContext.AltAminoAcids.length() - 1;
            sb.append(convertToTriLetters(proteinContext.AltAminoAcids.charAt(lastAaIndex)));
        }
        else
        {
            // use frameshift if the first different AA is not the stop
            for(int i = 0; i < min(proteinContext.RefAminoAcids.length(), proteinContext.AltAminoAcids.length()); ++i)
            {
                if(proteinContext.RefAminoAcids.charAt(i) != proteinContext.AltAminoAcids.charAt(i))
                {
                    if(proteinContext.AltAminoAcids.charAt(i) != proteinContext.AltAminoAcids.charAt(0))
                        sb.append(convertToTriLetters(proteinContext.AltAminoAcids.charAt(i)));

                    sb.append(HGVS_FRAMESHIFT);
                    return;
                }
            }
        }

        sb.append(HGVS_STOP_LOST);
    }

    private static void formStartLost(final ProteinContext proteinContext, final StringBuilder sb)
    {
        // no other details are required
        sb.append("Met1?");
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
