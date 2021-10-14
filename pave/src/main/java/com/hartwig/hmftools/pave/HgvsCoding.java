package com.hartwig.hmftools.pave;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseStrandBases;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.NON_CODING;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;

public final class HgvsCoding
{
    /* Rules and conventions
    - coding bases for -ve strand is reversed bases
    - start codon start is position 1
    - nucleotide -1 is the 1st base prior to the start codon
    - nucleotide *1 is the 1st base 3’ of the translation stop codon
    - use most 3’ position in case of homology for duplications and deletions
    - phased inframe variants should get the combined impact
    - building routine
        1. c. or n.
        2. - or * if 5' or 3' UTR
        3. coding base
        4. distance into intron (- if previous, + if next) if intronic
        3 & 4 repeated if more than 1 base variant
        5. mutation details - del, dup, ins or combo
        6. bases, reversed if -ve strand
    - for multi-base variants, distance into intron goes up if positive (ie into next intron), otherwise down
    - insertion as classified as duplication if microhomology matches the inserted bases
     */

    private static final String CODING_ID = "c.";
    private static final String NON_CODING_ID = "n.";

    public static final String HGVS_TYPE_DEL = "del";
    public static final String HGVS_TYPE_DUP = "dup";
    public static final String HGVS_TYPE_INS = "ins";

    public static final String HGVS_UNKNOWN = "unknown"; // misclassified or malformed string

    public static void set(final VariantData variant, final CodingContext codingContext)
    {
        codingContext.Hgvs = generate(variant, codingContext);
    }

    public static String generate(final VariantData variant, final CodingContext codingContext)
    {
        if(codingContext.RegionType == UPSTREAM) // undefined by the standard
            return "";

        StringBuilder sb = new StringBuilder();

        addTranscriptType(codingContext, sb);

        if(variant.isBaseChange())
        {
            formPointMutation(variant, codingContext, sb);
        }
        else if(variant.isDeletion())
        {
            formDeletion(variant, codingContext, sb);
        }
        else if(variant.isInsert())
        {
            formInsertion(variant, codingContext, sb);
        }

        return sb.toString();
    }

    private static void addTranscriptType(final CodingContext codingContext, final StringBuilder sb)
    {
        if(codingContext.CodingType == NON_CODING)
            sb.append(NON_CODING_ID);
        else
            sb.append(CODING_ID);
    }

    private static void addUtr(final CodingContext codingContext, final StringBuilder sb)
    {
        if(codingContext.CodingType == UTR_5P && codingContext.CodingBase > 1)
            sb.append('-');
        else if(codingContext.CodingType == UTR_3P)
            sb.append('*');
    }

    private static void addCodingBase(final CodingContext codingContext, int codingBase, final StringBuilder sb, boolean isSecond)
    {
        if(isSecond)
            sb.append('_');

        addUtr(codingContext, sb);

        if((codingContext.CodingType == UTR_5P || codingContext.CodingType == UTR_3P) && codingContext.CodingBase == 0)
            codingBase = 1;

        sb.append(codingBase);
    }

    private static void addIntronicPosition(int position, final StringBuilder sb)
    {
        if(position > 0)
            sb.append('+');
        else
            sb.append('-');

        sb.append(abs(position));
    }

    private static void addBases(final String bases, final CodingContext codingContext, final StringBuilder sb)
    {
        if(codingContext.Strand == POS_STRAND)
            sb.append(bases);
        else
            sb.append(reverseStrandBases(bases));
    }

    private static void formDeletion(final VariantData variant, final CodingContext codingContext, final StringBuilder sb)
    {
        // examples:
        // c.191delG
        // c.203_204delCG
        // c.1208_1210delAGA - shows the 3 deleted bases and their coding positions
        // c.60+5601_60+5610delTTTTTTTTTT - 10 bases from intron (CTTTTTTTTTT>C)
        // c.123-7190_123-7175delTGTGTGTGTGTGTGTG (+ve strand) neg distance drops the intronic position

        int delLength = abs(variant.baseDiff());
        int codingBase = codingContext.CodingBase;
        int nearestExon = codingContext.NearestExonDistance;

        // move past ref base
        if(codingContext.RegionType == EXONIC)
            ++codingBase;

        ++nearestExon;

        addCodingBase(codingContext, codingBase, sb, false);

        if(codingContext.RegionType == INTRONIC)
            addIntronicPosition(nearestExon, sb);

        if(delLength > 1)
        {
            if(codingContext.RegionType == EXONIC)
                codingBase += delLength - 1;

            addCodingBase(codingContext, codingBase, sb, true);

            if(codingContext.RegionType == INTRONIC)
            {
                nearestExon += delLength - 1;
                addIntronicPosition(nearestExon, sb);
            }
        }

        String deletedBases = variant.Ref.substring(1);
        sb.append(HGVS_TYPE_DEL);
        addBases(deletedBases, codingContext, sb);
    }

    public static boolean isDuplication(final VariantData variant)
    {
        return variant.isInsert() && variant.microhomology().equals(variant.Alt.substring(1));
    }

    private static void formInsertion(final VariantData variant, final CodingContext codingContext, final StringBuilder sb)
    {
        // coding: c.1033_1034insA
        // intronic post-exon: c.15+1619_15+1620insTTTGTT
        // 5'UTR: c.-23-304_-23-303insA
        // c.77-100_77-99insACACAC neg intronic distance from (T>TACACAC)

        if(isDuplication(variant))
        {
            formDuplication(variant, codingContext, sb);
            return;
        }

        String insertedBases = variant.Alt.substring(1);
        int codingBase = codingContext.CodingBase;
        int nearestExon = codingContext.NearestExonDistance;

        addCodingBase(codingContext, codingBase, sb, false);

        if(codingContext.RegionType == INTRONIC)
            addIntronicPosition(nearestExon, sb);

        if(codingContext.RegionType == EXONIC)
            codingBase += 1;

        addCodingBase(codingContext, codingBase, sb, true);

        if(codingContext.RegionType == INTRONIC)
        {
            nearestExon += 1;
            addIntronicPosition(nearestExon, sb);
        }

        sb.append(HGVS_TYPE_INS);
        addBases(insertedBases, codingContext, sb);
    }

    private static void formDuplication(final VariantData variant, final CodingContext codingContext, final StringBuilder sb)
    {
        // c.377+9769_377+9772dupAAAT
        // c.466+18711_466+18712dupGA
        // c.1440-82dupT - single location if a single duplicated base
        // 	n.210-64739_210-64734dupGTGTGT

        int codingBase = codingContext.CodingBase;
        int nearestExon = codingContext.NearestExonDistance;
        String insertedBases = variant.Alt.substring(1);
        int insertLength = insertedBases.length();

        if(insertLength == 1)
        {
            addCodingBase(codingContext, codingBase, sb, false);

            if(codingContext.RegionType == INTRONIC)
                addIntronicPosition(nearestExon, sb);
        }
        else
        {
            int baseShift = insertLength - 1;

            // on the negative strand the positions are flipped around
            int codingBaseStart;
            int codingBaseEnd;
            int intronBaseStart;
            int intronBaseEnd;

            if(codingContext.Strand == POS_STRAND)
            {
                codingBaseStart = codingBase;
                codingBaseEnd = codingContext.RegionType == EXONIC ? codingBase + baseShift : codingBase;

                intronBaseStart = nearestExon;
                intronBaseEnd = codingContext.RegionType == INTRONIC ? nearestExon + baseShift : nearestExon;
            }
            else
            {
                codingBaseEnd = codingBase;
                codingBaseStart = codingContext.RegionType == EXONIC ? codingBase - baseShift : codingBase;

                intronBaseEnd = nearestExon;
                intronBaseStart = codingContext.RegionType == INTRONIC ? nearestExon - baseShift : nearestExon;
            }

            if(codingContext.Strand == NEG_STRAND && codingContext.CodingType == UTR_5P && codingContext.RegionType == EXONIC)
            {
                // may apply to both strands
                codingBaseStart = codingBaseEnd + baseShift;
            }

            addCodingBase(codingContext, codingBaseStart, sb, false);

            if(codingContext.RegionType == INTRONIC)
                addIntronicPosition(intronBaseStart, sb);

            addCodingBase(codingContext, codingBaseEnd, sb, true);

            if(codingContext.RegionType == INTRONIC)
            {
                addIntronicPosition(intronBaseEnd, sb);
            }
        }

        sb.append(HGVS_TYPE_DUP);
        addBases(insertedBases, codingContext, sb);
    }

    private static void formPointMutation(final VariantData variant, final CodingContext codingContext, final StringBuilder sb)
    {
        // upstream: c.-14G>C
        // coding: c.76A>C
        // post-exon intronic: c.88+1G>T
        // pre-exon intronic: c.89-2A>C
        // MNVs intronic c.38-20810_38-20808delTTGinsACT
        // MNVs exonic: c.132_133delCCinsTT
        // MNVs 3'UTR: c.*907_*908delCCinsGT
        // MNVs 5'UTR: c.-83_-82delGGinsAA

        int codingBase = codingContext.CodingBase;
        int nearestExon = codingContext.NearestExonDistance;
        int varLength = variant.Ref.length();

        if(varLength == 1)
        {
            addCodingBase(codingContext, codingBase, sb, false);

            if(codingContext.RegionType == INTRONIC)
                addIntronicPosition(nearestExon, sb);

            addBases(variant.Ref, codingContext, sb);
            sb.append(">");
            addBases(variant.Alt, codingContext, sb);
        }
        else
        {
            int codingBaseStart = codingBase;
            int codingBaseEnd = codingContext.RegionType == EXONIC ? codingBase + varLength - 1 : codingBase;

            if(codingContext.CodingType == UTR_5P && codingContext.RegionType == EXONIC)
            {
                // swap and reduce by 1
                codingBaseStart = codingBaseEnd - 1;
                codingBaseEnd = codingBase - 1;
            }

            addCodingBase(codingContext, codingBaseStart, sb, false);

            if(codingContext.RegionType == INTRONIC)
                addIntronicPosition(nearestExon, sb);

            addCodingBase(codingContext, codingBaseEnd, sb, true);

            if(codingContext.RegionType == INTRONIC)
            {
                nearestExon += varLength - 1;
                addIntronicPosition(nearestExon, sb);
            }

            sb.append(HGVS_TYPE_DEL);
            addBases(variant.Ref, codingContext, sb);
            sb.append(HGVS_TYPE_INS);
            addBases(variant.Alt, codingContext, sb);
        }
    }

}
