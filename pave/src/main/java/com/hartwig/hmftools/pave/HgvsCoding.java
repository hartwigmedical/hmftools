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
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;

public final class HgvsCoding
{
    /* Rules and conventions
    - coding bases for -ve strand is reversed bases
    - no nucleotide 0 - example?
    - start codon start is position 1
    - nucleotide -1 is the 1st base prior to the start codon
    - nucleotide *1 is the 1st base 3’ of the translation stop codon
    - other rules:
        - use most 3’ position in case of homology for duplications and deletions
        - phased inframe variants should get the combined impact

    - questions
        - MNV treated as del-ins rather than >>
        - example of crossing splice boundary?

    - decision flow and inputs
        - 5' UTR - exonic bases prior to coding start
        - 3' UTR - exonic bases after to coding end (and record total coding bases)
        - upstream - number of bases before transcript start
        - bases before start-codon (starting at 1, denoted as -X)
        - bases after stop-codon (starting at 1, denoted as *X)
        - type of SNV - point, MNV (del-ins?), DEL or INS, Duplication (INS with full homology)
        - based into / before exon (choose lowest)


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
        StringBuilder sb = new StringBuilder();

        if(codingContext.CodingType == NON_CODING)
            sb.append(NON_CODING_ID);
        else
            sb.append(CODING_ID);

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

    private static void formDeletion(final VariantData variant, final CodingContext codingContext, final StringBuilder sb)
    {
        // coding: c.76_78delACT
        // intronic pre-exon: c.726-5537_726-5536delTT
        // 5'UTR c.-147-1093delA

        int nonCodingDistance = codingContext.NearestExonDistance;
        int delLength = abs(variant.baseDiff());

        if(codingContext.CodingType == CODING)
        {
            if(codingContext.RegionType == EXONIC)
            {
                sb.append(codingContext.CodingBase);
                sb.append("_");
                sb.append(codingContext.CodingBase + codingContext.DeletedCodingBases);
            }
            else
            {
                char prePostSign = nonCodingDistance < 0 ? '-' : '+';

                if(delLength > 1)
                {
                    sb.append(codingContext.CodingBase + delLength);
                    sb.append(prePostSign);
                    sb.append(abs(nonCodingDistance));
                    sb.append("_");
                }

                sb.append(codingContext.CodingBase);
                sb.append(prePostSign);
                sb.append(abs(nonCodingDistance));
            }

        }
        else if(codingContext.CodingType == UTR_5P || codingContext.RegionType == UPSTREAM)
        {
            sb.append("-");
            sb.append(abs(nonCodingDistance));
        }
        else if(codingContext.CodingType == UTR_3P)
        {
            sb.append("*");
            sb.append(abs(nonCodingDistance));
        }
        else
        {
            sb.append(HGVS_UNKNOWN);
        }

        sb.append(HGVS_TYPE_DEL);

        String ref = variant.Ref.substring(1);

        if(codingContext.Strand == POS_STRAND)
            ref = reverseStrandBases(ref);

        sb.append(ref);
    }

    private static void formInsertion(final VariantData variant, final CodingContext codingContext, final StringBuilder sb)
    {
        // coding: c.1033_1034insA
        // intronic post-exon: c.15+1619_15+1620insTTTGTT
        // 5'UTR: c.-23-304_-23-303insA


    }

    private static void formPointMutation(final VariantData variant, final CodingContext codingContext, final StringBuilder sb)
    {
        // upstream: c.-14G>C
        // coding: c.76A>C
        // post exon intronic: c.88+1G>T
        // pre exon intronic: c.89-2A>C
        int nonCodingDistance = codingContext.NearestExonDistance;

        if(codingContext.CodingType == CODING)
        {
            sb.append(codingContext.CodingBase);

            if(codingContext.RegionType == INTRONIC)
            {
                char prePostSign = nonCodingDistance < 0 ? '-' : '+';
                sb.append(prePostSign);
                sb.append(abs(nonCodingDistance));
            }
        }
        else
        {
            char prePostSign = (codingContext.CodingType == UTR_5P || codingContext.RegionType == UPSTREAM) ? '-' : '*';
            sb.append(prePostSign);
            sb.append(abs(nonCodingDistance));
        }

        String ref = codingContext.Strand == POS_STRAND ? variant.Ref : reverseStrandBases(variant.Ref);
        String alt = codingContext.Strand == POS_STRAND ? variant.Alt : reverseStrandBases(variant.Alt);

        // TODO - handle MNVs
        sb.append(ref);
        sb.append(">");
        sb.append(alt);
    }

}
