package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sv.LineElements.LINE_POLY_AT_REQ;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_FRAGS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.esvee.assembly.SequenceCompare.compareSequencesMismatchNoQuals;
import static com.hartwig.hmftools.esvee.caller.LineChecker.hasLineSequence;
import static com.hartwig.hmftools.esvee.caller.Variant.hasLength;
import static com.hartwig.hmftools.esvee.common.SvConstants.isSbx;

import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.sv.SvVcfTags;

import htsjdk.variant.variantcontext.Genotype;

public final class SeqTechUtils
{
    // filters

    // SBX-specific filters
    public static final int SBX_STRAND_BIAS_NON_BND_MIN_FRAGS = 10;
    public static final int SBX_HEURISTIC_THRESHOLD = 50;
    public static final int SBX_HEURISTIC_SGL_THRESHOLD = 30;

    public static final int SBX_HEURISTIC_SHORT_LENGTH = 300;
    public static final double SBX_HEURISTIC_ASM_LENGTH_FACTOR = 600;
    public static final int SBX_HEURISTIC_INEXACT_HOM_FACTOR = 2;
    public static final int SBX_HEURISTIC_INEXACT_HOM_MAX = 50;

    public static final int SBX_HEURISTIC_DUP_INS_LENGTH = 40;
    public static final int SBX_HEURISTIC_SHORT_DUP_INS_PENALTY = 20;

    public static final int SBX_HEURISTIC_INV_SHORT_LENGTH = 100;
    public static final int SBX_HEURISTIC_INV_INS_LENGTH = 20;

    public static final int SBX_HEURISTIC_LOCATION_JUNC_PENALTY = 30;
    public static final int SBX_HEURISTIC_BND_LOCATION_JUNC_PENALTY = 20;

    public static final int SBX_HEURISTIC_BND_INS_PENALTY = 20;
    public static final int SBX_HEURISTIC_BND_LINE_BONUS = 20;
    public static final int SBX_HEURISTIC_SGL_LINE_BONUS = 30;

    public static final double SBX_HEURISTIC_BND_LINE_PERC = 0.9;

    public static final List<String> SBX_INV_INSERT_MOTIFS = List.of(
            "ACCATCAATCCGGATGTATGCCGGATTGATGGT",
            "AGGAGTAACATCCATGTATGGGATGTTACTCCT",
            "ACCCGTAGTATTAATGTATGTAATACTACGGGT",
            "ACTATTGAAGGCTATGTATGAGCCTTCAATAGT",
            "AGGATTATGGCAAATGTATGTTGCCATAATCCT",
            "AGGCGAGTATTCCATACATGGAATACTCGCCT",
            "AGGTTATCAGCTTATGTATGAAGCTGATAACCT",
            "AGGGACATTACCAATGTATGTGGTAATGTCCCT",
            "ACAACAGCCGAAGATGTATGCTTCGGCTGTTGT");

    public static final int SBX_INV_INSERT_MOTIF_MIN_LENGTH = 29;
    public static final int SBX_INV_INSERT_MOTIF_MAX_LENGTH = 34;
    public static final int SBX_INV_INSERT_MOTIF_MAX_DIFFS = 3;

    protected static boolean isSbxStrandBias(final Variant var)
    {
        if(!isSbx())
            return false;

        if(var.type() != BND && var.type() != INV && var.type() != StructuralVariantType.SGL)
            return false;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(var.breakends()[se] == null)
                continue;

            Breakend breakend = var.breakends()[se];

            int splitFragments = 0;
            double maxStrandBias = 0;

            for(Genotype genotype : breakend.Context.getGenotypes())
            {
                splitFragments += breakend.fragmentCount(genotype, SPLIT_FRAGS);

                double strandBias = getGenotypeAttributeAsDouble(genotype, SvVcfTags.STRAND_BIAS, 0.5);
                double adjStrandBias = strandBias > 0.5 ? 1 - strandBias : strandBias;
                maxStrandBias = max(adjStrandBias, maxStrandBias);
            }

            if(maxStrandBias > 0)
                return false;

            if(var.type() != BND && splitFragments < SBX_STRAND_BIAS_NON_BND_MIN_FRAGS)
                return false;
        }

        return true;
    }

    protected static boolean isSbxArtefact(final Variant var)
    {
        if(!isSbx())
            return false;

        if(var.isHotspot())
            return false;

        double maxStrandBias = 0;
        int fullAssemblyLength = var.contextStart().getAttributeAsInt(ASM_LENGTH, 0);
        int inexactHomLength = 0;
        int localJunctionCount = 0;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(var.breakends()[se] == null)
                continue;

            Breakend breakend = var.breakends()[se];

            if(breakend.closeToOriginalJunction())
                ++localJunctionCount;

            inexactHomLength = max(inexactHomLength, breakend.InexactHomology.length());

            for(Genotype genotype : breakend.Context.getGenotypes())
            {
                double strandBias = getGenotypeAttributeAsDouble(genotype, SvVcfTags.STRAND_BIAS, 0.5);
                double adjStrandBias = strandBias > 0.5 ? 1 - strandBias : strandBias;
                maxStrandBias = max(adjStrandBias, maxStrandBias);
            }
        }

        boolean nonLineStranded = !var.isLineSite() && maxStrandBias == 0;

        double heuristicValue = var.qual() * min(fullAssemblyLength / SBX_HEURISTIC_ASM_LENGTH_FACTOR, 1);

        double heuristicThreshold = SBX_HEURISTIC_THRESHOLD;

        double inexactHomPenalty = min(inexactHomLength * SBX_HEURISTIC_INEXACT_HOM_FACTOR, SBX_HEURISTIC_INEXACT_HOM_MAX);

        int varLength = 0;

        if(hasLength(var.type()))
        {
            varLength = var.svLength();

            if(varLength > SBX_HEURISTIC_SHORT_LENGTH)
                return false;
        }

        switch(var.type())
        {
            case DEL:
                if(nonLineStranded)
                    return true;

                heuristicValue -= inexactHomPenalty;

                if(localJunctionCount < 2)
                    heuristicValue -= SBX_HEURISTIC_LOCATION_JUNC_PENALTY;

                if(var.isLineSite())
                    heuristicValue += SBX_HEURISTIC_SGL_LINE_BONUS;

                break;

            case DUP:
            case INS:
                if(nonLineStranded)
                    return true;

                heuristicValue -= inexactHomPenalty;

                if(localJunctionCount < 2)
                    heuristicValue -= SBX_HEURISTIC_LOCATION_JUNC_PENALTY;

                if(var.isLineSite())
                {
                    heuristicValue += SBX_HEURISTIC_SGL_LINE_BONUS;
                }
                else
                {
                    if(varLength < SBX_HEURISTIC_DUP_INS_LENGTH)
                        heuristicValue -= SBX_HEURISTIC_SHORT_DUP_INS_PENALTY;
                }

                break;

            case INV:
                if(varLength == 0)
                    return true;

                if(varLength < SBX_HEURISTIC_INV_SHORT_LENGTH && var.insertSequence().length() > SBX_HEURISTIC_INV_INS_LENGTH)
                    return true;

                break;

            case BND:

                if(localJunctionCount < 2)
                    heuristicValue -= SBX_HEURISTIC_BND_LOCATION_JUNC_PENALTY;

                boolean isLineInsert = hasLineSequence(var.insertSequence(), Orientation.FORWARD)
                        || hasLineSequence(var.insertSequence(), Orientation.REVERSE);

                if(isLineInsert)
                    heuristicValue += SBX_HEURISTIC_BND_LINE_BONUS;
                else if(var.insertSequence().length() >= LINE_POLY_AT_REQ)
                    heuristicValue -= SBX_HEURISTIC_BND_INS_PENALTY;

                break;

            case SGL:
                if(nonLineStranded)
                    return true;

                heuristicThreshold = SBX_HEURISTIC_SGL_THRESHOLD;

                if(localJunctionCount == 0)
                    heuristicValue -= SBX_HEURISTIC_LOCATION_JUNC_PENALTY;

                if(var.isLineSite())
                    heuristicValue += SBX_HEURISTIC_SGL_LINE_BONUS;

                break;
        }

        return heuristicValue < heuristicThreshold;
    }

    protected static boolean isSbxZeroLengthInversion(final Variant var)
    {
        if(!isSbx() || var.type() != INV)
            return false;

        if(var.svLength() != 0)
            return false;

        String insertSequence = var.insertSequence();

        if(insertSequence.length() == 0)
            return true;

        if(insertSequence.length() < SBX_INV_INSERT_MOTIF_MIN_LENGTH || insertSequence.length() > SBX_INV_INSERT_MOTIF_MAX_LENGTH)
            return false;

        String insSeqReversed = Nucleotides.reverseComplementBases(insertSequence);

        int minSequenceDiffs = -1;

        for(String motifSeq : SBX_INV_INSERT_MOTIFS)
        {
            int lengthDiff = insertSequence.length() - motifSeq.length();

            for(int i = 0; i < abs(lengthDiff); ++i)
            {
                int s1Offset = lengthDiff > 0 ? i : 0;
                int s2Offset = lengthDiff < 0 ? i : 0;
                int diffs = sequenceMismatches(insertSequence, motifSeq, s1Offset, s2Offset);
                diffs = min(diffs, sequenceMismatches(insSeqReversed, motifSeq, s1Offset, s2Offset));

                if(minSequenceDiffs < 0 || diffs < minSequenceDiffs)
                    minSequenceDiffs = diffs;
            }
        }

        return minSequenceDiffs <= SBX_INV_INSERT_MOTIF_MAX_DIFFS;
    }

    protected static int sequenceMismatches(final String s1, final String s2, int s1Offset, int s2Offset)
    {
        int s1IndexStart = s1Offset;
        int s1IndexEnd = s1.length() - 1;
        int s2IndexStart = s2Offset;
        int s2IndexEnd = s2.length() - 1;

        return compareSequencesMismatchNoQuals(
                s1.getBytes(), s1IndexStart, s1IndexEnd, s2.getBytes(), s2IndexStart, s2IndexEnd, SBX_INV_INSERT_MOTIF_MAX_DIFFS);
    }
}
