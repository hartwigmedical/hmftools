package com.hartwig.hmftools.common.gene;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import java.util.ArrayList;
import java.util.List;

public class TranscriptUtils
{
    public static int codingBasesToPhase(int codingBases, int startPhase)
    {
        if(startPhase == PHASE_2)
            codingBases += 1;
        else if(startPhase == PHASE_0)
            codingBases += 2;

        // coding base: 1 2 3 4 5 6 7 8 9 10
        // phase:       1 2 0 1 2 0 1 2 0
        if(codingBases >= 0)
            return codingBases % 3;

        // coding base: -6 -5 -4 -3 -2 -1  0
        // phase:        0  1  2  0  1  2  0

        int mod = abs(codingBases) % 3;

        if(mod == 0)
            return 0;

        return mod == 1 ? 2 : 1;
    }

    public static int tickPhaseForward(int phase, int increment) { return increment != 0 ? (phase + increment) % 3 : phase; }
    public static int tickPhaseForward(int phase) { return (phase + 1) % 3; }

    public static int calcExonicCodingPhase(final ExonData exon, int codingStart, int codingEnd, byte strand, int position)
    {
        // if coding has started in an earlier exon, then take the exon start phase and process new coding bases
        // ie if previous exon's end phase was 2, then start phase of this exon will also be 2, and first base will tick to phase 0
        // if coding begins with the first base, then if phase is not specified, then assume '1',
        // otherwise use the starting phase to adjust the exonic phase ie:
        // PredEndPhase = if ExonPhase < 0 then CodingBases %% 3 else ExonPhase+CodingBases %% 3

        if(!positionWithin(position, codingStart, codingEnd))
            return PHASE_NONE;

        // coding bases up to and including the specified position
        int codingBases = strand == POS_STRAND ? position - max(exon.Start, codingStart) + 1 : min(exon.End, codingEnd) - position + 1;

        if((strand == POS_STRAND && codingStart == exon.Start) || (strand == NEG_STRAND && codingEnd == exon.End))
        {
            // coding starts on the first base of this exon - take the initial phase and tick it forward
            // phasing starts on 1 if not specified (ie is -1)
            int startPhase = exon.PhaseStart == PHASE_NONE ? PHASE_0 : exon.PhaseStart;
            return (startPhase + codingBases) % 3;
        }
        else if((strand == POS_STRAND && codingStart < exon.Start) || (strand == NEG_STRAND && codingEnd > exon.End))
        {
            // coding started before this exon
            return (exon.PhaseStart + codingBases) % 3;
        }
        else
        {
            // coding starts within the exon - assumes start with phase 1
            return codingBases % 3;
        }
    }

    public static int codingBaseLength(final TranscriptData transData)
    {
        if(transData.nonCoding())
            return 0;

        int codingBases = 0;

        for(ExonData exon : transData.exons())
        {
            if(transData.CodingStart > exon.End)
                continue;

            if(transData.CodingEnd < exon.Start)
                break;

            codingBases += min(exon.End, transData.CodingEnd) - max(exon.Start, transData.CodingStart) + 1;
        }

        return codingBases;
    }

    public static int calcCodingStartPositionAdjustment(final TranscriptData transData, final ExonData exon)
    {
        if(transData.nonCoding())
            return 0;

        int exonPosition = transData.Strand == POS_STRAND ? exon.Start : exon.End;

        if(!positionWithin(exonPosition, transData.CodingStart, transData.CodingEnd))
            return 0;

        // finds the phased coding adjustment for from the first coding base
        int startPhase = calcExonicCodingPhase(exon, transData.CodingStart, transData.CodingEnd, transData.Strand, exonPosition);

        if(startPhase == PHASE_2)
            return transData.Strand == POS_STRAND ? 2 : -2;
        else if(startPhase == PHASE_0)
            return transData.Strand == POS_STRAND ? 1 : -1;
        else
            return 0; // already on the codon start
    }

    public static CodingBaseData calcCodingBases(final TranscriptData transData, int position)
    {
        // intronic phase can read from the preceding exon end phase (regardless of whether coding has begun
        // exonic phase depends on where coding begins - see function for details
        boolean inCodingRegion = false;
        boolean codingRegionEnded = false;

        CodingBaseData cbData = new CodingBaseData();

        if(transData.nonCoding())
            return cbData;

        int codingStart = transData.CodingStart;
        int codingEnd = transData.CodingEnd;

        boolean posInCodingRegion = positionWithin(position, codingStart, codingEnd);

        boolean isExonic = false;

        for (int i = 0; i < transData.exons().size(); ++i)
        {
            final ExonData exon = transData.exons().get(i);
            final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

            int exonStart = exon.Start;
            int exonEnd = exon.End;

            if(!isExonic && positionWithin(position, exonStart, exonEnd))
            {
                isExonic = true;

                // set phase
                if(posInCodingRegion)
                    cbData.Phase = calcExonicCodingPhase(exon, codingStart, codingEnd, transData.Strand, position);
            }

            // set phase for intronic positions
            if(position > exonEnd && nextExon != null && position < nextExon.Start)
            {
                if(transData.Strand == POS_STRAND)
                    cbData.Phase = exon.PhaseEnd;
                else
                    cbData.Phase = nextExon.PhaseEnd;
            }

            if (!inCodingRegion)
            {
                if (exonEnd >= codingStart)
                {
                    // coding region begins in this exon
                    inCodingRegion = true;

                    cbData.TotalCodingBases += exonEnd - codingStart + 1;

                    // check whether the position falls in this exon and if so before or after the coding start
                    if (position >= codingStart)
                    {
                        if (position < exonEnd)
                            cbData.CodingBases += position - codingStart + 1;
                        else
                            cbData.CodingBases += exonEnd - codingStart + 1;
                    }
                }
            }
            else if (!codingRegionEnded)
            {
                if (exonStart > codingEnd)
                {
                    codingRegionEnded = true;
                }
                else if (exonEnd >= codingEnd)
                {
                    // coding region ends in this exon
                    codingRegionEnded = true;

                    cbData.TotalCodingBases += codingEnd - exonStart + 1;

                    if (position >= exonStart)
                    {
                        if (position < codingEnd)
                            cbData.CodingBases += position - exonStart + 1;
                        else
                            cbData.CodingBases += codingEnd - exonStart + 1;
                    }
                }
                else
                {
                    // take all of the exon's bases
                    cbData.TotalCodingBases += exonEnd - exonStart + 1;

                    if (position >= exonStart)
                    {
                        if (position < exonEnd)
                            cbData.CodingBases += position - exonStart + 1;
                        else
                            cbData.CodingBases += exonEnd - exonStart + 1;
                    }
                }

                if(codingRegionEnded)
                    break;
            }
        }

        // adjust for strand of the gene
        if(transData.Strand == NEG_STRAND)
        {
            if(posInCodingRegion && isExonic)
                cbData.CodingBases = cbData.TotalCodingBases - cbData.CodingBases + 1;
            else
                cbData.CodingBases = cbData.TotalCodingBases - cbData.CodingBases;
        }

        if(!posInCodingRegion)
            cbData.Phase = PHASE_NONE;

        return cbData;
    }

    public static List<int[]> getCodingBaseRanges(final TranscriptData transData, int startPosition, boolean searchUp, int requiredBases)
    {
        // returns ranges of coding bases starting from a given position for a specified number of bases up or down from there
        final List<int[]> codingBaseRanges = new ArrayList<>();

        if(transData.nonCoding())
            return codingBaseRanges;

        if(startPosition < transData.CodingStart || startPosition > transData.CodingEnd)
            return codingBaseRanges;

        final List<ExonData> exonDataList = transData.exons();

        if(searchUp)
        {
            for (int i = 0; i < exonDataList.size(); ++i)
            {
                final ExonData exon = exonDataList.get(i);

                if(startPosition > exon.End)
                    continue;

                if(exon.Start > transData.CodingEnd)
                    break; // no more coding bases

                int exonBaseStart = max(exon.Start, startPosition);
                int exonBaseEnd = min(transData.CodingEnd, exon.End);

                int exonBaseCount = exonBaseEnd - exonBaseStart + 1;

                if(requiredBases >= exonBaseCount)
                {
                    // take them all
                    requiredBases -= exonBaseCount;
                }
                else
                {
                    exonBaseEnd = exonBaseStart + requiredBases - 1;
                    requiredBases = 0;
                }

                codingBaseRanges.add(new int[] {exonBaseStart, exonBaseEnd});

                if (requiredBases <= 0)
                    break;
            }
        }
        else
        {
            for(int i = exonDataList.size() - 1; i >= 0; --i)
            {
                final ExonData exon = exonDataList.get(i);

                if(startPosition < exon.Start)
                    continue;

                if(exon.End < transData.CodingStart)
                    break;

                int exonBaseEnd = min(exon.End, startPosition);
                int exonBaseStart = max(transData.CodingStart, exon.Start);

                int exonBaseCount = exonBaseEnd - exonBaseStart + 1;

                if(requiredBases >= exonBaseCount)
                {
                    requiredBases -= exonBaseCount;
                }
                else
                {
                    exonBaseStart = exonBaseEnd - requiredBases + 1;
                    requiredBases = 0;
                }

                codingBaseRanges.add(0, new int[] {exonBaseStart, exonBaseEnd});

                if (requiredBases <= 0)
                    break;
            }
        }

        return codingBaseRanges;
    }
}
