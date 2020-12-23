package com.hartwig.hmftools.common.fusion;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;

import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;

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

    public static int tickPhaseForward(int phase, int increment) { return (phase + increment) % 3; }
    public static int tickPhaseForward(int phase) { return (phase + 1) % 3; }

    /*
    public static int calcPositionPhasing(final Transcript transcript, boolean isUpstream)
    {
        int codingBases = transcript.codingBases();

        // factor in insert sequence for the upstream partner
        if(isUpstream && !transcript.gene().insertSequence().isEmpty())
        {
            codingBases += transcript.gene().insertSequence().length();
        }

        return codingBasesToPhase(codingBases);
    }
    */

    public static CodingBaseData calcCodingBases(final TranscriptData transData, int position)
    {
        boolean inCodingRegion = false;
        boolean codingRegionEnded = false;

        CodingBaseData cbData = new CodingBaseData();

        if(transData.CodingStart == null)
            return cbData;

        int codingStart = transData.CodingStart;
        int codingEnd = transData.CodingEnd;

        // by convention unless coding starts in the first base of an exon
        int startPhase = PHASE_1;
        boolean isExonic = false;

        for (ExonData exon : transData.exons())
        {
            int exonStart = exon.Start;
            int exonEnd = exon.End;

            if(!isExonic)
                isExonic = positionWithin(position, exonStart, exonEnd);

            if (!inCodingRegion)
            {
                if (exonEnd >= codingStart)
                {
                    // coding region begins in this exon
                    inCodingRegion = true;

                    if(transData.Strand == POS_STRAND && exon.Start == codingStart)
                        startPhase = tickPhaseForward(exon.PhaseStart);

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

                    if(transData.Strand == NEG_STRAND && exonEnd == codingEnd)
                        startPhase = tickPhaseForward(exon.PhaseStart);

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
            }
        }

        // adjust for strand of the gene
        boolean posInCodingRegion = positionWithin(position, codingStart, codingEnd);

        if(transData.Strand == NEG_STRAND)
        {
            if(posInCodingRegion && isExonic)
                cbData.CodingBases = cbData.TotalCodingBases - cbData.CodingBases + 1;
            else
                cbData.CodingBases = cbData.TotalCodingBases - cbData.CodingBases;
        }

        if(posInCodingRegion)
            cbData.Phase = codingBasesToPhase(cbData.CodingBases, startPhase);
        else
            cbData.Phase = PHASE_NONE;

        return cbData;
    }

}
