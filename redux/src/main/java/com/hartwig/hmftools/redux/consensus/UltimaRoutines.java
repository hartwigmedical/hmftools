package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.maxQual;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.minQual;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.HALF_PHRED_SCORE_SCALING;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULT_QUAL_TAG_DELIM;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractT0Values;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractTpValues;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.isHighBaseQual;
import static com.hartwig.hmftools.redux.consensus.BaseQualPair.NO_BASE;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;

import htsjdk.samtools.SAMRecord;

public final class UltimaRoutines
{
    public static void finaliseRead(final RefGenomeInterface refGenome, final SAMRecord record)
    {
        ConsensusType consensusType = UltimaBamUtils.deriveConsensusType(record);
        record.setAttribute(CONSENSUS_TYPE_ATTRIBUTE, consensusType.toString());
    }

    public static void setLowQualTag(final SAMRecord record)
    {
        /*
        byte[] tpValues = extractTpValues(record);
        byte[] t0Values = extractT0Values(record);

        byte[] bases = record.getReadBases();
        byte[] baseQuals = record.getBaseQualities();

        StringJoiner sj = null;

        byte hpBase = NO_BASE;
        byte hpQual = -1;
        int hpStartIndex = -1;

        for(int i = 0; i < baseQuals.length; ++i)
        {
            byte base = bases[i];
            byte qual = baseQuals[i];

            byte minHpQual = -1;

            // deletions
            byte t0Qual = t0Values[i];

            if(t0Qual > 0)
                minHpQual = t0Qual;

            // homopolymer adjustments
            byte tpQual = -1;

            if(hpBase != NO_BASE)
            {
                if(hpBase == base)
                {
                    // HP continues
                    tpQual = hpQual;
                }
                else
                {
                    // apply last HP qual then clear
                    int hpLength = i - hpStartIndex + 1;

                    if(belowLowQualThreshold(hpLength, hpQual))
                    {
                        if(sj == null)
                            sj = new StringJoiner(ULT_QUAL_TAG_DELIM);

                        if(hpLength > 1)
                            sj.add(format("%d-%d", hpStartIndex, i));
                        else
                            sj.add(String.valueOf(i));
                    }
                }
            }

            if(!isHighBaseQual(qual))
            {
                hpQual = qual;
                tpQual = qual;

                if(i < baseQuals.length - 1)
                {
                    if(bases[i + 1] == base)
                    {
                        // HP length > 1
                        if(tpQual - HALF_PHRED_SCORE_SCALING < 0)
                            tpQual = 0;
                        else
                            tpQual -= HALF_PHRED_SCORE_SCALING;
                    }
                }
            }

            if(tpQual >= 0 && t0Qual >= 0)
                minHpQual = minQual(t0Qual, tpQual);
            else if(tpQual >= 0)
                minHpQual = tpQual;
            if(t0Qual >= 0)
                minHpQual = t0Qual;
        }
        */
    }

    private static boolean belowLowQualThreshold(int homopolymerLength, byte qual)
    {
        if(homopolymerLength == 1)
            return qual < 25;
        else if(homopolymerLength == 2)
            return qual < 22;
        else if(homopolymerLength == 3)
            return qual < 19;
        else if(homopolymerLength <= 6)
            return qual < 16;
        else if(homopolymerLength < 8)
            return qual < 14;
        else
            return qual < 12;
    }
}
