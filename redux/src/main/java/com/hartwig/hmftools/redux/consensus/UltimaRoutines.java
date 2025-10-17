package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_TYPE_ATTRIBUTE;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.HALF_PHRED_SCORE_SCALING;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULT_QUAL_TAG;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.ULT_QUAL_TAG_DELIM;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.extractT0Values;
import static com.hartwig.hmftools.common.sequencing.UltimaBamUtils.isHighBaseQual;

import java.util.StringJoiner;

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
        setLowQualTag(record);
    }

    public static void setLowQualTag(final SAMRecord record)
    {
        byte[] t0Values = extractT0Values(record);

        byte[] bases = record.getReadBases();
        byte[] baseQuals = record.getBaseQualities();

        StringJoiner sj = null;

        int hpStartIndex = -1;
        int hpEndIndex = -1;
        boolean hpLowQual = false;
        int lowQualStartIndex = -1;

        for(int i = 0; i < baseQuals.length; ++i)
        {
            byte qual = baseQuals[i];

            boolean isLowQual = false;

            // deletions
            if(t0Values[i] < ULTIMA_HP_DELETION_LOW_QUAL_THRESHOLD)
                isLowQual = true;

            // homopolymer adjustments
            if(hpStartIndex >= 0 && i <= hpEndIndex)
            {
                // in an existing homopolymer with low qual
                isLowQual |= hpLowQual;
            }
            else
            {
                // assess a new homopolymer
                if(!isHighBaseQual(qual))
                {
                    byte base = bases[i];

                    // find the length of the HP
                    int hpLength = 1;

                    for(int j = i + 1; j < baseQuals.length; ++j)
                    {
                        if(bases[j] != base)
                            break;

                        ++hpLength;
                    }

                    byte tpQual = qual;

                    if(hpLength > 1)
                    {
                        if(tpQual - HALF_PHRED_SCORE_SCALING < 0)
                            tpQual = 0;
                        else
                            tpQual -= HALF_PHRED_SCORE_SCALING;
                    }

                    hpStartIndex = i;
                    hpEndIndex = i + hpLength - 1;

                    if(belowLowQualThreshold(hpLength, tpQual))
                    {
                        isLowQual = true;
                        hpLowQual = true;
                    }
                }
            }

            if(isLowQual)
            {
                if(lowQualStartIndex < 0)
                    lowQualStartIndex = i;

                if(i == baseQuals.length - 1)
                {
                    if(sj == null)
                        sj = new StringJoiner(ULT_QUAL_TAG_DELIM);

                    addLowQualTagEntry(sj, lowQualStartIndex, baseQuals.length - 1);
                }
            }
            else
            {
                if(lowQualStartIndex >= 0)
                {
                    int lowQualEndIndex = i - 1;

                    if(sj == null)
                        sj = new StringJoiner(ULT_QUAL_TAG_DELIM);

                    addLowQualTagEntry(sj, lowQualStartIndex, lowQualEndIndex);
                    lowQualStartIndex = -1;
                }
            }
        }

        if(sj != null)
            record.setAttribute(ULT_QUAL_TAG, sj.toString());
    }

    private static void addLowQualTagEntry(final StringJoiner sj, int lowQualStartIndex, int lowQualEndIndex)
    {
        if(lowQualEndIndex > lowQualStartIndex)
            sj.add(format("%d-%d", lowQualStartIndex, lowQualEndIndex));
        else
            sj.add(String.valueOf(lowQualStartIndex));
    }

    private static final byte ULTIMA_HP_DELETION_LOW_QUAL_THRESHOLD = 20;

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
