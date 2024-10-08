package com.hartwig.hmftools.chord.indel;

import com.hartwig.hmftools.chord.common.SmallVariant;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import org.jetbrains.annotations.NotNull;

public class IndelVariant
{
    public final SmallVariant mVariant;

    public final int mIndelLength;
    public final boolean mIsDeletion;

    public IndelVariant(SmallVariant variant)
    {
        if(!variant.isIndel())
            throw new IllegalStateException(variant + " is not an indel");

        mVariant = variant;

        int signedLength = mVariant.AltBases.length() - mVariant.RefBases.length();

        mIsDeletion = signedLength<0;
        mIndelLength = Math.abs(signedLength);
    }

    public String getFlankingRefBasesRight(@NotNull RefGenomeSource refGenome, int indelLengthUnits)
    {
        int start;
        int end;

        if(mIsDeletion)
        {
            start = mVariant.Position + 1 + mIndelLength;
            end = mVariant.Position + mIndelLength + mIndelLength*indelLengthUnits;
        }
        else
        {
            start = mVariant.Position + 1;
            end = mVariant.Position + mIndelLength*indelLengthUnits;
        }

        String chromosome = mVariant.Chromosome.toString();

        end = Math.min(end, refGenome.getChromosomeLength(chromosome));

        return new String(refGenome.getBases(chromosome, start, end));
    }

    public String getFlankingRefBasesLeft(@NotNull RefGenomeSource refGenome, int indelLengthUnits)
    {
        int start = mVariant.Position - mIndelLength*indelLengthUnits + 1;
        int end = mVariant.Position;

        return new String(refGenome.getBases(mVariant.Chromosome.toString(), start, end));
    }

    public String getIndelSequence()
    {
        String indelSequence;
        if(mIsDeletion)
            indelSequence = mVariant.RefBases.substring(1);
        else
            indelSequence = mVariant.AltBases.substring(1);

        return indelSequence;
    }

    public static int countHomologousBases(String indelSequence, String flankSequence)
    {
        int nBasesToScan = Math.min(indelSequence.length(), flankSequence.length());

        int nHomologousBases = 0;

        for(int baseIndex = 0; baseIndex<nBasesToScan; baseIndex++)
        {
            char indelBase = indelSequence.charAt(baseIndex);
            char flankBase = flankSequence.charAt(baseIndex);

            if(indelBase==flankBase)
                nHomologousBases++;
            else
                break;
        }

        return nHomologousBases;
    }

    private static final char REPEAT_UNIT_PLACEHOLDER = '#';

    public static int countRepeatUnits(String indelSequence, String rightFlankSequence)
    {
        // TODO: Replace indel sequence in flank sequence with Zs. Then count number of left aligned Zs

        // The indel sequence itself counts as one repeat
        String sequence = indelSequence + rightFlankSequence;

        String sequenceWithPlaceholder = sequence.replaceAll(indelSequence, String.valueOf(REPEAT_UNIT_PLACEHOLDER));
        int repeatUnitCount = 0;
        for(int i = 0; i < sequenceWithPlaceholder.length(); i++)
        {
            if(sequenceWithPlaceholder.charAt(i) == REPEAT_UNIT_PLACEHOLDER)
                repeatUnitCount++;
            else
                break;
        }

        return repeatUnitCount;
    }
}
