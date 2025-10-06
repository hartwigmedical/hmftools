package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.BASE_QUAL_MINIMUM;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.isUncertainBaseQual;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.maxQual;
import static com.hartwig.hmftools.common.sequencing.SequencingType.ILLUMINA;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;
import static com.hartwig.hmftools.redux.consensus.BaseQualPair.NO_BASE;
import static com.hartwig.hmftools.redux.consensus.ConsensusState.setReadPositionStartOffsets;
import static com.hartwig.hmftools.redux.consensus.IlluminaRoutines.isDualStrandAndIsFirstInPair;

import java.util.List;

import javax.annotation.Nullable;

import com.hartwig.hmftools.common.sequencing.SequencingType;

import htsjdk.samtools.SAMRecord;

public class BaseBuilder
{
    private final RefGenome mRefGenome;
    private final ConsensusStatistics mConsensusStats;
    private final SequencingType mSequencingType;

    // cached for the majority of successive reads being on the same chromosome, to protect against ref genome base requests beyond limits
    private int mChromosomeLength;

    public BaseBuilder(final RefGenome refGenome, final ConsensusStatistics consensusStats, final SequencingType sequencingType)
    {
        mRefGenome = refGenome;
        mConsensusStats = consensusStats;
        mSequencingType = sequencingType;
        mChromosomeLength = 0;
    }

    public void setChromosomLength(int chromosomeLength) { mChromosomeLength = chromosomeLength; }
    public int chromosomeLength() { return mChromosomeLength; }
    public RefGenome refGenome() { return mRefGenome; }
    public SequencingType sequencingType() { return mSequencingType; }

    public static final int INVALID_POSITION = -1;

    public void buildReadBases(final List<SAMRecord> reads, final ConsensusState consensusState)
    {
        int chromosomeLength = mChromosomeLength;
        if(chromosomeLength == 0)
            chromosomeLength = mRefGenome.getChromosomeLength(reads.get(0).getReferenceName());

        int baseLength = consensusState.Bases.length;
        int readCount = reads.size();
        String chromosome = reads.get(0).getContig();

        int[] readOffsets = new int[readCount];
        boolean[] isFirstInPair = new boolean[readCount];
        boolean isDualStrand = mSequencingType == ILLUMINA && isDualStrandAndIsFirstInPair(reads, isFirstInPair);

        for(int i = 0; i < readCount; ++i)
        {
            readOffsets[i] = reads.get(i).getReadBases().length - baseLength;
        }

        int[] readPositionStartOffsets = setReadPositionStartOffsets(reads, consensusState.UnclippedPosStart, true);

        byte[] locationBases = new byte[readCount];
        byte[] locationQuals = new byte[readCount];

        for(int baseIndex = 0; baseIndex < baseLength; ++baseIndex)
        {
            // check bases at this index - work on the premise that most bases will agree
            boolean hasMismatch = false;
            byte maxQual = -1;
            int uncertainQualCount = 0;

            byte firstBase = NO_BASE;

            for(int r = 0; r < readCount; ++r)
            {
                // on reverse strand, say base length = 10 (so 0-9 for longest read), if a read has length 8 then it will
                locationBases[r] = NO_BASE;

                int readIndex = baseIndex - readPositionStartOffsets[r];

                // if a read has no position offset then its read bases can be read with the same index as the consensus (template) read
                // if its read position offset is -ve, then it starts earlier and so its read index will be ahead to match up

                if(readIndex < 0 || readIndex >= reads.get(r).getReadBases().length)
                    continue;

                locationBases[r] = reads.get(r).getReadBases()[readIndex];
                locationQuals[r] = reads.get(r).getBaseQualities()[readIndex];

                if(firstBase == NO_BASE)
                    firstBase = locationBases[r];
                else
                    hasMismatch |= locationBases[r] != firstBase;

                maxQual = maxQual(locationQuals[r], maxQual);

                if(isUncertainBaseQual(locationQuals[r]))
                    ++uncertainQualCount;
            }

            if(!hasMismatch)
            {
                consensusState.Bases[baseIndex] = firstBase;

                if(uncertainQualCount >= 0.5 * locationQuals.length)
                    consensusState.BaseQualities[baseIndex] = BASE_QUAL_MINIMUM;
                else
                    consensusState.BaseQualities[baseIndex] = maxQual;
            }
            else
            {
                int basePosition = consensusState.UnclippedPosStart + baseIndex; // works since there are no indels in this routine

                if(basePosition < 1 || basePosition > chromosomeLength)
                    basePosition = INVALID_POSITION; // protect against over-runs from soft-clips - rare but possible

                BaseQualPair consensusBaseAndQual = determineBaseAndQual(
                        locationBases, locationQuals, chromosome, basePosition, isDualStrand, isFirstInPair);

                consensusState.Bases[baseIndex] = consensusBaseAndQual.Base;
                consensusState.BaseQualities[baseIndex] = consensusBaseAndQual.Qual;
            }
        }
    }

    public BaseQualPair determineBaseAndQual(
            final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position,
            boolean isDualStrand, @Nullable final boolean[] isFirstInPair)
    {
        if(mSequencingType == ILLUMINA)
        {
            if(isDualStrand && position != INVALID_POSITION)
            {
                // split the reads into 2 consensus reads and then compare
                return IlluminaRoutines.determineDualStrandBaseAndQual(
                        isFirstInPair, locationBases, locationQuals, chromosome, position, mRefGenome);
            }
            else
            {
                return IlluminaRoutines.determineBaseAndQual(locationBases, locationQuals, chromosome, position, mRefGenome);
            }
        }
        else if(mSequencingType == SBX)
        {
            return SbxRoutines.determineBaseAndQual(locationBases, locationQuals, chromosome, position, mRefGenome);
        }
        else
        {
            // TODO: decide for Ultima
            return BaseQualPair.INVALID;
        }
    }
}
