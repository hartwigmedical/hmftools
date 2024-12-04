package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;

import java.util.List;
import com.hartwig.hmftools.common.sequencing.SequencingType;

import htsjdk.samtools.SAMRecord;

public class BaseBuilder
{
    private final RefGenome mRefGenome;
    private final BaseBuilderConfig mBaseBuilderConfig;

    // cached for the majority of successive reads being on the same chromosome, to protect against ref genome base requests beyond limits
    private int mChromosomeLength;

    public BaseBuilder(final RefGenome refGenome, final SequencingType sequencingType, final ConsensusStatistics consensusStats)
    {
        mRefGenome = refGenome;
        mBaseBuilderConfig = BaseBuilderConfig.fromSequencingType(sequencingType, refGenome, consensusStats);
        mChromosomeLength = 0;
    }

    public void setChromosomLength(int chromosomeLength) { mChromosomeLength = chromosomeLength; }
    public int chromosomeLength() { return mChromosomeLength; }
    public RefGenome refGenome() { return mRefGenome; }

    public BaseBuilderConfig baseBuilderConfig() { return mBaseBuilderConfig; }

    public static final byte NO_BASE = 0;
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
        boolean[] isFirstInPair = null;
        boolean isDualStrand = false;
        if(mBaseBuilderConfig.PairedReads)
        {
            isFirstInPair = new boolean[readCount];
            isDualStrand = isDualStrandAndIsFirstInPair(reads, isFirstInPair);
        }

        for(int i = 0; i < readCount; ++i)
        {
            readOffsets[i] = reads.get(i).getReadBases().length - baseLength;
        }

        byte[] locationBases = new byte[readCount];
        byte[] locationQuals = new byte[readCount];

        for(int baseIndex = 0; baseIndex < baseLength; ++baseIndex)
        {
            // check bases at this index - work on the premise that most bases will agree
            boolean hasMismatch = false;
            int maxQual = 0;
            byte firstBase = NO_BASE;

            for(int r = 0; r < readCount; ++r)
            {
                // on reverse strand, say base length = 10 (so 0-9 for longest read), if a read has length 8 then it will
                SAMRecord read = reads.get(r);

                locationBases[r] = NO_BASE;

                int readIndex;
                if(consensusState.IsForward)
                {
                    readIndex = baseIndex;

                    if(readOffsets[r] != 0 && baseIndex >= read.getReadBases().length)
                        continue;
                }
                else
                {
                    readIndex = baseIndex + readOffsets[r];

                    if(readIndex < 0)
                        continue;
                }

                locationBases[r] = reads.get(r).getReadBases()[readIndex];
                locationQuals[r] = reads.get(r).getBaseQualities()[readIndex];

                if(firstBase == NO_BASE)
                    firstBase = locationBases[r];
                else
                    hasMismatch |= locationBases[r] != firstBase;

                maxQual = max(locationQuals[r], maxQual);
            }

            if(!hasMismatch && mBaseBuilderConfig.UseSimpleNoMismatchLogic)
            {
                consensusState.Bases[baseIndex] = firstBase;
                consensusState.BaseQualities[baseIndex] = (byte)maxQual;
            }
            else
            {
                int basePosition = consensusState.MinUnclippedPosStart + baseIndex;

                if(basePosition < 1 || basePosition > chromosomeLength)
                    basePosition = INVALID_POSITION; // protect against over-runs from soft-clips - rare but possible

                byte[] consensusBaseAndQual = mBaseBuilderConfig.determineBaseAndQual(
                        isDualStrand, isFirstInPair, locationBases, locationQuals, chromosome, basePosition);
                consensusState.Bases[baseIndex] = consensusBaseAndQual[0];
                consensusState.BaseQualities[baseIndex] = consensusBaseAndQual[1];
            }
        }
    }

    public static boolean isDualStrandAndIsFirstInPair(final List<SAMRecord> reads, final boolean[] isFirstInPairOut)
    {
        if(!reads.get(0).getReadPairedFlag())
            return false;

        boolean isDualStrand = false;

        for(int i = 0; i < reads.size(); ++i)
        {
            isFirstInPairOut[i] = reads.get(i).getFirstOfPairFlag();

            if(i > 0)
                isDualStrand |= isFirstInPairOut[0] != isFirstInPairOut[i];
        }

        return isDualStrand;
    }
}
