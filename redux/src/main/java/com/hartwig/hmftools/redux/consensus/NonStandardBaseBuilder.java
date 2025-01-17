package com.hartwig.hmftools.redux.consensus;

import java.util.List;

import com.hartwig.hmftools.common.sequencing.SequencingType;

import htsjdk.samtools.SAMRecord;

public abstract class NonStandardBaseBuilder
{
    protected final RefGenome mRefGenome;
    protected int mChromosomeLength;

    public NonStandardBaseBuilder(final RefGenome refGenome)
    {
        mRefGenome = refGenome;
        mChromosomeLength = 0;
    }

    public abstract void buildConsensusRead(final List<SAMRecord> reads, final ConsensusState consensusState, boolean hasIndels);

    public void setChromosomeLength(int chromosomeLength) { mChromosomeLength = chromosomeLength; }

    public static NonStandardBaseBuilder fromSequencingType(final SequencingType sequencingType, final RefGenome refGenome)
    {
        return null;
    }
}
