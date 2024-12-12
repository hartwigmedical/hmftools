package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.sequencing.SequencingType.BIOMODAL;
import static com.hartwig.hmftools.common.sequencing.SequencingType.SBX;

import com.hartwig.hmftools.common.sequencing.SequencingType;

import htsjdk.samtools.SAMRecord;

public abstract class BaseBuilderConfig
{
    public final boolean PairedReads;
    public final boolean UseSimpleNoMismatchLogic;

    protected final RefGenome mRefGenome;
    protected final ConsensusStatistics mConsensusStats;

    public BaseBuilderConfig(final RefGenome refGenome, boolean pairedReads, boolean useSimpleNoMismatchLogic,
            final ConsensusStatistics consensusStats)
    {
        mRefGenome = refGenome;
        PairedReads = pairedReads;
        UseSimpleNoMismatchLogic = useSimpleNoMismatchLogic;
        mConsensusStats = consensusStats;
    }

    public byte[] preprocessReadBases(final SAMRecord read)
    {
        return read.getReadBases();
    }

    public abstract byte[] determineBaseAndQual(boolean isDualStrand, final boolean[] isFirstInPair, final byte[] locationBases,
            final byte[] locationQuals, final String chromosome, int position);

    public static BaseBuilderConfig fromSequencingType(final SequencingType sequencingType, final RefGenome refGenome,
            final ConsensusStatistics consensusStats)
    {
        if(sequencingType == SBX)
            return new SBXBaseBuilderConfig(refGenome, consensusStats);

        if(sequencingType == BIOMODAL)
            return new BiomodalBaseBuilderConfig(refGenome, consensusStats);

        return new DefaultBaseBuilderConfig(refGenome, consensusStats);
    }
}
