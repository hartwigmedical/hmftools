package com.hartwig.hmftools.redux.common;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.addConsensusReadAttribute;
import static com.hartwig.hmftools.common.bam.UmiReadType.DUAL;
import static com.hartwig.hmftools.common.bam.UmiReadType.SINGLE;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.ReadInfo.readToString;

import java.util.List;
import java.util.function.BiConsumer;

import com.hartwig.hmftools.common.bam.UmiReadType;
import com.hartwig.hmftools.redux.consensus.ConsensusReadInfo;
import com.hartwig.hmftools.redux.consensus.ConsensusReads;

import htsjdk.samtools.SAMRecord;

public abstract class DuplicateGroup
{
    protected final String mUmiId; // the UMI if enabled

    protected SAMRecord mConsensusRead;
    private SAMRecord mPrimaryRead; // if no consensus is formed, the selected primary read
    private boolean mDualStrand;

    public DuplicateGroup(final String id)
    {
        mUmiId = id;

        mConsensusRead = null;
        mPrimaryRead = null;
        mDualStrand = false;
    }

    public abstract int readCount();
    public abstract List<SAMRecord> reads();

    public abstract boolean readIsLower();

    public SAMRecord consensusRead()
    {
        return mConsensusRead;
    }

    public void setPrimaryRead(final SAMRecord read)
    {
        mPrimaryRead = read;
    }

    public boolean isPrimaryRead(final SAMRecord read)
    {
        return mPrimaryRead == read;
    }

    public String umiId()
    {
        return mUmiId;
    }

    public void registerDualStrand()
    {
        mDualStrand = true;
    }

    public boolean hasDualStrand()
    {
        return mDualStrand;
    }

    public abstract String consensusCoordinatesKey();

    public void formConsensusRead(final ConsensusReads consensusReads)
    {
        List<SAMRecord> reads = reads();
        boolean readIsLower = readIsLower();

        try
        {
            ConsensusReadInfo consensusReadInfo = consensusReads.createConsensusRead(reads, readIsLower, mUmiId);

            // set consensus read attributes
            int firstInPairCount = (int) reads.stream().filter(x -> !x.getReadPairedFlag() || x.getFirstOfPairFlag()).count();
            int readCount = reads.size();
            boolean isPrimaryGroup = firstInPairCount >= readCount / 2;

            if(!isPrimaryGroup)
            {
                firstInPairCount = readCount - firstInPairCount; // adjusted so both reads report the same ratio
            }

            UmiReadType umiReadType = mDualStrand ? DUAL : SINGLE;

            addConsensusReadAttribute(consensusReadInfo.ConsensusRead, readCount, firstInPairCount, umiReadType);

            mConsensusRead = consensusReadInfo.ConsensusRead;
        }
        catch(Exception e)
        {
            RD_LOGGER.error("error forming consensus: {}", toString());

            for(SAMRecord read : reads)
            {
                RD_LOGGER.error("read: {}", readToString(read));
            }

            e.printStackTrace();
        }
    }

    public abstract void processReads(final BiConsumer<SAMRecord, FragmentCoords> consumer);
}
