package com.hartwig.hmftools.redux.jitter;

import static com.hartwig.hmftools.common.bam.ConsensusType.NONE;
import static com.hartwig.hmftools.redux.jitter.JitterAnalyserConstants.MIN_FLANKING_BASE_MATCHES;

import com.hartwig.hmftools.common.bam.ConsensusType;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

// we need to look through the cigar to decide whether this read has matches
// we can either
// 1. check the M elements and see if they add up
// 2. check the D elements to see if any was deleted
// 3. we must check the I elements in case the polymer was lengthened.
class MicrosatelliteRead
{
    public static final Logger sLogger = LogManager.getLogger(MicrosatelliteRead.class);

    private static final ThreadLocal<MicrosatelliteRead> THREAD_INSTANCE = new ThreadLocal<>()
    {
        @Override
        protected MicrosatelliteRead initialValue()
        {
            return new MicrosatelliteRead();
        }
    };

    private RefGenomeMicrosatellite mRefGenomeMicrosatellite;
    private ConsensusType mConsensusType;

    public boolean shouldDropRead;

    private int numAligned;
    private int numInserted;
    private int numDeleted;

    private int numMatchedBefore;
    private int numMatchedAfter;

    private MicrosatelliteRead()
    {
        clear();
    }

    private void clear()
    {
        mRefGenomeMicrosatellite = null;
        mConsensusType = null;

        shouldDropRead = false;

        numAligned = 0;
        numInserted = 0;
        numDeleted = 0;

        numMatchedBefore = 0;
        numMatchedAfter = 0;
    }

    private void analyse(final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record, @Nullable final ConsensusMarker consensusMarker)
    {
        clear();

        mRefGenomeMicrosatellite = refGenomeMicrosatellite;

        // this read needs to wholly contain the homopolymer to be counted
        if(record.getAlignmentStart() > refGenomeMicrosatellite.referenceStart()
        || record.getAlignmentEnd() < refGenomeMicrosatellite.referenceEnd())
        {
            shouldDropRead = true;
        }
        else
        {
            traverseCigar(record);
        }

        // check that we got the flanking bases aligned
        if(numMatchedBefore < MIN_FLANKING_BASE_MATCHES || numMatchedAfter < MIN_FLANKING_BASE_MATCHES)
        {
            shouldDropRead = true;
            // sLogger.trace("read: {} not enough flanking bases aligned(before: {}, after: {}), dropping read",
            //         record, numMatchedBefore, numMatchedAfter);
        }

        if(!shouldDropRead)
        {
            if(consensusMarker == null)
                mConsensusType = NONE;
            else
                mConsensusType = consensusMarker.consensusType(refGenomeMicrosatellite, record);
        }
        else
        {
            mConsensusType = NONE;
        }

        // do some validation and logging
        if(!shouldDropRead)
        {
            int readRepeatLength = numAligned + numInserted;
            if(readRepeatLength != refGenomeMicrosatellite.baseLength() - numDeleted + numInserted)
            {
                sLogger.error("read({}) {}, incorrect read repeat length({}) numAligned({}) numInserted({}) numDeleted({}) " +
                                "homopolymer({})",
                        record, record.getCigarString(), readRepeatLength, numAligned, numInserted, numDeleted, refGenomeMicrosatellite.genomeRegion);
            }
        }
    }

    public ConsensusType consensusType() { return mConsensusType; }

    public int readRepeatLength()
    {
        Validate.isTrue(!shouldDropRead);
        // calculate the repeat length
        int readRepeatLength = numAligned + numInserted;

        return readRepeatLength;
    }

    public int numRepeatUnits()
    {
        return readRepeatLength() / mRefGenomeMicrosatellite.unit.length;
    }

    public int jitter() { return numRepeatUnits() - mRefGenomeMicrosatellite.numRepeat; }

    @Nullable
    public static MicrosatelliteRead from(
            final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record, @Nullable final ConsensusMarker consensusMarker)
    {
        MicrosatelliteRead instance = THREAD_INSTANCE.get();
        instance.analyse(refGenomeMicrosatellite, record, consensusMarker);
        if(instance.consensusType() == null)
            return null;

        return instance;
    }

    private void handleAlignment(final CigarElement e, final int startRefPos)
    {
        // check if this alignment spans the repeat
        // TODO: check for substitution??
        int numAlignedToMs = Math.min(startRefPos + e.getLength(), mRefGenomeMicrosatellite.referenceEnd() + 1) -
                Math.max(startRefPos, mRefGenomeMicrosatellite.referenceStart());

        if(numAlignedToMs > 0)
        {
            // this is in the repeat section
            numAligned += numAlignedToMs;
        }

        int numAlignedToFlankingStart = Math.min(startRefPos + e.getLength(), mRefGenomeMicrosatellite.referenceStart()) -
                Math.max(startRefPos, mRefGenomeMicrosatellite.referenceStart() - MIN_FLANKING_BASE_MATCHES);

        if(numAlignedToFlankingStart > 0)
        {
            numMatchedBefore += numAlignedToFlankingStart;
        }

        int numAlignedToFlankingEnd =
                Math.min(startRefPos + e.getLength(), mRefGenomeMicrosatellite.referenceEnd() + 1 + MIN_FLANKING_BASE_MATCHES) -
                        Math.max(startRefPos, mRefGenomeMicrosatellite.referenceEnd() + 1);

        if(numAlignedToFlankingEnd > 0)
        {
            numMatchedAfter += numAlignedToFlankingEnd;
        }
    }

    private void handleInsert(final SAMRecord record, final CigarElement e, final int readIndex, final int refPos)
    {
        if(refPos >= mRefGenomeMicrosatellite.referenceStart() && refPos <= mRefGenomeMicrosatellite.referenceEnd() + 1)
        {
            sLogger.trace("read: {} inserted {} bases", record, e.getLength());

            // check the repeat unit
            for(int i = 0; i < e.getLength(); ++i)
            {
                // the inserted bases must be a multiple of the repeat unit
                if(record.getReadBases()[readIndex + i] != mRefGenomeMicrosatellite.unit[i % mRefGenomeMicrosatellite.unit.length])
                {
                    shouldDropRead = true;
                    // should drop this base
                    sLogger.trace("read: {} inserted bases: {} vs ms unit: {} mismatch, dropping read",
                            record, record.getReadString().substring(readIndex, readIndex + e.getLength()),
                            mRefGenomeMicrosatellite.unitString());
                    break;
                }
            }

            numInserted += e.getLength();
        }
        else if(refPos + MIN_FLANKING_BASE_MATCHES >= mRefGenomeMicrosatellite.referenceStart() &&
                refPos - MIN_FLANKING_BASE_MATCHES <= mRefGenomeMicrosatellite.referenceEnd() + 1)
        {
            // there is an insert very close to the homopolymer, drop this read
            shouldDropRead = true;
        }
    }

    private void handleDelete(final CigarElement e, final int startRefPos)
    {
        int endRefPos = startRefPos + e.getLength() - 1;
        if(startRefPos >= mRefGenomeMicrosatellite.referenceStart() && endRefPos <= mRefGenomeMicrosatellite.referenceEnd())
        {
            // the whole delete is inside the repeat, this is nice simple case
            numDeleted += e.getLength();
        }
        else if((startRefPos < mRefGenomeMicrosatellite.referenceStart() && endRefPos >= mRefGenomeMicrosatellite.referenceStart() - 1) ||
                (startRefPos <= mRefGenomeMicrosatellite.referenceEnd() + 1 && endRefPos > mRefGenomeMicrosatellite.referenceEnd()))
        {
            // if the delete cross over the start or the end, drop the read
            // also drop if the delete is just before the start of the polymer
            shouldDropRead = true;
        }
    }

    private void handleLeftSoftClip(final SAMRecord record)
    {
        // drop this read completely if the soft clip is near the repeat
        if(record.getAlignmentStart() >= mRefGenomeMicrosatellite.referenceStart())
        {
            shouldDropRead = true;
        }
    }

    private void handleRightSoftClip(int startRefPosition)
    {
        // drop this read completely if the soft clip is inside the repeat
        if(startRefPosition <= mRefGenomeMicrosatellite.referenceEnd())
        {
            shouldDropRead = true;
        }
    }

    private void traverseCigar(final SAMRecord record)
    {
        final Cigar cigar = record.getCigar();

        int readIndex = 0;
        int refBase = record.getAlignmentStart();

        for(int i = 0; i < cigar.numCigarElements(); i++)
        {
            final CigarElement e = cigar.getCigarElement(i);
            switch(e.getOperator())
            {
                case H:
                    break; // ignore hard clips - no need to skip either bases or positions
                case P:
                    break; // ignore pads
                case S:
                    if(i == 0)
                    {
                        handleLeftSoftClip(record);
                    }
                    else if(i == cigar.numCigarElements() - 1)
                    {
                        handleRightSoftClip(refBase);
                    }
                    readIndex += e.getLength();
                    break; // soft clip read bases
                case N:
                    //handleSkippedReference(record, e, readIndex, refBase);
                    refBase += e.getLength();
                    break;  // reference skip
                case D:
                    handleDelete(e, refBase);
                    refBase += e.getLength();
                    break;
                case I:
                    handleInsert(record, e, readIndex, refBase);
                    readIndex += e.getLength();
                    break;
                case M:
                case EQ:
                case X:
                    handleAlignment(e, refBase);
                    readIndex += e.getLength();
                    refBase += e.getLength();
                    break;
                default:
                    throw new IllegalStateException("Case statement didn't deal with op: " + e.getOperator() + " in CIGAR: " + cigar);
            }
        }
    }
}
