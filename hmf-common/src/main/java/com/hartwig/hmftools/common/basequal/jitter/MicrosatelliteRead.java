package com.hartwig.hmftools.common.basequal.jitter;

import static com.hartwig.hmftools.common.basequal.jitter.JitterAnalyserConstants.MIN_FLANKING_BASE_MATCHES;

import org.apache.commons.lang3.Validate;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

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

    final RefGenomeMicrosatellite refGenomeMicrosatellite;
    boolean shouldDropRead = false;
    int numAligned = 0;
    int numInserted = 0;
    int numDeleted = 0;

    int numMatchedBefore = 0;
    int numMatchedAfter = 0;

    private MicrosatelliteRead(final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record)
    {
        this.refGenomeMicrosatellite = refGenomeMicrosatellite;

        // this read needs to wholly contain the homopolymer to be counted
        if(record.getAlignmentStart() > refGenomeMicrosatellite.referenceStart() ||
        record.getAlignmentEnd() < refGenomeMicrosatellite.referenceEnd())
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
            sLogger.trace("read: {} not enough flanking bases aligned(before: {}, after: {}), dropping read",
                    record, numMatchedBefore, numMatchedAfter);
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

    public int readRepeatLength()
    {
        Validate.isTrue(!shouldDropRead);
        // calculate the repeat length
        int readRepeatLength = numAligned + numInserted;

        // check that deleted makes sense
        // Validate.isTrue(readRepeatLength == refGenomeHomopolymer.numRepeat - numDeleted + numInserted);

        return readRepeatLength;
    }

    public int numRepeatUnits()
    {
        return readRepeatLength() / refGenomeMicrosatellite.unit.length;
    }

    public static MicrosatelliteRead from(final RefGenomeMicrosatellite refGenomeMicrosatellite, final SAMRecord record)
    {
        return new MicrosatelliteRead(refGenomeMicrosatellite, record);
    }

    public void handleAlignment(final SAMRecord record, final CigarElement e, final int startReadIndex, final int startRefPos)
    {
        // check if this alignment spans the repeat
        // TODO: check for substitution??
        int numAlignedToMs = Math.min(startRefPos + e.getLength(), refGenomeMicrosatellite.referenceEnd() + 1) -
                             Math.max(startRefPos, refGenomeMicrosatellite.referenceStart());

        if(numAlignedToMs > 0)
        {
            // this is in the repeat section
            numAligned += numAlignedToMs;
        }

        int numAlignedToFlankingStart = Math.min(startRefPos + e.getLength(), refGenomeMicrosatellite.referenceStart()) -
                                        Math.max(startRefPos, refGenomeMicrosatellite.referenceStart() - MIN_FLANKING_BASE_MATCHES);

        if(numAlignedToFlankingStart > 0)
        {
            numMatchedBefore += numAlignedToFlankingStart;
        }

        int numAlignedToFlankingEnd = Math.min(startRefPos + e.getLength(), refGenomeMicrosatellite.referenceEnd() + 1 + MIN_FLANKING_BASE_MATCHES) -
                Math.max(startRefPos, refGenomeMicrosatellite.referenceEnd() + 1);

        if(numAlignedToFlankingEnd > 0)
        {
            numMatchedAfter += numAlignedToFlankingEnd;
        }
    }

    public void handleInsert(final SAMRecord record, final CigarElement e, final int readIndex, final int refPos)
    {
        if(refPos >= refGenomeMicrosatellite.referenceStart() && refPos <= refGenomeMicrosatellite.referenceEnd() + 1)
        {
            sLogger.trace("read: {} inserted {} bases", record, e.getLength());

            // check the repeat unit
            for(int i = 0; i < e.getLength(); ++i)
            {
                // the inserted bases must be a multiple of the repeat unit
                if(record.getReadBases()[readIndex + i] != refGenomeMicrosatellite.unit[ i % refGenomeMicrosatellite.unit.length ])
                {
                    shouldDropRead = true;
                    // should drop this base
                    sLogger.trace("read: {} inserted bases: {} vs ms unit: {} mismatch, dropping read",
                            record, record.getReadString().substring(readIndex, readIndex + e.getLength()),
                            refGenomeMicrosatellite.unitString());
                    break;
                }
            }

            numInserted += e.getLength();
        }
        else if(refPos + MIN_FLANKING_BASE_MATCHES >= refGenomeMicrosatellite.referenceStart() &&
                refPos - MIN_FLANKING_BASE_MATCHES <= refGenomeMicrosatellite.referenceEnd() + 1)
        {
            // there is an insert very close to the homopolymer, drop this read
            shouldDropRead = true;
        }
    }

    public void handleDelete(final SAMRecord record, final CigarElement e, final int readIndex, final int startRefPos)
    {
        int endRefPos = startRefPos + e.getLength() - 1;
        if(startRefPos >= refGenomeMicrosatellite.referenceStart() && endRefPos <= refGenomeMicrosatellite.referenceEnd())
        {
            // the whole delete is inside the repeat, this is nice simple case
            numDeleted += e.getLength();
        }
        else if((startRefPos < refGenomeMicrosatellite.referenceStart() && endRefPos >= refGenomeMicrosatellite.referenceStart() - 1) ||
                (startRefPos <= refGenomeMicrosatellite.referenceEnd() + 1 && endRefPos > refGenomeMicrosatellite.referenceEnd()))
        {
            // if the delete cross over the start or the end, drop the read
            // also drop if the delete is just before the start of the polymer
            shouldDropRead = true;
        }
    }

    public void handleLeftSoftClip(final SAMRecord record, final CigarElement element)
    {
        // drop this read completely if the soft clip is near the repeat
        if(record.getAlignmentStart() >= refGenomeMicrosatellite.referenceStart())
        {
            shouldDropRead = true;
        }
    }

    public void handleRightSoftClip(final SAMRecord record, final CigarElement element, int startReadIndex, int startRefPosition)
    {
        // drop this read completely if the soft clip is inside the repeat
        if(startRefPosition <= refGenomeMicrosatellite.referenceEnd())
        {
            shouldDropRead = true;
        }
    }

    void traverseCigar(final SAMRecord record)
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
                        handleLeftSoftClip(record, e);
                    }
                    else if(i == cigar.numCigarElements() - 1)
                    {
                        handleRightSoftClip(record, e, readIndex, refBase);
                    }
                    readIndex += e.getLength();
                    break; // soft clip read bases
                case N:
                    //handleSkippedReference(record, e, readIndex, refBase);
                    refBase += e.getLength();
                    break;  // reference skip
                case D:
                    handleDelete(record, e, readIndex, refBase);
                    refBase += e.getLength();
                    break;
                case I:
                    handleInsert(record, e, readIndex, refBase);
                    readIndex += e.getLength();
                    break;
                case M:
                case EQ:
                case X:
                    handleAlignment(record, e, readIndex, refBase);
                    readIndex += e.getLength();
                    refBase += e.getLength();
                    break;
                default:
                    throw new IllegalStateException("Case statement didn't deal with op: " + e.getOperator() + "in CIGAR: " + cigar);
            }
        }
    }
}
