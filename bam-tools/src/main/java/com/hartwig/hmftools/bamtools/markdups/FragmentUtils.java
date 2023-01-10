package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.BmConfig.BM_LOGGER;
import static com.hartwig.hmftools.bamtools.markdups.FragmentCoordinates.NO_COORDS;
import static com.hartwig.hmftools.bamtools.markdups.FragmentCoordinates.formCoordinate;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.NONE;
import static com.hartwig.hmftools.bamtools.markdups.FragmentStatus.CANDIDATE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.samtools.SamRecordUtils.orientation;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.List;

import htsjdk.samtools.SAMRecord;

public class FragmentUtils
{
    public static int getUnclippedPosition(final SAMRecord read)
    {
        int position;

        if(orientation(read) == POS_ORIENT)
        {
            position = read.getAlignmentStart();
            if(read.getCigar().isLeftClipped())
                position -= read.getCigar().getFirstCigarElement().getLength();
        }
        else
        {
            position = read.getAlignmentEnd();
            if(read.getCigar().isRightClipped())
                position += read.getCigar().getLastCigarElement().getLength();
        }

        return position;
    }

    public static FragmentCoordinates getFragmentCoordinates(final List<SAMRecord> reads)
    {
        SAMRecord firstRead = null;
        SAMRecord mateRead = null;

        for(SAMRecord read : reads)
        {
            if(read.getSupplementaryAlignmentFlag())
                continue;

            if(firstRead == null)
            {
                firstRead = read;
            }
            else
            {
                mateRead = read;
                break;
            }
        }

        if(firstRead.getReadUnmappedFlag() && mateRead != null)
        {
            firstRead = mateRead;
            mateRead = null;
        }

        boolean readForwardStrand = orientation(firstRead) == POS_ORIENT;

        int readCoordinate = firstRead.getCigar() != null ?
                getUnclippedPosition(firstRead) : getUnclippedPosition(firstRead.getAlignmentStart(), firstRead.getCigarString(), readForwardStrand);

        int readStrandPosition = readForwardStrand ? readCoordinate : -readCoordinate;
        String readCoordStr = formCoordinate(firstRead.getReferenceName(), readCoordinate, readForwardStrand);

        if(!firstRead.getReadPairedFlag() || firstRead.getReadUnmappedFlag() || firstRead.getMateUnmappedFlag())
            return new FragmentCoordinates(readCoordStr, readStrandPosition);

        if(mateRead == null && !firstRead.hasAttribute(MATE_CIGAR_ATTRIBUTE))
            return new FragmentCoordinates(readCoordStr, readStrandPosition, true);

        boolean mateForwardStrand = !firstRead.getMateNegativeStrandFlag();
        int mateCoordinate;

        if(mateRead != null)
        {
            mateRead.getReferenceName();
            mateCoordinate = mateRead.getCigar() != null ?
                    getUnclippedPosition(mateRead) : getUnclippedPosition(mateRead.getAlignmentStart(), mateRead.getCigarString(), mateForwardStrand);
        }
        else
        {
            String mateCigar = firstRead.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
            mateCoordinate = getUnclippedPosition(firstRead.getMateAlignmentStart(), mateCigar, mateForwardStrand);
        }

        int mateStrandPosition = mateForwardStrand ? mateCoordinate : -mateCoordinate;
        String mateCoordStr = formCoordinate(firstRead.getMateReferenceName(), mateCoordinate, mateForwardStrand);

        boolean readLowerPos;
        if(firstRead.getReferenceIndex() == firstRead.getMateReferenceIndex())
        {
            readLowerPos = readCoordinate <= mateCoordinate;
        }
        else
        {
            readLowerPos = firstRead.getReferenceIndex() < firstRead.getMateReferenceIndex();
        }

        return readLowerPos ?
                new FragmentCoordinates(readCoordStr + "_" + mateCoordStr, readStrandPosition)
                : new FragmentCoordinates(mateCoordStr + "_" + readCoordStr, mateStrandPosition);
    }

    public static FragmentStatus calcFragmentStatus(final Fragment first, final Fragment second)
    {
        if(first.unpaired() != second.unpaired())
            return NONE;

        if(!first.coordinates().Incomplete && !second.coordinates().Incomplete)
            return first.coordinates().Key.equals(second.coordinates().Key) ? DUPLICATE : NONE;

        if(first.initialPosition() != second.initialPosition())
            return NONE;

        // mate start positions must be within close proximity
        SAMRecord firstRead = first.reads().get(0);
        SAMRecord secondRead = second.reads().get(0);

        if(!firstRead.getMateReferenceName().equals(secondRead.getMateReferenceName()))
            return NONE;

        if(firstRead.getMateNegativeStrandFlag() != secondRead.getMateNegativeStrandFlag())
            return NONE;

        return abs(firstRead.getMateAlignmentStart() - secondRead.getMateAlignmentStart()) < firstRead.getReadLength()
                ? CANDIDATE : NONE;
    }

    private static final String CHR_PARTITION_DELIM = "_";

    public static String formChromosomePartition(final String chromosome, int position, int partitionSize)
    {
        int partition = position / partitionSize;
        return chromosomeIndicator(chromosome) + partition;
    }

    public static String chromosomeIndicator(final String chromosome)
    {
        return chromosome + CHR_PARTITION_DELIM;
    }

    public static String readToString(final SAMRecord read)
    {
        return format("id(%s) coords(%s:%d-%d) cigar(%s) mate(%s:%d) flags(%d)",
                read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(),
                read.getCigarString(), read.getMateReferenceName(), read.getMateAlignmentStart(), read.getFlags());
    }

    public static int getUnclippedPosition(final int readStart, final String cigarStr, final boolean forwardStrand)
    {
        int currentPosition = readStart;
        int elementLength = 0;

        for(int i = 0; i < cigarStr.length(); ++i)
        {
            char c = cigarStr.charAt(i);
            boolean isAddItem = (c == 'D' || c == 'M' || c == 'S' || c == 'N');

            if(isAddItem)
            {
                if(forwardStrand)
                {
                    // back out the left clip if present
                    return c == 'S' ? readStart - elementLength : readStart;
                }

                if(c == 'S' && readStart == currentPosition)
                {
                    // ignore left-clip when getting reverse strand position
                }
                else
                {
                    currentPosition += elementLength;
                }

                elementLength = 0;
                continue;
            }

            int digit = c - '0';
            if (digit >= 0 && digit <= 9)
            {
                elementLength = elementLength * 10 + digit;
            }
            else
            {
                elementLength = 0;
            }
        }

        // always pointing to the start of the next element, so need to move back a base
        return currentPosition - 1;
    }
}