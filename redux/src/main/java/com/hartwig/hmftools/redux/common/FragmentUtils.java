package com.hartwig.hmftools.redux.common;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.getFivePrimeUnclippedPosition;
import static com.hartwig.hmftools.redux.common.FragmentCoordinates.formCoordinate;
import static com.hartwig.hmftools.redux.common.FragmentCoordinates.formKey;
import static com.hartwig.hmftools.redux.common.FragmentStatus.DUPLICATE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.NONE;
import static com.hartwig.hmftools.redux.common.FragmentStatus.CANDIDATE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MATE_CIGAR_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.orientation;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMRecord;

public class FragmentUtils
{
    public static FragmentCoordinates getFragmentCoordinates(final List<SAMRecord> reads, final boolean useMateCigar)
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
                getFivePrimeUnclippedPosition(firstRead) :
                getFivePrimeUnclippedPosition(firstRead.getAlignmentStart(), firstRead.getCigarString(), readForwardStrand);

        int readStrandPosition = readForwardStrand ? readCoordinate : -readCoordinate;
        String readCoordStr = formCoordinate(firstRead.getReferenceName(), readCoordinate, readForwardStrand);

        if(!firstRead.getReadPairedFlag() || firstRead.getReadUnmappedFlag() || firstRead.getMateUnmappedFlag())
        {
            // include the fragment length
            return new FragmentCoordinates(formKey(readCoordStr, firstRead.getInferredInsertSize()), readStrandPosition, true);
        }

        if(mateRead == null)
        {
            if(!useMateCigar || !firstRead.hasAttribute(MATE_CIGAR_ATTRIBUTE))
            {
                // the fragment orientation will  be accurately set once both reads are collated
                return new FragmentCoordinates(readCoordStr, readStrandPosition, firstRead.getFirstOfPairFlag(), true);
            }
        }

        boolean mateForwardStrand = !firstRead.getMateNegativeStrandFlag();

        int mateCoordinate;

        if(mateRead != null)
        {
            mateRead.getReferenceName();
            mateCoordinate = mateRead.getCigar() != null ?
                    getFivePrimeUnclippedPosition(mateRead) :
                    getFivePrimeUnclippedPosition(mateRead.getAlignmentStart(), mateRead.getCigarString(), mateForwardStrand);
        }
        else
        {
            String mateCigar = firstRead.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
            mateCoordinate = getFivePrimeUnclippedPosition(firstRead.getMateAlignmentStart(), mateCigar, mateForwardStrand);
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

        boolean lowerReadFirst = readLowerPos ? firstRead.getFirstOfPairFlag() : !firstRead.getFirstOfPairFlag();

        return readLowerPos ?
                new FragmentCoordinates(formKey(readCoordStr, mateCoordStr), readStrandPosition, lowerReadFirst)
                : new FragmentCoordinates(formKey(mateCoordStr, readCoordStr), mateStrandPosition, lowerReadFirst);
    }

    public static FragmentStatus calcFragmentStatus(final Fragment first, final Fragment second, boolean requireOrientationMatch)
    {
        if(first.unpaired() != second.unpaired())
            return NONE;

        if(!first.coordinates().Incomplete && !second.coordinates().Incomplete)
            return first.coordinates().matches(second.coordinates(), requireOrientationMatch) ? DUPLICATE : NONE;

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
}