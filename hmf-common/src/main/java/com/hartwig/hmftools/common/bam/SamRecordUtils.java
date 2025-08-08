package com.hartwig.hmftools.common.bam;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.CigarUtils.getReadBoundaryPosition;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;

import static htsjdk.samtools.CigarOperator.D;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

public final class SamRecordUtils
{
    public static final String SUPPLEMENTARY_ATTRIBUTE = SAMTag.SA.name();
    public static final String MATE_CIGAR_ATTRIBUTE = SAMTag.MC.name();
    public static final String MATE_QUALITY_ATTRIBUTE = SAMTag.MQ.name();
    public static final String NUM_MUTATONS_ATTRIBUTE = SAMTag.NM.name();
    public static final String MISMATCHES_AND_DELETIONS_ATTRIBUTE = SAMTag.MD.name();
    public static final String SECONDARY_ATTRIBUTE = SAMTag.HI.name();
    public static final String ALIGNMENT_SCORE_ATTRIBUTE = SAMTag.AS.name();
    public static final String READ_GROUP_ATTRIBUTE = SAMTag.RG.name();

    public static final String XS_ATTRIBUTE = "XS";


    // Redux tags
    public static final String CONSENSUS_READ_ATTRIBUTE = "CR";
    public static final String CONSENSUS_TYPE_ATTRIBUTE = "UT"; // originally UMI type - single, dual or none
    public static final String UMI_ATTRIBUTE = "UI"; // the UMI group ID
    public static final String CONSENSUS_INFO_DELIM = ";";
    public static final String BASE_MODIFICATIONS_ATTRIBUTE = "MM";

    public static final String UNMAP_ATTRIBUTE = "UM"; // the read was unmapped
    public static final String UNMAPP_COORDS_DELIM = ":";

    public static final String NO_CHROMOSOME_NAME = "*";
    public static final String NO_CIGAR = "*";
    public static final int NO_CHROMOSOME_INDEX = -1;
    public static final int NO_POSITION = 0;
    public static final int INVALID_READ_INDEX = -1;

    public static final int UNSET_COUNT = -1;

    public static final int PHRED_OFFSET = 33;

    public static final Logger SAM_LOGGER = LogManager.getLogger(SamRecordUtils.class);

    // convenience methods to avoid triggering a crash if unpaired
    public static boolean firstInPair(final SAMRecord record)
    {
        return !record.getReadPairedFlag() || record.getFirstOfPairFlag();
    }

    public static boolean secondInPair(final SAMRecord record)
    {
        return record.getReadPairedFlag() && record.getSecondOfPairFlag();
    }

    public static boolean mateNegativeStrand(final SAMRecord record)
    {
        return record.getReadPairedFlag() && record.getMateNegativeStrandFlag();
    }

    public static int inferredInsertSize(final SAMRecord record)
    {
        if(record.getReadPairedFlag())
            return record.getInferredInsertSize();

        int insertSize = record.getInferredInsertSize();
        return insertSize > 0 ? insertSize : record.getReadBases().length;
    }

    public static int inferredInsertSizeAbs(final SAMRecord record) { return abs(inferredInsertSize(record)); }

    public static boolean properPair(final SAMRecord record)
    {
        return record.getReadPairedFlag() && record.getProperPairFlag();
    }

    public static boolean mateUnmapped(final SAMRecord record)
    {
        return record.getReadPairedFlag() && record.getMateUnmappedFlag();
    }

    public static byte orientation(final SAMRecord read) { return !read.getReadNegativeStrandFlag() ? ORIENT_FWD : ORIENT_REV; }

    public static boolean isFlagSet(final int flags, final SAMFlag flag) { return (flags & flag.intValue()) != 0; }

    public static int getBaseQuality(final char quality)
    {
        return quality - PHRED_OFFSET;
    }

    public static void addConsensusReadAttribute(
            final SAMRecord record, int readCount, int firstInPairCount, final ConsensusType consensusType, int pcrClusterCount)
    {
        String consensusReadValue = pcrClusterCount == UNSET_COUNT
                ? format("%d;%d", readCount, firstInPairCount)
                : format("%d;%d;%d", readCount, firstInPairCount, pcrClusterCount);

        record.setAttribute(CONSENSUS_READ_ATTRIBUTE, consensusReadValue);
        record.setAttribute(CONSENSUS_TYPE_ATTRIBUTE, consensusType.toString());
    }

    @VisibleForTesting
    public static void addConsensusReadAttribute(final SAMRecord record, int readCount, int firstInPairCount, final ConsensusType consensusType)
    {
        addConsensusReadAttribute(record, readCount, firstInPairCount, consensusType, UNSET_COUNT);
    }

    public static ConsensusType extractConsensusType(final SAMRecord record)
    {
        String consensusTypeStr = record.getStringAttribute(CONSENSUS_TYPE_ATTRIBUTE);

        if(consensusTypeStr != null)
            return ConsensusType.valueOf(consensusTypeStr);
        else
            return ConsensusType.NONE;
    }

    public static int getMateAlignmentEnd(final SAMRecord read)
    {
        String mateCigarStr = read.getStringAttribute(MATE_CIGAR_ATTRIBUTE);
        if(mateCigarStr == null || mateCigarStr.equals(NO_CIGAR))
            return NO_POSITION;

        return getMateAlignmentEnd(read.getMateAlignmentStart(), mateCigarStr);
    }

    public static int getMateAlignmentEnd(final int readStart, @NotNull final String cigarStr)
    {
        return getReadBoundaryPosition(readStart, cigarStr, false, false);
    }

    public static int getFivePrimeUnclippedPosition(final int readStart, @NotNull final String cigarStr, final boolean forwardStrand)
    {
        return getReadBoundaryPosition(readStart, cigarStr, forwardStrand, true);
    }

    public static int getThreePrimeUnclippedPosition(final int readStart, @NotNull final String cigarStr, final boolean forwardStrand)
    {
        return getReadBoundaryPosition(readStart, cigarStr, !forwardStrand, true);
    }

    public static int getFivePrimeUnclippedPosition(final SAMRecord read) { return getUnclippedPosition(read, true); }

    public static int getThreePrimeUnclippedPosition(final SAMRecord read) { return getUnclippedPosition(read, false); }

    public static int getUnclippedPosition(final SAMRecord read, boolean fivePrime)
    {
        // returns the 5' or 3' position of the read, factoring in any soft-clipped bases
        int position;

        if((orientation(read) == ORIENT_FWD) == fivePrime)
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

    public static List<int[]> generateMappedCoords(final Cigar cigar, int posStart)
    {
        final List<int[]> mappedCoords = Lists.newArrayList();

        // first establish whether the read is split across 2 distant regions, and if so which it maps to
        int posOffset = 0;
        boolean continueRegion = false;

        for(CigarElement element : cigar.getCigarElements())
        {
            if(element.getOperator() == CigarOperator.S)
            {
                // nothing to skip
            }
            else if(element.getOperator() == D)
            {
                posOffset += element.getLength();
                continueRegion = true;
            }
            else if(element.getOperator() == CigarOperator.I)
            {
                // nothing to skip
                continueRegion = true;
            }
            else if(element.getOperator() == CigarOperator.N)
            {
                posOffset += element.getLength();
                continueRegion = false;
            }
            else if(element.getOperator() == CigarOperator.M)
            {
                int readStartPos = posStart + posOffset;
                int readEndPos = readStartPos + element.getLength() - 1;

                if(continueRegion && !mappedCoords.isEmpty())
                {
                    int[] lastRegion = mappedCoords.get(mappedCoords.size() - 1);
                    lastRegion[SE_END] = readEndPos;
                }
                else
                {
                    mappedCoords.add(new int[] { readStartPos, readEndPos });
                }

                posOffset += element.getLength();
                continueRegion = false;
            }
        }

        return mappedCoords;
    }

    public static String getOrientationString(final SAMRecord read)
    {
        if(!read.getReadPairedFlag())
            return "";

        if(read.getReadUnmappedFlag() || read.getMateUnmappedFlag())
            return "";

        String firstStr;
        String secondStr;
        if(read.getFirstOfPairFlag())
        {
            firstStr = read.getReadNegativeStrandFlag() ? "R" : "F";
            secondStr = read.getMateNegativeStrandFlag() ? "R" : "F";
        }
        else
        {
            firstStr = read.getMateNegativeStrandFlag() ? "R" : "F";
            secondStr = read.getReadNegativeStrandFlag() ? "R" : "F";
        }

        return format("%s1%s2", firstStr, secondStr);
    }

    public static String readToString(final SAMRecord read)
    {
        return format("id(%s) coords(%s:%d-%d) cigar(%s) mate(%s:%d) flags(%d)",
                read.getReadName(), read.getContig(), read.getAlignmentStart(), read.getAlignmentEnd(),
                read.getCigarString(), read.getMateReferenceName(), read.getMateAlignmentStart(), read.getFlags());
    }
}
