package com.hartwig.hmftools.common.samtools;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;

import static htsjdk.samtools.CigarOperator.D;

import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public final class SamRecordUtils
{
    public static final int PHRED_OFFSET = 33;

    public static final String SUPPLEMENTARY_ATTRIBUTE = "SA";
    public static final String SECONDARY_ATTRIBUTE = "HI";
    public static final String NO_MATE_CHROMOSOME = "*";

    public static final Logger SAM_LOGGER = LogManager.getLogger(SamRecordUtils.class);

    public static int getBaseQuality(final char quality)
    {
        return quality - PHRED_OFFSET;
    }

    public static int getBaseQuality(final SAMRecord record, int readPosition)
    {
        return getAvgBaseQuality(record, readPosition, 1);
    }

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

    public static boolean properPair(final SAMRecord record)
    {
        return record.getReadPairedFlag() && record.getProperPairFlag();
    }

    public static boolean mateUnmapped(final SAMRecord record)
    {
        return record.getReadPairedFlag() && record.getMateUnmappedFlag();
    }

    public static int getAvgBaseQuality(final SAMRecord record, int readPosition, int length)
    {
        assert (readPosition > 0);

        int score = 0;
        final String baseQualities = record.getBaseQualityString();
        for(int index = readPosition - 1; index < Math.min(readPosition - 1 + length, baseQualities.length()); index++)
        {
            int baseScore = getBaseQuality(baseQualities.charAt(index));
            score += baseScore;
        }
        return score / length;
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
}
