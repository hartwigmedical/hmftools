package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.common.purple.PurpleCopyNumber.buildChromosomeMap;
import static com.hartwig.hmftools.purple.copynumber.LohCalcData.LOH_NONE;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;

public final class LohCalcs
{
    public static final double MIN_LOH_CN = 0.5;
    public static final int LOH_SHORT_INSERT_LENGTH = 1000;

    public static LohCalcData calcLohData(final List<PurpleCopyNumber> copyNumbers)
    {
        Map<String,List<PurpleCopyNumber>> chrCopyNumberMap = buildChromosomeMap(copyNumbers);

        LohCalcData totalLohData = new LohCalcData();

        for(Map.Entry<String,List<PurpleCopyNumber>> entry : chrCopyNumberMap.entrySet())
        {
            totalLohData.add(calcLohSegments(entry.getKey(), entry.getValue(), LOH_SHORT_INSERT_LENGTH, false, true));
        }

        return totalLohData;
    }

    public static LohCalcData calcLohSegments(
            final String chromosomeStr, final List<PurpleCopyNumber> copyNumbers, double minLohLength, boolean allowWholeArm,
            boolean countByArm)
    {
        HumanChromosome chromosome = HumanChromosome.fromString(chromosomeStr);

        if(chromosome.isAllosome())
            return LOH_NONE;

        int lohCount = 0;
        int lohBaseCount = 0;
        int totalBaseCount = 0;
        boolean inLohSegment = false;
        int lohSegmentStart = 0;
        boolean beforeCentromere = false;
        boolean hasNonLohSegment = false;

        for(int i = 0; i < copyNumbers.size(); ++i)
        {
            PurpleCopyNumber copyNumber = copyNumbers.get(i);

            if(copyNumber.segmentStartSupport() == SegmentSupport.TELOMERE)
            {
                // reset at start of new chromosome
                inLohSegment = false;
                beforeCentromere = true;
            }
            else if(copyNumber.segmentStartSupport() == SegmentSupport.CENTROMERE)
            {
                beforeCentromere = false;
            }

            if(beforeCentromere && chromosome.hasShortArm())
                continue;

            totalBaseCount += copyNumber.length();

            PurpleCopyNumber nextCopyNumber = i < copyNumbers.size() - 1 ? copyNumbers.get(i + 1) : null;

            double minAlleleCn = copyNumber.minorAlleleCopyNumber();
            boolean isLohSegment = minAlleleCn < MIN_LOH_CN;
            hasNonLohSegment |= !isLohSegment;

            boolean isLohNextSegment = nextCopyNumber != null && nextCopyNumber.minorAlleleCopyNumber() < MIN_LOH_CN;

            if(!inLohSegment && isLohSegment)
            {
                inLohSegment = true;
                lohSegmentStart = copyNumber.start();
            }

            boolean endOfArm = copyNumber.segmentEndSupport() == SegmentSupport.TELOMERE
                    || (countByArm && copyNumber.segmentEndSupport() == SegmentSupport.CENTROMERE);

            boolean checkLohEnd = false;

            if(inLohSegment)
            {
                if(endOfArm)
                {
                    checkLohEnd = true;
                }
                else if(!isLohNextSegment)
                {
                    if(nextCopyNumber != null && nextCopyNumber.length() <= LOH_SHORT_INSERT_LENGTH)
                        checkLohEnd = false;
                    else
                        checkLohEnd = true;
                }
            }

            if(checkLohEnd)
            {
                int endPosition = copyNumber.end();
                int segmentLength = endPosition - lohSegmentStart + 1;

                if(segmentLength > minLohLength)
                {
                    /*
                    PPL_LOGGER.trace("LOH segment({}:{}-{}) length({})",
                            copyNumber.chromosome(), lohSegmentStart, copyNumber.end(), segmentLength);
                    */

                    ++lohCount;
                    lohBaseCount += segmentLength;
                }

                inLohSegment = false;
            }
        }

        if(allowWholeArm)
        {
            if(lohCount == 0)
                return LOH_NONE;

            if(!hasNonLohSegment)
                lohCount = 1; // convert an LOH across both arms to just a single LOH

            return new LohCalcData(lohBaseCount, totalBaseCount, lohCount);
        }
        else
        {
            return hasNonLohSegment ? new LohCalcData(lohBaseCount, totalBaseCount, lohCount) : LOH_NONE;

        }
    }
}
