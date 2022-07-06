package com.hartwig.hmftools.purple.tools;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

public class HrdDetection
{
    private final int mSegmentBreakLength;
    private final int mSegmentImbalanceLength;
    private final int mLohMinLength;
    private final double mPloidyModifier;

    private static final int DEFAULT_SEGMENT_IMBALANCE_LENGTH = 11_000_000;
    private static final int DEFAULT_LOH_MIN_LENGTH = 15_000_000;
    private static final int DEFAULT_SEGMENT_BREAK_LENGTH = 10_000_000;
    private static final double DEFAULT_PLOIDY_FACTOR = 15.5;

    private static final int LOH_SHORT_INSERT_LENGTH = 1000;
    private static final double MIN_LOH_CN = 0.5;
    private static final double MAX_COPY_NUM_DIFF = 0.5;
    private static final double MAX_COPY_NUM_DIFF_PERC = 0.15;

    private static final List<String> LOH_SHORT_ARM_CHROMOSOMES = Lists.newArrayList("13", "14", "15", "21", "22");
    private static final List<String> LOH_IGNORE_CHROMOSOMES = Lists.newArrayList("X", "Y");

    public HrdDetection()
    {
        this(DEFAULT_SEGMENT_BREAK_LENGTH, DEFAULT_SEGMENT_IMBALANCE_LENGTH, DEFAULT_LOH_MIN_LENGTH, DEFAULT_PLOIDY_FACTOR);
    }

    public HrdDetection(final int segmentBreakLength, final int segmentImbalanceLength, final int lohMinLength, final double ploidyModifier)
    {
        mSegmentBreakLength = segmentBreakLength;
        mSegmentImbalanceLength = segmentImbalanceLength;
        mLohMinLength = lohMinLength;
        mPloidyModifier = ploidyModifier;
    }

    public HrdData calculateHrdData(final List<PurpleCopyNumber> copyNumbers, final double samplePloidy)
    {
        int lohSegments = calcLohSegments(copyNumbers);

        int unbalancedSegments = calcSegmentImbalances(copyNumbers);

        double segmentBreaks = calcSegmentBreaks(copyNumbers, samplePloidy);

        return new HrdData(lohSegments, unbalancedSegments, segmentBreaks);
    }

    public int calcLohSegments(final List<PurpleCopyNumber> copyNumbers)
    {
        int lohCount = 0;
        boolean inLohSegment = false;
        int lohSegmentStart = 0;
        boolean beforeCentromere = false;

        for(int i = 0; i < copyNumbers.size(); ++i)
        {
            PurpleCopyNumber copyNumber = copyNumbers.get(i);

            if(LOH_IGNORE_CHROMOSOMES.contains(RefGenomeFunctions.stripChrPrefix(copyNumber.chromosome())))
                continue;

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

            if(beforeCentromere && LOH_SHORT_ARM_CHROMOSOMES.contains(RefGenomeFunctions.stripChrPrefix(copyNumber.chromosome())))
                continue;

            PurpleCopyNumber nextCopyNumber = i < copyNumbers.size() - 1 ? copyNumbers.get(i + 1) : null;

            double minAlleleCn = copyNumber.minorAlleleCopyNumber();
            boolean isLohSegment = minAlleleCn < MIN_LOH_CN;

            boolean isLohNextSegment = nextCopyNumber != null && nextCopyNumber.minorAlleleCopyNumber() < MIN_LOH_CN;


            if(!inLohSegment && isLohSegment)
            {
                inLohSegment = true;
                lohSegmentStart = copyNumber.start();
            }

            boolean endOfArm = copyNumber.segmentEndSupport() == SegmentSupport.TELOMERE
                    || copyNumber.segmentEndSupport() == SegmentSupport.CENTROMERE;

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
                // int endPosition = endOfArm ? copyNumber.end() : copyNumber.start();
                int endPosition = copyNumber.end();
                int segmentLength = endPosition - lohSegmentStart - 1;

                if(segmentLength >= mLohMinLength)
                {
                    PPL_LOGGER.trace("LOH segment({}:{}-{}) length({})",
                            copyNumber.chromosome(), lohSegmentStart, copyNumber.end(), segmentLength);
                    ++lohCount;
                }

                inLohSegment = false;
            }
        }

        return lohCount;
    }

    public int calcSegmentImbalances(final List<PurpleCopyNumber> copyNumbers)
    {
        int imbalances = 0;

        int telomereStartImbalanceStart = 0;
        int postCentromenereImbalanceStart = 0;
        boolean isPostCentromere = false;

        for(PurpleCopyNumber copyNumber : copyNumbers)
        {
            boolean hasImbalance = hasImbalance(copyNumber);

            // consider the starting position
            if(copyNumber.segmentStartSupport() == SegmentSupport.TELOMERE)
            {
                isPostCentromere = false;
                postCentromenereImbalanceStart = 0;

                if(hasImbalance)
                {
                    telomereStartImbalanceStart = 1;
                }
            }
            else if(copyNumber.segmentStartSupport() == SegmentSupport.CENTROMERE)
            {
                // cannot start at the centromere
                isPostCentromere = true;
            }
            else if(isPostCentromere && hasImbalance && postCentromenereImbalanceStart == 0)
            {
                postCentromenereImbalanceStart = copyNumber.start();
            }

            // consider the end position
            if(copyNumber.segmentEndSupport() == SegmentSupport.CENTROMERE)
            {
                if(telomereStartImbalanceStart > 0 && !hasImbalance)
                {
                    // check the previous segment(s)
                    int segmentLength = copyNumber.start() - postCentromenereImbalanceStart;

                    if(segmentLength >= mSegmentImbalanceLength && !hasImbalance)
                    {
                        PPL_LOGGER.trace("imbalance from telomere segment({}:{}-{}) length({})",
                                copyNumber.chromosome(), postCentromenereImbalanceStart, copyNumber.end(), segmentLength);

                        ++imbalances;
                    }
                }

                // cannot extend to the centromere
                telomereStartImbalanceStart = 0;
            }
            else if(copyNumber.segmentEndSupport() == SegmentSupport.TELOMERE)
            {
                if(postCentromenereImbalanceStart > 0 && hasImbalance)
                {
                    int segmentLength = copyNumber.end() - postCentromenereImbalanceStart;

                    if(segmentLength >= mSegmentImbalanceLength)
                    {
                        PPL_LOGGER.trace("imbalance to telomere segment({}:{}-{}) length({})",
                                copyNumber.chromosome(), postCentromenereImbalanceStart, copyNumber.end(), segmentLength);

                        ++imbalances;
                    }
                }
            }
            else
            {
                if(telomereStartImbalanceStart > 0)
                {
                    int segmentLength = copyNumber.end() - postCentromenereImbalanceStart;

                    if(segmentLength >= mSegmentImbalanceLength && !hasImbalance)
                    {
                        PPL_LOGGER.trace("imbalance from telomere segment({}:{}-{}) length({})",
                                copyNumber.chromosome(), postCentromenereImbalanceStart, copyNumber.end(), segmentLength);

                        ++imbalances;
                        telomereStartImbalanceStart = 0;
                        continue;
                    }

                    if(hasImbalance)
                        continue;

                    if(copyNumber.length() <= LOH_SHORT_INSERT_LENGTH)
                        continue;

                    // too short so cancel this possible segment
                    telomereStartImbalanceStart = 0;
                }
                else if(postCentromenereImbalanceStart > 0 && !hasImbalance)
                {
                    if(copyNumber.length() <= LOH_SHORT_INSERT_LENGTH)
                        continue;

                    postCentromenereImbalanceStart = 0;
                }
            }
        }

        return imbalances;
    }

    private static boolean hasImbalance(final PurpleCopyNumber copyNumber)
    {
        double minCopyNumber = copyNumber.minorAlleleCopyNumber();
        double maxCopyNumber = copyNumber.majorAlleleCopyNumber();

        return !copyNumbersEqual(minCopyNumber, maxCopyNumber);
    }

    public double calcSegmentBreaks(final List<PurpleCopyNumber> copyNumbers, final double ploidy)
    {
        int segmentBreaks = 0;

        double cnBaseTotal = 0;
        int segmentStart = 0;
        double previousSegmentCn = 0;
        boolean hasPreviousSegment = false;

        for(PurpleCopyNumber copyNumber : copyNumbers)
        {
            boolean startNewSegment = false;

            if(copyNumber.segmentStartSupport() == SegmentSupport.TELOMERE
            || copyNumber.segmentStartSupport() == SegmentSupport.CENTROMERE)
            {
                startNewSegment = true;
            }
            else if(copyNumber.segmentEndSupport() == SegmentSupport.CENTROMERE
            || copyNumber.segmentEndSupport() == SegmentSupport.TELOMERE)
            {
                // don't compare the reigons straddling the centromere or onto the next chromosome
                cnBaseTotal = 0;
            }
            else
            {
                int previousSegmentEnd = copyNumber.start() - 1;
                int segmentLength = previousSegmentEnd - segmentStart;

                if(segmentLength >= mSegmentBreakLength)
                {
                    double segmentCn = cnBaseTotal / segmentLength;
                    if(hasPreviousSegment)
                    {
                        if(!copyNumbersEqual(segmentCn, previousSegmentCn))
                        {
                            PPL_LOGGER.trace(format("chr(%s) segment break: previous(cn=%.2f) current(len=%d cn=%.2f)",
                                    copyNumber.chromosome(), previousSegmentCn, segmentLength, segmentCn));

                            ++segmentBreaks;
                        }
                    }

                    startNewSegment = true;
                    previousSegmentCn = segmentCn;
                }

                if(cnBaseTotal > 0)
                {
                    previousSegmentCn = cnBaseTotal / segmentLength;
                }
            }

            if(startNewSegment)
            {
                segmentStart = copyNumber.start();
                cnBaseTotal = 0;
                hasPreviousSegment = true;
            }

            cnBaseTotal += copyNumber.length() * copyNumber.averageTumorCopyNumber();

        }

        double adjustedSegmentBreaks = segmentBreaks - mPloidyModifier * ploidy;

        return adjustedSegmentBreaks;
    }

    private static boolean copyNumbersEqual(double cn1, double cn2)
    {
        double copyNumDiff = abs(cn2 - cn1);
        double copyNumDiffPerc = copyNumDiff / max(abs(cn1), abs(cn2));

        if (copyNumDiff > MAX_COPY_NUM_DIFF && copyNumDiffPerc > MAX_COPY_NUM_DIFF_PERC)
            return false;

        return true;
    }

}
