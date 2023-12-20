package com.hartwig.hmftools.purple.hrd;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.CHR_PREFIX;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.purple.PurpleCopyNumber.buildChromosomeMap;
import static com.hartwig.hmftools.common.purple.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.SegmentSupport.TELOMERE;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.copynumber.LohCalcs.LOH_SHORT_INSERT_LENGTH;
import static com.hartwig.hmftools.purple.hrd.HrdStatus.HRD_DEFICIENT_HIGH_CONFIDENCE;
import static com.hartwig.hmftools.purple.hrd.HrdStatus.HRD_DEFICIENT_LOW_CONFIDENCE;
import static com.hartwig.hmftools.purple.hrd.HrdStatus.HRD_PROFICIENT_HIGH_CONFIDENCE;
import static com.hartwig.hmftools.purple.hrd.HrdStatus.HRD_PROFICIENT_LOW_CONFIDENCE;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.purple.copynumber.LohCalcData;
import com.hartwig.hmftools.purple.copynumber.LohCalcs;

public class HrdDetection
{
    private final int mSegmentImbalanceLength;
    private final int mLohMinLength;

    private static final int DEFAULT_SEGMENT_IMBALANCE_LENGTH = 11_000_000;
    private static final int SEGMENT_IMBALANCE_MIN_SEGMENT_LENGTH = 100;

    private static final int DEFAULT_LOH_MIN_LENGTH = 15_000_000;

    private static final int DEFAULT_SEGMENT_BREAK_MIN_LENGTH = 100_000;
    private static final int DEFAULT_SEGMENT_BREAK_LENGTH = 3_000_000;

    private static final int HRD_STATUS_LOH_CUTOFF = 8;
    private static final int HRD_STATUS_IMBALANCE_CUTOFF = 7;
    private static final int HRD_STATUS_SEGMENT_BREAKS_CUTOFF = 35;

    private static final double MAX_COPY_NUM_DIFF = 0.5;
    private static final double MAX_COPY_NUM_DIFF_PERC = 0.2;

    public HrdDetection()
    {
        this(DEFAULT_SEGMENT_IMBALANCE_LENGTH, DEFAULT_LOH_MIN_LENGTH);
    }

    public HrdDetection(final int segmentImbalanceLength, final int lohMinLength)
    {
        mSegmentImbalanceLength = segmentImbalanceLength;
        mLohMinLength = lohMinLength;
    }

    public static HrdStatus determineHrdStatus(final HrdData hrdData)
    {
        /*
            If 2 of 3 measures exceed their cutoffs, it is HRD DEFICIENT (LOW CONFIDENCE)
            If 3 of 3 measures exceed their cutoffs, it is HRD DEFICIENT (HIGH CONFIDENCE)

            If 2 of 3 measures equal or are below their cutoffs, it is HRD PROFICIENT (LOW CONFIDENCE)
            If 3 of 3 measures equal or are below their cutoffs, it is HRD PROFICIENT (HIGH CONFIDENCE)
         */

        int negativeCount = 0;

        if(hrdData.LohSegments > HRD_STATUS_LOH_CUTOFF)
            ++negativeCount;

        if(hrdData.SegmentImbalances > HRD_STATUS_IMBALANCE_CUTOFF)
            ++negativeCount;

        if(hrdData.SegmentBreaks > HRD_STATUS_SEGMENT_BREAKS_CUTOFF)
            ++negativeCount;

        if(negativeCount == 3)
            return HRD_DEFICIENT_HIGH_CONFIDENCE;
        else if(negativeCount == 2)
            return HRD_DEFICIENT_LOW_CONFIDENCE;
        else if(negativeCount == 1)
            return HRD_PROFICIENT_LOW_CONFIDENCE;
        else
            return HRD_PROFICIENT_HIGH_CONFIDENCE;
    }

    public HrdData calculateHrdData(final List<PurpleCopyNumber> copyNumbers)
    {
        if(copyNumbers.isEmpty())
            return new HrdData(0, 0, 0);

        String exampleChr = copyNumbers.get(0).chromosome();
        RefGenomeVersion refGenomeVersion = exampleChr.startsWith(CHR_PREFIX) ? V38 : V37;

        // build a chr-segment map
        Map<String,List<PurpleCopyNumber>> chrCopyNumberMap = buildChromosomeMap(copyNumbers);

        // filter copy numbers to required types
        for(List<PurpleCopyNumber> chrCopyNumbers : chrCopyNumberMap.values())
        {
            int index = 0;
            while(index < chrCopyNumbers.size())
            {
                PurpleCopyNumber copyNumber = chrCopyNumbers.get(index);

                if(copyNumber.method() != CopyNumberMethod.BAF_WEIGHTED || HumanChromosome.fromString(copyNumber.chromosome()).isAllosome())
                {
                    chrCopyNumbers.remove(index);
                }
                else
                {
                    ++index;
                }
            }
        }

        int totalLohSegments = 0;
        int totalUnbalancedSegments = 0;
        int totalSegmentBreaks = 0;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            if(chromosome.isAllosome())
                continue;

            String chrStr = refGenomeVersion.versionedChromosome(chromosome.toString());
            List<PurpleCopyNumber> chrCopyNumbers = chrCopyNumberMap.get(chrStr);

            if(chrCopyNumbers == null || chrStr.isEmpty())
                continue;

            totalLohSegments += calcLohSegments(chrStr, chrCopyNumbers);

            totalUnbalancedSegments += calcSegmentImbalances(chrStr, chrCopyNumbers);

            totalSegmentBreaks += calcSegmentBreaks(chrStr, chrCopyNumbers);
        }

        return new HrdData(totalLohSegments, totalUnbalancedSegments, totalSegmentBreaks);
    }

    public int calcLohSegments(final String chromosome, final List<PurpleCopyNumber> copyNumbers)
    {
        LohCalcData lohCalcData = LohCalcs.calcLohSegments(chromosome, copyNumbers, mLohMinLength, false);

        return lohCalcData.Segments;
    }

    public int calcSegmentImbalances(final String chromosome, final List<PurpleCopyNumber> copyNumbers)
    {
        int imbalances = 0;

        // min CN segment length 100 - skip or remove?

        List<PurpleCopyNumber> localCopyNumbers = copyNumbers.stream()
                .filter(x -> x.end() - x.start() + 1 > SEGMENT_IMBALANCE_MIN_SEGMENT_LENGTH).collect(Collectors.toList());

        int copyNumberCount = localCopyNumbers.size();

        for(int i = 0; i <= 1; ++i)
        {
            boolean pArm = (i == 0);

            int index = pArm ? 0 : copyNumberCount - 1;

            int imbalanceStartPosition = -1;

            while(index >= 0 && index < copyNumberCount)
            {
                PurpleCopyNumber copyNumber = localCopyNumbers.get(index);

                // break on centromere
                if(pArm && copyNumber.segmentStartSupport() == CENTROMERE)
                    break;
                else if(!pArm && copyNumber.segmentEndSupport() == CENTROMERE)
                    break;

                boolean hasImbalance = hasImbalance(copyNumber);

                if(copyNumber.segmentStartSupport() == TELOMERE || copyNumber.segmentEndSupport() == TELOMERE)
                {
                    // must start with an imbalance
                    if(!hasImbalance)
                        break;

                    imbalanceStartPosition = pArm ? copyNumber.start() : copyNumber.end();
                }
                else
                {
                    if(!hasImbalance)
                    {
                        // imbalance has finished, check it is long enough
                        int segmentStart;
                        int segmentEnd;

                        if(pArm)
                        {
                            segmentStart = imbalanceStartPosition;
                            segmentEnd = copyNumber.start(); // ended on the previous segment
                        }
                        else
                        {
                            segmentStart = copyNumber.end();
                            segmentEnd = imbalanceStartPosition;
                        }

                        int segmentLength = segmentEnd - segmentStart;

                        if(segmentLength >= mSegmentImbalanceLength)
                        {
                            PPL_LOGGER.trace("imbalance on chr({}) arm({}) segment({}-{}) length({})",
                                    chromosome, pArm ? "P" : "Q", segmentStart, segmentEnd, segmentLength);

                            ++imbalances;
                        }

                        break;
                    }
                }

                if(pArm)
                    ++index;
                else
                    --index;
            }
        }

        /*
        int telomereStartImbalanceStart = 0;
        int postCentromenereImbalanceStart = 0;
        boolean isPostCentromere = false;

        for(PurpleCopyNumber copyNumber : localCopyNumbers)
        {
            boolean hasImbalance = hasImbalance(copyNumber);

            // consider the starting position
            if(copyNumber.segmentStartSupport() == TELOMERE)
            {
                isPostCentromere = false;
                postCentromenereImbalanceStart = 0;

                if(hasImbalance)
                {
                    telomereStartImbalanceStart = 1;
                }
            }
            else if(copyNumber.segmentStartSupport() == CENTROMERE)
            {
                // cannot start at the centromere and if not already
                isPostCentromere = true;
            }
            else if(isPostCentromere && hasImbalance && postCentromenereImbalanceStart == 0)
            {
                postCentromenereImbalanceStart = copyNumber.start();
            }

            // consider the end position
            if(copyNumber.segmentEndSupport() == CENTROMERE)
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
            else if(copyNumber.segmentEndSupport() == TELOMERE)
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
        */

        return imbalances;
    }

    private static boolean hasImbalance(final PurpleCopyNumber copyNumber)
    {
        double minCopyNumber = copyNumber.minorAlleleCopyNumber();
        double maxCopyNumber = copyNumber.majorAlleleCopyNumber();

        return !copyNumbersEqualAbsAndRelative(minCopyNumber, maxCopyNumber);
    }

    public double calcSegmentBreaks(final String chromosome, final List<PurpleCopyNumber> copyNumbers)
    {
        if(copyNumbers.isEmpty())
            return 0;

        int segmentBreaks = 0;

        // covert to a simpler and modifiable form, and merge small segments

        CopyNumberSegment currentSegment = null;
        List<CopyNumberSegment> segments = Lists.newArrayList();

        // first create segments if the length exceeds the minimum 100K length
        for(PurpleCopyNumber copyNumber : copyNumbers)
        {
            if(currentSegment == null)
            {
                currentSegment = new CopyNumberSegment(copyNumber);
                segments.add(currentSegment);
            }
            else
            {
                if(copyNumber.end() - currentSegment.start() < DEFAULT_SEGMENT_BREAK_MIN_LENGTH
                || copyNumber.length() < DEFAULT_SEGMENT_BREAK_MIN_LENGTH)
                {
                    // extend the current segment
                    currentSegment.setEnd(copyNumber.end());
                }
                else
                {
                    currentSegment = new CopyNumberSegment(copyNumber);
                    segments.add(currentSegment);
                }
            }
        }

        // next merge segments with small CN transitions
        int index = 0;
        while(index < segments.size() - 1)
        {
            CopyNumberSegment segment = segments.get(index);
            CopyNumberSegment nextSegment = segments.get(index + 1);

            if(!copyNumbersEqualAbs(segment.CopyNumber, nextSegment.CopyNumber))
            {
                ++index;
                continue;
            }

            // merge and average these segments
            mergeSegments(segment, nextSegment);
            segments.remove(index + 1);

            // stay with current
        }

        // now iteratively merge regions less than 3MB if the surrounding regions have no net change
        while(segments.size() > 1)
        {
            CopyNumberSegment lowerSegment = null;
            CopyNumberSegment upperSegment = null;
            CopyNumberSegment shortedSegment = null;

            for(int i = 0; i < segments.size(); ++i)
            {
                CopyNumberSegment segment = segments.get(i);

                if(segment.length() > DEFAULT_SEGMENT_BREAK_LENGTH)
                    continue;

                CopyNumberSegment prevSegment = i > 0 ? segments.get(i - 1) : null;
                CopyNumberSegment nextSegment = i < segments.size() - 1 ? segments.get(i + 1) : null;

                if(shortedSegment == null || shortedSegment.length() > segment.length())
                {
                    shortedSegment = segment;
                    lowerSegment = prevSegment;
                    upperSegment = nextSegment;
                }
            }

            if(shortedSegment == null) // no segments to merge
                break;

            segments.remove(shortedSegment);

            // if the 2 surrounding regions are within CN threshold then merge, else assign the removed region equally amongst the surrounding regions
            if(upperSegment != null && lowerSegment != null)
            {
                if(copyNumbersEqualAbs(upperSegment.CopyNumber, lowerSegment.CopyNumber))
                {
                    mergeSegments(lowerSegment, upperSegment);
                    segments.remove(upperSegment);
                }
                else
                {
                    lowerSegment.setEnd(shortedSegment.end());
                }
            }
            else if(lowerSegment != null)
            {
                lowerSegment.setEnd(shortedSegment.end());
                lowerSegment.Segments[SE_END] = shortedSegment.segmentSupportEnd();
            }
            else if(upperSegment != null)
            {
                upperSegment.setStart(shortedSegment.start());
                upperSegment.Segments[SE_START] = shortedSegment.segmentSupportStart();
            }
        }

        // now count up the breaks
        for(int i = 0; i < segments.size() - 1; ++i)
        {
            CopyNumberSegment segment = segments.get(i);
            CopyNumberSegment nextSegment = segments.get(i + 1);

            if(!copyNumbersEqualAbs(segment.CopyNumber, nextSegment.CopyNumber))
            {
                ++segmentBreaks;
            }
        }

        PPL_LOGGER.debug("chr({}) segmentBreaks({}})", chromosome, segmentBreaks);

        if(PPL_LOGGER.isTraceEnabled())
        {
            for(CopyNumberSegment segment : segments)
            {
                PPL_LOGGER.trace(format("segment: %s", segment));
            }
        }

        return segmentBreaks;
    }

    private static void mergeSegments(final CopyNumberSegment lower, final CopyNumberSegment upper)
    {
        if(lower.start() > upper.start())
            return;

        long combinedLength = lower.length() + upper.length();
        double totalCopyNumber = lower.length() * lower.CopyNumber + upper.length() * upper.CopyNumber;
        double averageCopyNumber = totalCopyNumber / combinedLength;
        lower.setEnd(upper.end());
        lower.Segments[SE_END] = upper.segmentSupportEnd();
        lower.CopyNumber = averageCopyNumber;
    }

    private static boolean copyNumbersEqualAbs(double cn1, double cn2)
    {
        double copyNumDiff = abs(cn2 - cn1);
        return copyNumDiff < MAX_COPY_NUM_DIFF;
    }

    private static boolean copyNumbersEqualAbsAndRelative(double cn1, double cn2)
    {
        double copyNumDiff = abs(cn2 - cn1);
        double copyNumDiffPerc = copyNumDiff / max(abs(cn1), abs(cn2));

        if(copyNumDiff > MAX_COPY_NUM_DIFF && copyNumDiffPerc > MAX_COPY_NUM_DIFF_PERC)
            return false;

        return true;
    }

    private class CopyNumberSegment
    {
        public final ChrBaseRegion Region;
        public double CopyNumber;
        public final SegmentSupport[] Segments;

        public CopyNumberSegment(final PurpleCopyNumber copyNumber)
        {
            Region = new ChrBaseRegion(copyNumber.chromosome(), copyNumber.start(), copyNumber.end());
            CopyNumber = copyNumber.averageTumorCopyNumber();
            Segments = new SegmentSupport[] { copyNumber.segmentStartSupport(), copyNumber.segmentEndSupport() };
        }

        public String chromosome() { return Region.Chromosome; }

        public int start() { return Region.start(); }
        public void setStart(int pos) { Region.setStart(pos); }

        public int end() { return Region.end(); }
        public void setEnd(int pos) { Region.setEnd(pos); }

        public int length() { return Region.baseLength(); }

        public SegmentSupport segmentSupportStart() { return Segments[SE_START]; }
        public SegmentSupport segmentSupportEnd() { return Segments[SE_END]; }

        public String toString()
        {
            return format("region(%s) segs(%s-%s) copyNumber(%.1f) length(%d)",
                Region.toString(), Segments[SE_START], Segments[SE_END], CopyNumber, length());
        }
    }

}
