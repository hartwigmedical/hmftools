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
import static com.hartwig.hmftools.common.purple.HrdStatus.HRD_DEFICIENT_HIGH_CONFIDENCE;
import static com.hartwig.hmftools.common.purple.HrdStatus.HRD_DEFICIENT_LOW_CONFIDENCE;
import static com.hartwig.hmftools.common.purple.HrdStatus.HRD_PROFICIENT_HIGH_CONFIDENCE;
import static com.hartwig.hmftools.common.purple.HrdStatus.HRD_PROFICIENT_LOW_CONFIDENCE;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.HrdData;
import com.hartwig.hmftools.common.purple.HrdStatus;
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

    private static final int HRD_STATUS_LOH_CUTOFF = 7;
    private static final int HRD_STATUS_IMBALANCE_CUTOFF = 2;
    private static final int HRD_STATUS_SEGMENT_BREAKS_CUTOFF = 24;

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

    public static HrdStatus determineHrdStatus(final int lohSegments, final int segmentImbalances, final int segmentBreaks)
    {
        /*
            If 2 of 3 measures exceed their cutoffs, it is HRD DEFICIENT (LOW CONFIDENCE)
            If 3 of 3 measures exceed their cutoffs, it is HRD DEFICIENT (HIGH CONFIDENCE)

            If 2 of 3 measures equal or are below their cutoffs, it is HRD PROFICIENT (LOW CONFIDENCE)
            If 3 of 3 measures equal or are below their cutoffs, it is HRD PROFICIENT (HIGH CONFIDENCE)
         */

        int negativeCount = 0;

        if(lohSegments > HRD_STATUS_LOH_CUTOFF)
            ++negativeCount;

        if(segmentImbalances > HRD_STATUS_IMBALANCE_CUTOFF)
            ++negativeCount;

        if(segmentBreaks > HRD_STATUS_SEGMENT_BREAKS_CUTOFF)
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
            return new HrdData(0, 0, 0, HrdStatus.UNKNOWN);

        String exampleChr = copyNumbers.get(0).chromosome();
        RefGenomeVersion refGenomeVersion = exampleChr.startsWith(CHR_PREFIX) ? V38 : V37;

        // build a chr-segment map
        Map<String,List<PurpleCopyNumber>> chrCopyNumberMap = buildChromosomeMap(copyNumbers);

        int totalLohSegments = 0;
        int totalUnbalancedSegments = 0;
        int totalSegmentBreaks = 0;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            if(chromosome.isAllosome())
                continue;

            String chrStr = refGenomeVersion.versionedChromosome(chromosome.toString());

            // filter copy numbers to required types
            List<PurpleCopyNumber> chrCopyNumbers = chrCopyNumberMap.get(chrStr);

            if(chrCopyNumbers == null)
                continue;

            chrCopyNumbers = chrCopyNumbers.stream().filter(x -> x.method() == CopyNumberMethod.BAF_WEIGHTED).collect(Collectors.toList());

            if(chrCopyNumbers.isEmpty())
                continue;

            totalLohSegments += calcLohSegments(chrStr, chrCopyNumbers);

            totalUnbalancedSegments += calcSegmentImbalances(chrStr, chrCopyNumbers);

            totalSegmentBreaks += calcSegmentBreaks(chrStr, chrCopyNumbers);
        }

        HrdStatus hrdStatus = determineHrdStatus(totalLohSegments, totalUnbalancedSegments, totalSegmentBreaks);

        return new HrdData(totalLohSegments, totalUnbalancedSegments, totalSegmentBreaks, hrdStatus);
    }

    public int calcLohSegments(final String chromosome, final List<PurpleCopyNumber> copyNumbers)
    {
        LohCalcData lohCalcData = LohCalcs.calcLohSegments(chromosome, copyNumbers, mLohMinLength, false, false);

        return lohCalcData.Segments;
    }

    public int calcSegmentImbalances(final String chromosome, final List<PurpleCopyNumber> copyNumbers)
    {
        int imbalances = 0;

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

        PPL_LOGGER.trace("chr({}) segmentBreaks({})", chromosome, segmentBreaks);

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
