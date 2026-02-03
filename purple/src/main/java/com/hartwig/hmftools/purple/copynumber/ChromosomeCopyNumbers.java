package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.common.segmentation.Arm.P;
import static com.hartwig.hmftools.common.segmentation.Arm.Q;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

public class ChromosomeCopyNumbers
{
    private record Statistics(double mean, double median, double min, double max) {}

    private final Collection<PurpleCopyNumber> mRawCopyNumbers;
    private final ChrArmLocator mChrArmLocator;

    public ChromosomeCopyNumbers(Collection<PurpleCopyNumber> rawCopyNumbers, final ChrArmLocator chrArmLocator)
    {
        mRawCopyNumbers = rawCopyNumbers;
        mChrArmLocator = chrArmLocator;
    }

    public List<ChromosomeArmCopyNumber> data()
    {
        List<ChromosomeArmCopyNumber> result = new ArrayList<>();
        List<PurpleCopyNumber> currentArmCopyNumbers = new ArrayList<>();
        HumanChromosome currentChromosome = null;
        Arm currentArm = null;

        // The raw copy numbers are always in order (sorted by chromosome and location) and cover every arm of every chromosome
        for(PurpleCopyNumber purpleCopyNumber : mRawCopyNumbers)
        {
            final boolean switchChromosome = !Objects.equals(purpleCopyNumber.chr(), currentChromosome);
            Arm segmentArm = armForSegment(purpleCopyNumber);
            boolean switchArm = segmentArm != currentArm;
            if(switchChromosome || switchArm)
            {
                if(!currentArmCopyNumbers.isEmpty())
                {
                    Statistics statistics = calculateStats(currentArmCopyNumbers);
                    result.add(new ChromosomeArmCopyNumber(currentChromosome, currentArm, statistics.mean, statistics.median, statistics.min, statistics.max));
                }
                currentChromosome = purpleCopyNumber.chr();
                currentArm = segmentArm;
                currentArmCopyNumbers = new ArrayList<>();
            }
            currentArmCopyNumbers.add(purpleCopyNumber);
        }
        if(!currentArmCopyNumbers.isEmpty())
        {
            Statistics statistics = calculateStats(currentArmCopyNumbers);
            result.add(new ChromosomeArmCopyNumber(currentChromosome, currentArm, statistics.mean, statistics.median, statistics.min, statistics.max));
        }
        return result.stream().filter(ChromosomeArmCopyNumber::includeInReport).toList();
    }

    private Arm armForSegment(final PurpleCopyNumber purpleCopyNumber)
    {
        ChrArm startArm = mChrArmLocator.map(purpleCopyNumber.chromosome(), purpleCopyNumber.start());
        ChrArm endArm = mChrArmLocator.map(purpleCopyNumber.chromosome(), purpleCopyNumber.start());
        if(startArm.arm() != endArm.arm())
        {
            PPL_LOGGER.warn("Segment spans multiple chromosome arms: {}-{} on {}",
                    purpleCopyNumber.start(), purpleCopyNumber.end(), purpleCopyNumber.chromosome());
        }
        return startArm.arm() == P ? P : Q;
    }

    private Statistics calculateStats(final List<PurpleCopyNumber> rawCopyNumbers)
    {
        double totalLength = 0.0;
        double totalWeight = 0.0;
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        SortedSet<PurpleCopyNumber> rawSorted = new TreeSet<>(copyNumberComparator());
        for(PurpleCopyNumber purpleCopyNumber : rawCopyNumbers)
        {
            double length = purpleCopyNumber.length();
            totalLength += length;
            totalWeight += length * purpleCopyNumber.averageTumorCopyNumber();
            rawSorted.add(purpleCopyNumber);
            min = Math.min(min, purpleCopyNumber.averageTumorCopyNumber());
            max = Math.max(max, purpleCopyNumber.averageTumorCopyNumber());
        }
        double weightedMean = totalWeight / totalLength;
        int halfLength = (int) (totalLength / 2);
        double median = 0;
        int runningLength = 0;
        for(PurpleCopyNumber purpleCopyNumber : rawSorted)
        {
            runningLength += purpleCopyNumber.length();
            if(runningLength > halfLength)
            {
                median = purpleCopyNumber.averageTumorCopyNumber();
                break;
            }
        }
        return new Statistics(weightedMean, median, min, max);
    }

    private static Comparator<PurpleCopyNumber> copyNumberComparator()
    {
        return (o1, o2) ->
        {
            int byCN = Double.compare(o1.averageTumorCopyNumber(), o2.averageTumorCopyNumber());
            if(byCN == 0)
            {
                return Double.compare(o1.start(), o2.start());
            }
            return byCN;
        };
    }
}
