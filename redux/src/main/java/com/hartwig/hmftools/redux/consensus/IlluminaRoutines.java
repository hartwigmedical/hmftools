package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.BASE_QUAL_MINIMUM;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;
import static com.hartwig.hmftools.redux.consensus.BaseQualPair.NO_BASE;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;

import htsjdk.samtools.SAMRecord;

public final class IlluminaRoutines
{
    private static BaseQualPair createBaseQualPair(final byte base, final double qual)
    {
        // scale calculated quals to standard values
        return new BaseQualPair(base, BaseQualAdjustment.adjustBaseQual(qual));
    }

    public static BaseQualPair determineDualStrandBaseAndQual(
            final boolean[] isFirstInPair, final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position,
            final RefGenome refGenome)
    {
        int firstInPairCount = 0;
        int secondInPairCount = 0;

        for(int i = 0; i < isFirstInPair.length; ++i)
        {
            if(isFirstInPair[i])
                ++firstInPairCount;
            else
                ++secondInPairCount;
        }

        byte[] locationBasesFirst = new byte[firstInPairCount];
        byte[] locationQualsFirst = new byte[firstInPairCount];

        byte[] locationBasesSecond = new byte[secondInPairCount];
        byte[] locationQualsSecond = new byte[secondInPairCount];

        int firstIndex = 0;
        int secondIndex = 0;

        for(int i = 0; i < locationBases.length; ++i)
        {
            if(isFirstInPair[i])
            {
                locationBasesFirst[firstIndex] = locationBases[i];
                locationQualsFirst[firstIndex] = locationQuals[i];
                ++firstIndex;
            }
            else
            {
                locationBasesSecond[secondIndex] = locationBases[i];
                locationQualsSecond[secondIndex] = locationQuals[i];
                ++secondIndex;
            }
        }

        BaseQualPair firstBaseAndQual = determineBaseAndQual(locationBasesFirst, locationQualsFirst, chromosome, position, refGenome);
        BaseQualPair secondBaseAndQual = determineBaseAndQual(locationBasesSecond, locationQualsSecond, chromosome, position, refGenome);

        if(!firstBaseAndQual.isValid())
            return secondBaseAndQual;

        if(!secondBaseAndQual.isValid())
            return firstBaseAndQual;

        if(firstBaseAndQual.Base == secondBaseAndQual.Base)
        {
            return createBaseQualPair(firstBaseAndQual.Base, max(firstBaseAndQual.Qual, secondBaseAndQual.Qual));
        }

        // TODO: decide whether to track this info or now
        // mConsensusStats.registerDualStrandMismatchReadGroup(readCount);

        byte refBase = refGenome.getRefBase(chromosome, position);
        boolean firstIsRef = firstBaseAndQual.Base == refBase;
        boolean secondIsRef = secondBaseAndQual.Base == refBase;

        if(!firstIsRef && !secondIsRef)
        {
            byte maxBase;
            int maxQual;
            int differingQual;
            if(firstBaseAndQual.Qual >= secondBaseAndQual.Qual)
            {
                maxBase = firstBaseAndQual.Base;
                maxQual = firstBaseAndQual.Qual;
                differingQual = secondBaseAndQual.Qual;
            }
            else
            {
                maxBase = secondBaseAndQual.Base;
                maxQual = secondBaseAndQual.Qual;
                differingQual = firstBaseAndQual.Qual;
            }

            byte qual = (byte) max(0, maxQual - differingQual);
            return createBaseQualPair(maxBase, qual);
        }

        int refQual;
        int differingQual;
        if(firstIsRef)
        {
            refQual = firstBaseAndQual.Qual;
            differingQual = secondBaseAndQual.Qual;
        }
        else
        {
            refQual = secondBaseAndQual.Qual;
            differingQual = firstBaseAndQual.Qual;
        }

        byte qual = (byte) max(BASE_QUAL_MINIMUM, refQual - differingQual);
        return createBaseQualPair(refBase, qual);
    }

    public static BaseQualPair checkCommonBaseAndQual(
            final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position, final RefGenome refGenome)
    {
        if(locationBases.length == 1)
        {
            // early exit for dual strand with a single read on one side - a very common scenario
            return new BaseQualPair(locationBases[0], locationQuals[0]);
        }

        // most common scenario is 2 reads with differing bases
        if(locationBases.length == 2)
        {
            int minQual = min(locationQuals[0], locationQuals[1]);
            int maxQual = max(locationQuals[0], locationQuals[1]);
            double calcQual = (double)maxQual * max(BASE_QUAL_MINIMUM, maxQual - minQual) / maxQual;

            if(locationQuals[0] > locationQuals[1])
                return createBaseQualPair(locationBases[0], calcQual);

            if(locationQuals[1] > locationQuals[0])
                return createBaseQualPair(locationBases[1], calcQual);

            // select whichever matches the ref otherwise the first
            if(chromosome != null && position != INVALID_POSITION)
            {
                byte refBase = refGenome.getRefBase(chromosome, position);

                if(locationBases[0] == refBase)
                    return createBaseQualPair(locationBases[0], calcQual);

                if(locationBases[1] == refBase)
                    return createBaseQualPair(locationBases[1], calcQual);
            }
            else
            {
                return createBaseQualPair(locationBases[0], calcQual);
            }
        }

        return BaseQualPair.INVALID;
    }

    public static BaseQualPair determineBaseAndQual(
            final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position, final RefGenome refGenome)
    {
        BaseQualPair baseQualPair = checkCommonBaseAndQual(locationBases, locationQuals, chromosome, position, refGenome);

        if(baseQualPair != BaseQualPair.INVALID)
            return baseQualPair;

        List<Byte> distinctBases = Lists.newArrayListWithCapacity(4);
        List<Integer> qualTotals = Lists.newArrayListWithCapacity(4);
        List<Integer> maxQuals = Lists.newArrayListWithCapacity(4);

        for(int i = 0; i < locationBases.length; ++i)
        {
            if(locationBases[i] == NO_BASE)
                continue;

            boolean found = false;

            for(int j = 0; j < distinctBases.size(); ++j)
            {
                if(distinctBases.get(j) == locationBases[i])
                {
                    int qualTotal = qualTotals.get(j) + locationQuals[i];
                    qualTotals.set(j, qualTotal);
                    maxQuals.set(j, max(maxQuals.get(j), locationQuals[i]));
                    found = true;
                    break;
                }
            }

            if(!found)
            {
                distinctBases.add(locationBases[i]);
                qualTotals.add((int)locationQuals[i]);
                maxQuals.add((int)locationQuals[i]);
            }
        }

        if(distinctBases.isEmpty())
        {
            return BaseQualPair.INVALID;
        }

        byte maxBase = distinctBases.get(0);
        boolean maxIsRef = false;
        int maxQualTotal = qualTotals.get(0);

        for(int i = 1; i < distinctBases.size(); ++i)
        {
            if(qualTotals.get(i) > maxQualTotal)
            {
                maxQualTotal = qualTotals.get(i);
                maxBase = distinctBases.get(i);
            }
            else if(chromosome != null && qualTotals.get(i) >= maxQualTotal && !maxIsRef && position != INVALID_POSITION)
            {
                // chromosome will be null for unmapped reads
                byte refBase = refGenome.getRefBase(chromosome, position);

                if(maxBase == refBase)
                {
                    maxIsRef = true;
                }
                else if(distinctBases.get(i) == refBase)
                {
                    maxQualTotal = qualTotals.get(i);
                    maxBase = distinctBases.get(i);
                    maxIsRef = true;
                }
            }
        }

        // collect base quals matching the selected base to find the median
        List<Integer> selectBaseQuals = Lists.newArrayList();

        for(int i = 0; i < locationBases.length; ++i)
        {
            if(locationBases[i] == maxBase)
                selectBaseQuals.add((int)locationQuals[i]);
        }

        Collections.sort(selectBaseQuals);

        int medianBaseQualIndex = selectBaseQuals.size() / 2;
        int medianBaseQual = selectBaseQuals.get(medianBaseQualIndex);

        int differingQual = 0;

        for(int i = 0; i < distinctBases.size(); ++i)
        {
            if(distinctBases.get(i) != maxBase)
                differingQual += qualTotals.get(i);
        }

        double calcQual = (double)medianBaseQual * max(BASE_QUAL_MINIMUM, maxQualTotal - differingQual) / maxQualTotal;

        return createBaseQualPair(maxBase, calcQual);
    }

    public static boolean isDualStrandAndIsFirstInPair(final List<SAMRecord> reads, final boolean[] isFirstInPairOut)
    {
        if(!reads.get(0).getReadPairedFlag())
            return false;

        boolean isDualStrand = false;

        for(int i = 0; i < reads.size(); ++i)
        {
            isFirstInPairOut[i] = reads.get(i).getFirstOfPairFlag();

            if(i > 0)
                isDualStrand |= isFirstInPairOut[0] != isFirstInPairOut[i];
        }

        return isDualStrand;
    }
}
