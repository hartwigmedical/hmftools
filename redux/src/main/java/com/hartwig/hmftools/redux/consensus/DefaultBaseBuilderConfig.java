package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.ceil;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.BASE_QUAL_MINIMUM;
import static com.hartwig.hmftools.common.qual.BaseQualAdjustment.adjustBaseQual;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.NO_BASE;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

public class DefaultBaseBuilderConfig extends BaseBuilderConfig
{
    public DefaultBaseBuilderConfig(final RefGenome refGenome, final ConsensusStatistics consensusStats)
    {
        super(refGenome, true, true, consensusStats);
    }

    protected DefaultBaseBuilderConfig(final RefGenome refGenome, boolean pairedReads, boolean useSimpleNoMismatchLogic,
            final ConsensusStatistics consensusStats)
    {
        super(refGenome, pairedReads, useSimpleNoMismatchLogic, consensusStats);
    }

    @Override
    public byte[] determineBaseAndQual(final boolean isDualStrand, final boolean[] isFirstInPair, final byte[] locationBases,
            final byte[] locationQuals, final String chromosome, final int position)
    {
        byte[] consensusBaseAndQual;
        if(isDualStrand && position != INVALID_POSITION)
        {
            // split the reads into 2 consensus reads and then compare
            consensusBaseAndQual = determineDualStrandBaseAndQual(isFirstInPair, locationBases, locationQuals, chromosome, position, BASE_QUAL_MINIMUM);
        }
        else
        {
            consensusBaseAndQual = determineBaseAndQual(locationBases, locationQuals, chromosome, position, BASE_QUAL_MINIMUM);
        }

        consensusBaseAndQual[1] = adjustBaseQual(consensusBaseAndQual[1]);
        return consensusBaseAndQual;
    }

    private byte[] determineDualStrandBaseAndQual(
            final boolean[] isFirstInPair, final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position, byte minBaseQual)
    {
        // TODO: Review this for same bug.
        int readCount = isFirstInPair.length;

        int firstInPairCount = 0;
        int secondInPairCount = 0;

        for(int i = 0; i < isFirstInPair.length; ++i)
        {
            if(isFirstInPair[i])
            {
                ++firstInPairCount;
            }
            else
            {
                ++secondInPairCount;
            }
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

        byte[] firstBaseAndQual = determineBaseAndQual(locationBasesFirst, locationQualsFirst, chromosome, position, minBaseQual);
        byte[] secondBaseAndQual = determineBaseAndQual(locationBasesSecond, locationQualsSecond, chromosome, position, minBaseQual);

        if(firstBaseAndQual[0] == NO_BASE)
        {
            return secondBaseAndQual;
        }

        if(secondBaseAndQual[0] == NO_BASE)
        {
            return firstBaseAndQual;
        }

        if(firstBaseAndQual[0] == secondBaseAndQual[0])
        {
            byte qual = (byte) max(firstBaseAndQual[1], secondBaseAndQual[1]);
            return new byte[] { firstBaseAndQual[0], qual };
        }

        mConsensusStats.registerDualStrandMismatchReadGroup(readCount);

        byte refBase = mRefGenome.getRefBase(chromosome, position);
        boolean firstIsRef = firstBaseAndQual[0] == refBase;
        boolean secondIsRef = secondBaseAndQual[0] == refBase;

        if(!firstIsRef && !secondIsRef)
        {
            byte maxBase;
            int maxQual;
            int differingQual;
            if(firstBaseAndQual[1] >= secondBaseAndQual[1])
            {
                maxBase = firstBaseAndQual[0];
                maxQual = firstBaseAndQual[1];
                differingQual = secondBaseAndQual[1];
            }
            else
            {
                maxBase = secondBaseAndQual[0];
                maxQual = secondBaseAndQual[1];
                differingQual = firstBaseAndQual[1];
            }

            byte qual = (byte) max(0, maxQual - differingQual);
            return new byte[] { maxBase, qual };
        }

        int refQual;
        int differingQual;
        if(firstIsRef)
        {
            refQual = firstBaseAndQual[1];
            differingQual = secondBaseAndQual[1];
        }
        else
        {
            refQual = secondBaseAndQual[1];
            differingQual = firstBaseAndQual[1];
        }

        byte qual = (byte) max(BASE_QUAL_MINIMUM, refQual - differingQual);
        return new byte[] { refBase, qual };
    }

    protected byte[] determineBaseAndQual(final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position, byte minBaseQual)
    {
        if(locationBases.length == 1)
        {
            // early exit for dual strand with a single read on one side - a very common scenario
            return new byte[] { locationBases[0], locationQuals[0] };
        }

        List<Byte> distinctBases = Lists.newArrayListWithCapacity(4);
        List<Integer> qualTotals = Lists.newArrayListWithCapacity(4);
        List<Integer> maxQuals = Lists.newArrayListWithCapacity(4);
        for(int i = 0; i < locationBases.length; ++i)
        {
            if(locationBases[i] == NO_BASE)
            {
                continue;
            }

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
                qualTotals.add((int) locationQuals[i]);
                maxQuals.add((int) locationQuals[i]);
            }
        }

        if(distinctBases.isEmpty())
        {
            return new byte[] { NO_BASE, 0 };
        }

        // chromosome will be null for unmapped reads
        Byte refBase = chromosome != null && position != INVALID_POSITION ? mRefGenome.getRefBase(chromosome, position) : null;
        byte maxBase = distinctBases.get(0);
        int maxQualTotal = qualTotals.get(0);
        for(int i = 1; i < distinctBases.size(); ++i)
        {
            if(qualTotals.get(i) > maxQualTotal)
            {
                maxQualTotal = qualTotals.get(i);
                maxBase = distinctBases.get(i);
            }
            else if(qualTotals.get(i) == maxQualTotal && refBase != null && refBase.equals(distinctBases.get(i)))
            {
                maxBase = refBase;
            }
        }

        // collect base quals matching the selected base to find the median
        List<Integer> selectBaseQuals = Lists.newArrayList();

        for(int i = 0; i < locationBases.length; ++i)
        {
            if(locationBases[i] == maxBase)
            {
                selectBaseQuals.add((int) locationQuals[i]);
            }
        }

        Collections.sort(selectBaseQuals);

        int medianBaseQualIndex = selectBaseQuals.size() / 2;
        int medianBaseQual = selectBaseQuals.get(medianBaseQualIndex);

        int differingQual = 0;

        for(int i = 0; i < distinctBases.size(); ++i)
        {
            if(distinctBases.get(i) != maxBase)
            {
                differingQual += qualTotals.get(i);
            }
        }

        double calcQual = (double) medianBaseQual * max(minBaseQual, maxQualTotal - differingQual) / maxQualTotal;

        return new byte[] { maxBase, (byte) round(ceil(calcQual)) };
    }
}
