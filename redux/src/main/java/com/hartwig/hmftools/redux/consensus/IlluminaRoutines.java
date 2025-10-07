package com.hartwig.hmftools.redux.consensus;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.UNSET_COUNT;
import static com.hartwig.hmftools.common.collect.Cluster.clusterCount;
import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.BASE_QUAL_MINIMUM;
import static com.hartwig.hmftools.common.sequencing.IlluminaBamUtils.getReadNameAttributes;
import static com.hartwig.hmftools.redux.consensus.BaseBuilder.INVALID_POSITION;
import static com.hartwig.hmftools.redux.consensus.BaseQualPair.NO_BASE;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.redux.BaseQualAdjustment;
import com.hartwig.hmftools.common.sequencing.IlluminaBamUtils;
import com.hartwig.hmftools.redux.duplicate.DuplicateGroup;

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

    private static BaseQualPair checkCommonBaseAndQual(
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
            int maxQual = max(locationQuals[0], locationQuals[1]);

            if(locationBases[0] == locationBases[1])
                return createBaseQualPair(locationBases[0], maxQual);

            int minQual = min(locationQuals[0], locationQuals[1]);
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

    private static final Set<Integer> OPTICAL_DUPLICATE_TILE_DIFFERENCE = Sets.newHashSet(0, 1, 999, 1000, 1001);
    private static final int OPTICAL_DUPLICATE_DISTANCE_THRESHOLD = 2_500;

    public static int calculatePCRClusterCount(final DuplicateGroup duplicateGroup)
    {
        List<IlluminaBamUtils.IlluminaReadNameAttributes> readNameAttributes = Lists.newArrayList();

        for(SAMRecord read : duplicateGroup.allReads())
        {
            IlluminaBamUtils.IlluminaReadNameAttributes attributes = getReadNameAttributes(read.getReadName());
            if(attributes == null)
                return UNSET_COUNT;

            readNameAttributes.add(attributes);
        }

        if(readNameAttributes.size() == 2)
        {
            IlluminaBamUtils.IlluminaReadNameAttributes readNameAttributes1 = readNameAttributes.get(0);
            IlluminaBamUtils.IlluminaReadNameAttributes readNameAttributes2 = readNameAttributes.get(1);
            int tileDifference = abs(readNameAttributes1.tileNumber() - readNameAttributes2.tileNumber());

            if(readNameAttributes1.laneKey().equals(readNameAttributes2.laneKey())
            && OPTICAL_DUPLICATE_TILE_DIFFERENCE.contains(tileDifference))
            {
                return 1;
            }

            return 2;
        }

        Map<String, List<IlluminaBamUtils.TileCoord>> tileCoordsByTile = Maps.newHashMap();

        for(IlluminaBamUtils.IlluminaReadNameAttributes attributes : readNameAttributes)
        {
            String tileKey = attributes.tileKey();
            IlluminaBamUtils.TileCoord tileCoord = attributes.tileCoord();
            tileCoordsByTile.computeIfAbsent(tileKey, key -> Lists.newArrayList());
            tileCoordsByTile.get(tileKey).add(tileCoord);
        }

        int pcrClusterCount = 0;
        for(List<IlluminaBamUtils.TileCoord> tileCoords : tileCoordsByTile.values())
        {
            if(tileCoords.size() == 1)
            {
                pcrClusterCount++;
                continue;
            }

            if(tileCoords.size() == 2)
            {
                IlluminaBamUtils.TileCoord tileCoord1 = tileCoords.get(0);
                IlluminaBamUtils.TileCoord tileCoord2 = tileCoords.get(1);

                if(tileCoord1.distance(tileCoord2) <= OPTICAL_DUPLICATE_DISTANCE_THRESHOLD)
                    pcrClusterCount++;
                else
                    pcrClusterCount += 2;

                continue;
            }

            pcrClusterCount += clusterCount(tileCoords, IlluminaBamUtils.TileCoord::distance, OPTICAL_DUPLICATE_DISTANCE_THRESHOLD);
        }

        return pcrClusterCount;
    }
}
