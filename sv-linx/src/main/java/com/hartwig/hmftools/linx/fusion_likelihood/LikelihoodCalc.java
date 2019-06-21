package com.hartwig.hmftools.linx.fusion_likelihood;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.hasAnyPhaseMatch;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.regionsPhaseMatched;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class LikelihoodCalc
{
    private static final Logger LOGGER = LogManager.getLogger(LikelihoodCalc.class);

    public static long calcNonProximateLikelihood(final GeneRangeData geneUp, final GeneRangeData geneDown)
    {
        // calculate the overlap area for 2 non-proximate genes
        long totalOverlap = 0;

        for (GenePhaseRegion regionUp : geneUp.getPhaseRegions())
        {
            // the downstream gene of the potential fusion cannot be non-coding
            for (GenePhaseRegion regionDown : geneDown.getPhaseRegions())
            {
                if(hasAnyPhaseMatch(regionUp, regionDown, false)
                || regionsPhaseMatched(regionUp, regionDown))
                {
                    long regionOverlap = regionUp.length() * regionDown.length();
                    totalOverlap += regionOverlap;
                }
            }
        }

        return totalOverlap;
    }

    public static Map<Integer, Long> calcOverlapBucketAreas(
            final List<Long> bucketLengths, GeneRangeData lowerGene, GeneRangeData upperGene,
            GenePhaseRegion lowerRegion, GenePhaseRegion upperRegion, boolean isDel)
    {
        Map<Integer, Long> bucketOverlapCounts = Maps.newHashMap();

        if(lowerRegion.length() < 0 || upperRegion.length() < 0)
        {
            LOGGER.warn("negative region lengths");
        }

        for (int i = 0; i < bucketLengths.size() - 1; ++i)
        {
            long minBucketLen = bucketLengths.get(i);
            long maxBucketLen = bucketLengths.get(i + 1);

            // first check whether the bucket can link these 2 genes
            if(lowerRegion.start() + minBucketLen >= upperRegion.end() || lowerRegion.end() + maxBucketLen <= upperRegion.start())
            {
                continue;
            }

            /* Example scenarios:
                - gene ranges 100-200 and 300-400
                - if min bucket length is less than 100 (ie the distance between the lower gene end and upper gene start) then no restriction
                - if max bucket length is greater than 300 (ie the distance from lower gene start to upper gene end) then no restriction
                - in this case the base overlap is 100 x 100, so 10K - the number of different ways a unique SV can fuse these 2 genes
                - if the min bucket length were increased say to 225, then the overlap ranges is limited
                - if the max bucket length were a limiting factor at say 275 then 100 + 275 only reaches 375
                - if the max bucket length = 150 then 100 + 150 doesn't reach the start of the upper gene, so the base overlap region must start from 150+
            */

            long baseOverlapArea = 0;

            if(minBucketLen <= upperRegion.start() - lowerRegion.end() && maxBucketLen >= upperRegion.end() - lowerRegion.start())
            {
                // no restriction on the overlap
                baseOverlapArea = lowerRegion.length() * upperRegion.length();
            }
            else
            {
                long lowerStart, lowerEnd, upperStart;

                long upperEnd = min(upperRegion.end(), lowerRegion.end() + maxBucketLen);

                if(lowerRegion.start() + minBucketLen > upperRegion.start())
                {
                    // min bucket length limits bases in the upper region
                    lowerStart = lowerRegion.start();
                    upperStart = lowerRegion.start() + minBucketLen;

                    if(lowerRegion.end() + minBucketLen > upperRegion.end())
                    {
                        lowerEnd = upperRegion.end() - minBucketLen;
                    }
                    else
                    {
                        lowerEnd = lowerRegion.end();
                    }
                }
                else
                {
                    // no min length restrictions, so just check if the max length is restrictive
                    upperStart = upperRegion.start();
                    lowerEnd = lowerRegion.end();

                    if(lowerRegion.start() + maxBucketLen < upperRegion.start())
                    {
                        lowerStart = upperRegion.start() - maxBucketLen;
                    }
                    else
                    {
                        lowerStart = lowerRegion.start();
                    }
                }

                for (long base = lowerStart; base <= lowerEnd; ++base)
                {
                    long overlap = min(upperEnd, base + maxBucketLen) - max(upperStart, base + minBucketLen);

                    if(overlap > 0)
                        baseOverlapArea += overlap;

                    // is there any early exit here once overlap starts to be negative?
                }
            }

            setBucketLengthData(isDel ? lowerGene.getDelFusionBaseCounts() : lowerGene.getDupFusionBaseCounts(), i, baseOverlapArea);

            if(!upperGene.GeneData.GeneId.equals(lowerGene.GeneData.GeneId))
            {
                // avoid double counting same gene fusions
                setBucketLengthData(isDel ? upperGene.getDelFusionBaseCounts() : upperGene.getDupFusionBaseCounts(), i, baseOverlapArea);
            }

            bucketOverlapCounts.put(i, baseOverlapArea);
        }

        return bucketOverlapCounts;

    }

    public static void setBucketLengthData(Map<Integer,Long> countsData, int bucketIndex, long newCounts)
    {
        // initialise the array if empty
        Long bucketCount = countsData.get(bucketIndex);

        if(bucketCount == null)
        {
            countsData.put(bucketIndex, newCounts);
        }
        else
        {
            countsData.put(bucketIndex, bucketCount + newCounts);
        }
    }

}
