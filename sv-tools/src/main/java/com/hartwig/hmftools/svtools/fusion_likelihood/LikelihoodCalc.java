package com.hartwig.hmftools.svtools.fusion_likelihood;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.svtools.fusion_likelihood.CohortExpFusions.PRE_GENE_3P_DISTANCE;
import static com.hartwig.hmftools.svtools.fusion_likelihood.FusionLikelihood.FLC_LOGGER;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseRegion.hasAnyPhaseMatch;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseRegion.haveOverlap;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseRegion.regionsPhaseMatched;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

public class LikelihoodCalc
{
    public static Map<Integer,Long> calcOverlapBucketAreas(
            final List<Integer> bucketLengths, final List<RegionAllocator> regionAllocators,
            GeneRangeData lowerGene, GeneRangeData upperGene,
            GenePhaseRegion lowerRegion, GenePhaseRegion upperRegion, boolean isDel)
    {
        Map<Integer,Long> bucketOverlapCounts = Maps.newHashMap();

        if(lowerRegion.length() < 0 || upperRegion.length() < 0)
        {
            FLC_LOGGER.warn("negative region lengths");
            return bucketOverlapCounts;
        }

        for (int i = 0; i < bucketLengths.size() - 1; ++i)
        {
            int minBucketLen = bucketLengths.get(i);
            int maxBucketLen = bucketLengths.get(i + 1);

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
                if(regionAllocators != null)
                {
                    baseOverlapArea += regionAllocators.get(i).allocateBases(
                            lowerRegion.start(), lowerRegion.end(), upperRegion.start(), upperRegion.end(),
                            minBucketLen, maxBucketLen, true);
                }
                else
                {
                    baseOverlapArea = lowerRegion.length() * upperRegion.length();
                }
            }
            else
            {
                int lowerStart, lowerEnd, upperStart;

                int upperEnd = min(upperRegion.end(), lowerRegion.end() + maxBucketLen);

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

                int actUpperStart = max(upperStart, lowerStart + minBucketLen);
                int actUpperEnd = min(upperEnd, lowerEnd + maxBucketLen);

                if(regionAllocators != null)
                {
                    baseOverlapArea += regionAllocators.get(i).allocateBases(
                            lowerStart, lowerEnd, actUpperStart, actUpperEnd, minBucketLen, maxBucketLen, true);
                }
                else
                {
                    // per-base allocation without memory
                    for (int base = lowerStart; base <= lowerEnd; ++base)
                    {
                        int overlap = min(upperEnd, base + maxBucketLen) - max(upperStart, base + minBucketLen);

                        if (overlap > 0)
                            baseOverlapArea += overlap;

                        // is there any early exit here once overlap starts to be negative?
                    }
                }
            }

            setBucketLengthData(isDel ? lowerGene.getDelFusionBaseCounts() : lowerGene.getDupFusionBaseCounts(), i, baseOverlapArea);

            if(!upperGene.GeneData.GeneId.equals(lowerGene.GeneData.GeneId))
            {
                // avoid double-counting same gene fusions
                setBucketLengthData(isDel ? upperGene.getDelFusionBaseCounts() : upperGene.getDupFusionBaseCounts(), i, baseOverlapArea);
            }

            if(baseOverlapArea > 0)
                bucketOverlapCounts.put(i, baseOverlapArea);
        }

        return bucketOverlapCounts;
    }

    public static long calcGeneOverlapAreas(GeneRangeData geneUp, GeneRangeData geneDown,
            boolean isDel, int minBucketLen, int maxBucketLen, final RegionAllocator regionAllocator)
    {
        int upGeneTransStart = geneUp.getPhaseRegions().stream().mapToInt(x -> x.start()).min().orElse(0);
        int upGeneTransEnd = geneUp.getPhaseRegions().stream().mapToInt(x -> x.end()).max().orElse(0);
        int downGeneTransStart = geneDown.getPhaseRegions().stream().mapToInt(x -> x.start()).min().orElse(0);
        int downGeneTransEnd = geneDown.getPhaseRegions().stream().mapToInt(x -> x.end()).max().orElse(0);

        int upGeneStart = geneUp.GeneData.GeneStart;
        int upGeneEnd = geneUp.GeneData.GeneEnd;
        int downGeneStart = geneDown.GeneData.GeneStart;
        int downGeneEnd = geneDown.GeneData.GeneEnd;

        int strand = geneDown.GeneData.Strand;
        boolean sameGene = geneUp.GeneData.GeneId.equals(geneDown.GeneData.GeneId);

        if(!sameGene)
        {
            if (strand == 1)
                downGeneStart -= PRE_GENE_3P_DISTANCE;
            else
                downGeneEnd += PRE_GENE_3P_DISTANCE;
        }

        upGeneStart = min(upGeneStart, upGeneTransStart);
        upGeneEnd = max(upGeneEnd, upGeneTransEnd);
        downGeneStart = min(downGeneStart, downGeneTransStart);
        downGeneEnd = max(downGeneEnd, downGeneTransEnd);

        int lowerRegionStart, lowerRegionEnd, upperRegionStart, upperRegionEnd;

        if(sameGene)
        {
            upperRegionStart = lowerRegionStart = downGeneStart;
            upperRegionEnd = lowerRegionEnd = downGeneEnd;
        }
        else
        {
            // remove any overlaps by the DEL or DUP context
            if((isDel && strand == 1) || (!isDel && strand == -1))
            {
                lowerRegionStart = upGeneStart;
                lowerRegionEnd = upGeneEnd;
                upperRegionStart = max(downGeneStart, lowerRegionStart);
                upperRegionEnd = downGeneEnd;
            }
            else
            {
                upperRegionStart = upGeneStart;
                upperRegionEnd = upGeneEnd;
                lowerRegionStart = downGeneStart;
                lowerRegionEnd = min(downGeneEnd, upperRegionEnd);
            }
        }

        regionAllocator.reset();

        if(upperRegionEnd - upperRegionStart < 0 || lowerRegionEnd - lowerRegionStart < 0)
            return 0;

        // first check whether the bucket can link these 2 genes
        if(lowerRegionStart + minBucketLen >= upperRegionEnd || lowerRegionEnd + maxBucketLen <= upperRegionStart)
        {
            return 0;
        }

        long baseOverlapArea = 0;

        if(minBucketLen <= upperRegionStart - lowerRegionEnd && maxBucketLen >= upperRegionEnd - lowerRegionStart)
        {
            // no restriction on the overlap
            baseOverlapArea = regionAllocator.allocateBases(
                    lowerRegionStart, lowerRegionEnd, upperRegionStart, upperRegionEnd,
                    minBucketLen, maxBucketLen, true);
        }
        else
        {
            int lowerStart, lowerEnd, upperStart;

            int upperEnd = min(upperRegionEnd, lowerRegionEnd + maxBucketLen);

            if(lowerRegionStart + minBucketLen > upperRegionStart)
            {
                // min bucket length limits bases in the upper region
                lowerStart = lowerRegionStart;
                upperStart = lowerRegionStart + minBucketLen;

                if(lowerRegionEnd + minBucketLen > upperRegionEnd)
                {
                    lowerEnd = upperRegionEnd - minBucketLen;
                }
                else
                {
                    lowerEnd = lowerRegionEnd;
                }
            }
            else
            {
                // no min length restrictions, so just check if the max length is restrictive
                upperStart = upperRegionStart;
                lowerEnd = lowerRegionEnd;

                if(lowerRegionStart + maxBucketLen < upperRegionStart)
                {
                    lowerStart = upperRegionStart - maxBucketLen;
                }
                else
                {
                    lowerStart = lowerRegionStart;
                }
            }

            int actUpperStart = max(upperStart, lowerStart + minBucketLen);
            int actUpperEnd = min(upperEnd, lowerEnd + maxBucketLen);

            baseOverlapArea += regionAllocator.allocateBases(
                    lowerStart, lowerEnd, actUpperStart, actUpperEnd, minBucketLen, maxBucketLen, true);
        }

        return baseOverlapArea;
    }

    public static void setBucketLengthData(Map<Integer,Long> countsData, int bucketIndex, long newCounts)
    {
        if(newCounts <= 0)
            return;

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

    public static void reportGeneOverlaps(Map<String, List<GeneRangeData>> chrGeneDataMap)
    {
        // int blockSize = 1000;
        // int bucketMax = 250000; // to cover the longest arm
        // int[] overlapBuckets = new int[bucketMax];

        // exclude regions where one or both are not protein coding

        for (Map.Entry<String, List<GeneRangeData>> entry : chrGeneDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            List<GeneRangeData> geneRangeList = entry.getValue();

            if(geneRangeList.isEmpty())
                continue;

            FLC_LOGGER.info("calculating genic overlap for chromosome({}) geneCount({})", chromosome, geneRangeList.size());

            for(int s = 0; s <= 1; ++s)
            {
                int strand = (s == 0) ? 1 : -1;

                ChromosomeArm currentArm = geneRangeList.get(0).Arm;
                int totalGenicRegion = 0;

                for(int i = 0; i < geneRangeList.size(); ++i)
                {
                    GeneRangeData gene = geneRangeList.get(i);

                    if(gene.getPhaseRegions().isEmpty())
                        continue;

                    if(gene.GeneData.Strand != strand)
                        continue;

                    if(currentArm != gene.Arm)
                    {
                        // logArmStats(chromosome, currentArm, strand, overlapBuckets, totalGenicRegion, blockSize);

                        currentArm = gene.Arm;
                        totalGenicRegion = 0;
                    }

                    /*
                    int geneStart = gene.getPhaseRegions().get(0).start();
                    int geneEnd = gene.getPhaseRegions().get(gene.getPhaseRegions().size() - 1).end();

                    totalGenicRegion += (geneEnd - geneStart);

                    int lowerIndex = (int)floor(geneStart/(double)blockSize);
                    int upperIndex = (int)floor(geneEnd/(double)blockSize);

                    if(upperIndex < bucketMax)
                    {
                        for (int i = lowerIndex; i <= upperIndex; ++i)
                        {
                            ++overlapBuckets[i];
                        }
                    }
                    */

                    for(int j = i+1; j < geneRangeList.size(); ++j)
                    {
                        GeneRangeData gene2 = geneRangeList.get(j);

                        if (gene2.getPhaseRegions().isEmpty())
                            continue;

                        if (gene2.GeneData.Strand != strand)
                            continue;

                        if (!gene2.Arm.equals(gene.Arm))
                            break;

                        int totalOverlap = 0;

                        for (GenePhaseRegion region1 : gene.getPhaseRegions())
                        {
                            if(!region1.proteinCoding())
                                continue;

                            for (GenePhaseRegion region2 : gene2.getPhaseRegions())
                            {
                                if(!region2.proteinCoding())
                                    continue;

                                if (!haveOverlap(region1, region2, 0))
                                    continue;

                                int overlapStart = max(region1.start(), region2.start());
                                int overlapEnd = min(region1.end(), region2.end());

                                if(overlapEnd > overlapStart)
                                    totalOverlap += (overlapEnd - overlapStart);
                            }
                        }

                        if(totalOverlap > 0)
                        {
                            // Time,GeneId,GeneName,Strand,GeneStart,GeneEnd,GeneLength,RegionCount,CodingBases,ProteinCoding,
                            // OtherGeneId,OtherGeneName,OtherGeneStart,OtherGeneEnd,OtherGeneLength,OtherRegionCount,OtherCodingBases,TotalOverlap
                            FLC_LOGGER.info("GENE_OVERLAP: {},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                                    gene.GeneData.GeneId, gene.GeneData.GeneName, gene.GeneData.Strand,
                                    gene.GeneData.GeneStart, gene.GeneData.GeneEnd, gene.GeneData.GeneEnd - gene.GeneData.GeneStart,
                                    gene.getPhaseRegions().size(), gene.phasedRegionTotal(),
                                    gene2.GeneData.GeneId, gene2.GeneData.GeneName,
                                    gene2.GeneData.GeneStart, gene2.GeneData.GeneEnd, gene2.GeneData.GeneEnd - gene2.GeneData.GeneStart,
                                    gene2.getPhaseRegions().size(), gene2.phasedRegionTotal(), totalOverlap);
                        }
                    }
                }

                // logArmStats(chromosome, currentArm, strand, overlapBuckets, totalGenicRegion, blockSize);
            }
        }
    }

    private static void logArmStats(final String chromosome, final String arm, int strand,
            final int[] overlapBuckets, int totalGenicRegion, int blockSize)
    {
        int overlapRegions = 0;

        for(int i = 0; i < overlapBuckets.length;++i)
        {
            if(overlapBuckets[i] > 1)
                overlapRegions += blockSize;

            overlapBuckets[i] = 0;
        }

        FLC_LOGGER.info("chromosome({}) arm({}) strand({}) genicRegion({}) overlapRegion({}) percent({})",
                chromosome, arm, strand, totalGenicRegion, overlapRegions,
                String.format("%.1f", overlapRegions/(double)totalGenicRegion*100));
    }

    public static int calcNonProximateLikelihood(final GeneRangeData geneUp, final GeneRangeData geneDown)
    {
        // calculate the overlap area for 2 non-proximate genes
        int totalOverlap = 0;

        for (GenePhaseRegion regionUp : geneUp.getPhaseRegions())
        {
            // the downstream gene of the potential fusion cannot be non-coding
            for (GenePhaseRegion regionDown : geneDown.getPhaseRegions())
            {
                if(hasAnyPhaseMatch(regionUp, regionDown, false)
                        || regionsPhaseMatched(regionUp, regionDown))
                {
                    int regionOverlap = regionUp.length() * regionDown.length();
                    totalOverlap += regionOverlap;
                }
            }
        }

        return totalOverlap;
    }


}
