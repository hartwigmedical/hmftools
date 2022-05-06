package com.hartwig.hmftools.cup.somatics;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.formSamplePath;
import static com.hartwig.hmftools.cup.common.CupConstants.GEN_POS_CN_ADJUST_MAX;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.cup.ref.RefDataConfig;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public final class CopyNumberProfile
{
    public static double[] normaliseGenPosCountsByCopyNumber(double samplePloidy, final double[] counts, final double[] copyNumber)
    {
        final double[] adjustedCounts = new double[counts.length];

        for(int i = 0; i < counts.length; ++i)
        {
            double adjust = samplePloidy / max(0.01, copyNumber[i]);
            adjustedCounts[i] = counts[i] * min(adjust, GEN_POS_CN_ADJUST_MAX);
        }

        return adjustedCounts;
    }

    public static void buildCopyNumberProfile(
            final List<String> sampleIds, final RefDataConfig config, final Matrix copyNumberMatrix,
            final Map<String,Integer> sampleIndexMap, final PositionFrequencies posFrequencies)
    {
        CUP_LOGGER.debug("building copy-number profile reference data");

        for(String sampleId : sampleIds)
        {
            int sampleIndex = sampleIndexMap.get(sampleId);

            double[] cnProfile = null;

            if(config.DbAccess != null)
            {
                cnProfile = extractCopyNumberProfile(sampleId, config.DbAccess, null, posFrequencies);
            }
            else
            {
                final String purpleCopyNumberFile = formSamplePath(config.SampleCopyNumberFile, sampleId);
                cnProfile = extractCopyNumberProfile(sampleId, null, purpleCopyNumberFile, posFrequencies);
            }

            copyNumberMatrix.setCol(sampleIndex, cnProfile);
        }
    }

    public static double[] extractCopyNumberProfile(
            final String sampleId, final DatabaseAccess dbAccess, final String copyNumberFile, final PositionFrequencies posFrequencies)
    {
        try
        {
            final List<PurpleCopyNumber> copyNumbers = dbAccess != null ?
                    dbAccess.readCopynumbers(sampleId) : PurpleCopyNumberFile.read(copyNumberFile);

            return extractCopyNumberProfile(copyNumbers, posFrequencies);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error(" failed to retrieve CN data for sample({})", sampleId, e.toString());
            return new double[posFrequencies.getBucketCount()];
        }
    }

    public static double[] extractCopyNumberProfile(final List<PurpleCopyNumber> copyNumbers, final PositionFrequencies posFrequencies)
    {
        int totalBuckets = posFrequencies.getBucketCount();
        int bucketSize = posFrequencies.getBucketSize();
        double[] cnProfile = new double[posFrequencies.getBucketCount()];

        int currentBucketIndex = 0;
        double currentBucketCnTotal = 0;
        String currentChromosome = "";
        int chrStartBucket = 0;
        int chrEndBucket = 0;
        int bucketPosStart = 0;
        int bucketPosEnd = 0;

        for(PurpleCopyNumber copyNumber : copyNumbers)
        {
            String chromosome = copyNumber.chromosome();
            int cnRegionStart = copyNumber.start();
            int cnRegionEnd = copyNumber.end();
            double regionCN = copyNumber.averageTumorCopyNumber();
            boolean runsToTelomere = copyNumber.segmentEndSupport() == SegmentSupport.TELOMERE;

            // scenarios:
            // new chromosome purple record - fill in buckets up to the one that straddles past the end of the record
            // purple record runs out to telomere, fill in same CN for all remaining buckets on this chromosome
            // next purple record - fill in proportion of previous bucket, then continue
            if(!currentChromosome.equals(chromosome))
            {
                currentChromosome = chromosome;
                chrStartBucket = posFrequencies.getBucketIndex(chromosome, 1);
                currentBucketIndex = chrStartBucket;

                String nextChr = nextChromosome(chromosome);

                if(nextChr != null && posFrequencies.isValidChromosome(nextChr))
                    chrEndBucket = posFrequencies.getBucketIndex(nextChr, 1) - 1;
                else
                    chrEndBucket = totalBuckets - 1;

                bucketPosStart = 1;
                bucketPosEnd = bucketPosStart + bucketSize - 1;

                currentBucketCnTotal = 0;
            }
            else
            {
                // check if new purple record now goes past the current bucket
                if(cnRegionEnd >= bucketPosEnd)
                {
                    int remainingBucketBases = bucketPosEnd - cnRegionStart + 1;
                    currentBucketCnTotal += remainingBucketBases * regionCN;
                    double calcCn = currentBucketCnTotal / bucketSize;
                    cnProfile[currentBucketIndex] = calcCn;

                    ++currentBucketIndex;
                    bucketPosStart += bucketSize;
                    bucketPosEnd = bucketPosStart + bucketSize - 1;
                    currentBucketCnTotal = 0;
                }
            }

            if(runsToTelomere)
            {
                for(int bucket = currentBucketIndex; bucket <= chrEndBucket; ++bucket)
                {
                    cnProfile[bucket] = regionCN;
                }
            }
            else if(cnRegionEnd <= bucketPosEnd)
            {
                // purple region finishes before bucket
                currentBucketCnTotal += (cnRegionEnd - max(bucketPosStart, cnRegionStart) + 1) * regionCN;
            }
            else
            {
                // purple region extends beyond this bucket
                int cnRegionBucketEnd = posFrequencies.getBucketIndex(chromosome, cnRegionEnd);

                for(; currentBucketIndex < cnRegionBucketEnd; ++currentBucketIndex)
                {
                    cnProfile[currentBucketIndex] = regionCN;
                    bucketPosStart += bucketSize;
                    bucketPosEnd = bucketPosStart + bucketSize - 1;
                }

                // last bucket is only partially covered
                int newBucketBases = cnRegionEnd - bucketPosStart;
                currentBucketCnTotal += newBucketBases * regionCN;
            }
        }

        return cnProfile;
    }

    private static String nextChromosome(final String chromosome)
    {
        boolean getNext = false;

        for(HumanChromosome humanChromosome : HumanChromosome.values())
        {
            final String chr = humanChromosome.toString();

            if(getNext)
                return chr;

            if(chr.equals(chromosome))
                getNext = true;
        }

        return null;
    }
}
