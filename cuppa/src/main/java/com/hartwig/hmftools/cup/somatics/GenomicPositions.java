package com.hartwig.hmftools.cup.somatics;

import static com.hartwig.hmftools.common.sigs.SnvSigUtils.variantContext;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.AID_APOBEC_TRINUCLEOTIDE_CONTEXTS;
import static com.hartwig.hmftools.cup.somatics.AidApobecStatus.ALL;
import static com.hartwig.hmftools.cup.somatics.AidApobecStatus.FALSE_ONLY;
import static com.hartwig.hmftools.cup.somatics.AidApobecStatus.TRUE_ONLY;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.cup.common.SampleData;

public final class GenomicPositions
{
    public static void extractPositionFrequencyCounts(
            final List<SomaticVariant> variants, final PositionFrequencies posFrequencies, AidApobecStatus aidApobecStatus)
    {
        posFrequencies.clear();

        for(final SomaticVariant variant : variants)
        {
            if(variant.Type != VariantType.SNP)
                continue;

            // exclude male chromosome since is then unhelpful for multi-gender cancer types
            if(variant.Chromosome.equals("Y") || variant.Chromosome.equals("chrY"))
                continue;

            if(variant.TrinucleotideContext.contains("N"))
                continue;

            if(aidApobecStatus != ALL)
            {
                String bucketName = variantContext(variant.Ref, variant.Alt, variant.TrinucleotideContext);
                boolean isAA = AID_APOBEC_TRINUCLEOTIDE_CONTEXTS.contains(bucketName);

                if((aidApobecStatus == TRUE_ONLY && !isAA) || (aidApobecStatus == FALSE_ONLY && isAA))
                    continue;
            }

            if(!posFrequencies.isValidChromosome(variant.Chromosome))
            {
                CUP_LOGGER.warn("variant chr({}) position({}) cannot map to genomic position", variant.Chromosome, variant.Position);
                continue;
            }

            posFrequencies.addPosition(variant.Chromosome, variant.Position);
        }
    }

    @Deprecated
    public static Matrix convertSomaticVariantsToPosFrequencies(
            final String sampleId, final List<SomaticVariant> variants, final Map<String,Integer> samplePosFreqIndex,
            final PositionFrequencies posFrequencies, AidApobecStatus aidApobecStatus)
    {
        posFrequencies.clear();
        extractPositionFrequencyCounts(variants, posFrequencies, aidApobecStatus);

        final Matrix matrix = new Matrix(1, posFrequencies.getCounts().length);

        matrix.setRow(0, posFrequencies.getCounts());
        samplePosFreqIndex.put(sampleId, 0);

        return matrix;
    }

    public static void excludeChromosome(final Matrix matrix, final PositionFrequencies posFrequencies, final String chromosome)
    {
        int chromosomeStartIndex = posFrequencies.getBucketIndex(chromosome, 0);

        String nextChromosome = "";
        for(int i = 0; i < HumanChromosome.values().length; ++i)
        {
            if(HumanChromosome.values()[i].toString().equals(chromosome))
            {
                if(i < HumanChromosome.values().length - 1)
                    nextChromosome = HumanChromosome.values()[i + 1].toString();

                break;
            }
        }

        int chromosomeEndIndex;

        if(!nextChromosome.isEmpty())
        {
            chromosomeEndIndex = posFrequencies.getBucketIndex(nextChromosome, 0) - 1;
        }
        else
        {
            chromosomeEndIndex = posFrequencies.getBucketCount() - 1;
        }

        final double[][] data = matrix.getData();

        for(int b = chromosomeStartIndex; b <= chromosomeEndIndex; ++b)
        {
            for(int i = 0; i < matrix.Rows; ++i)
            {
                data[i][b] = 0;
            }
        }
    }

    @Deprecated
    public static Matrix buildCancerMatrix(
            final Matrix samplePosFreqCounts, final Map<String,Integer> sampleIndexMap,
            final List<String> cancerTypes, final Map<String,List<SampleData>> refCancerSampleData, int maxSampleCount)
    {
        // positions in the columns as per the sample matrix
        Matrix cancerGenPosMatrix = new Matrix(cancerTypes.size(), samplePosFreqCounts.Cols);
        final double[][] matrixData = cancerGenPosMatrix.getData();

        for(int i = 0; i < cancerTypes.size(); ++i)
        {
            String cancerType = cancerTypes.get(i);
            List<SampleData> samples = refCancerSampleData.get(cancerType);

            if(samples == null)
            {
                CUP_LOGGER.error("cancerType({}) missing ref samples", cancerType);
                return null;
            }

            for(final SampleData sample : samples)
            {
                if(!sampleIndexMap.containsKey(sample.Id))
                {
                    CUP_LOGGER.error("gen-pos sample index missing sample({})", sample.Id);
                    return null;
                }

                int sampleIndex = sampleIndexMap.get(sample.Id);
                double[] sampleCounts = samplePosFreqCounts.getRow(sampleIndex);

                if(sampleCounts == null)
                    continue;

                double sampleTotal = sumVector(sampleCounts);

                double reductionFactor = sampleTotal > maxSampleCount ? maxSampleCount / sampleTotal : 1.0;

                for(int b = 0; b < sampleCounts.length; ++b)
                {
                    matrixData[i][b] += reductionFactor * sampleCounts[b];
                }
            }
        }

        return cancerGenPosMatrix;
    }
}
