package com.hartwig.hmftools.cup.somatics;

import static com.hartwig.hmftools.common.sigs.SnvSigUtils.contextFromVariant;
import static com.hartwig.hmftools.common.sigs.SnvSigUtils.variantContext;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupConstants.AID_APOBEC_TRINUCLEOTIDE_CONTEXTS;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_OTHER;
import static com.hartwig.hmftools.cup.somatics.AidApobecStatus.ALL;
import static com.hartwig.hmftools.cup.somatics.AidApobecStatus.FALSE_ONLY;
import static com.hartwig.hmftools.cup.somatics.AidApobecStatus.TRUE_ONLY;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
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

    public static Matrix convertSomaticVariantsToPosFrequencies(
            final String sampleId, final List<SomaticVariant> variants, final Map<String,Integer> samplePosFreqIndex,
            final PositionFrequencies posFrequencies, AidApobecStatus aidApobecStatus)
    {
        posFrequencies.clear();
        extractPositionFrequencyCounts(variants, posFrequencies, aidApobecStatus);

        final Matrix matrix = new Matrix(posFrequencies.getCounts().length, 1);

        matrix.setCol(0, posFrequencies.getCounts());
        samplePosFreqIndex.put(sampleId, 0);

        return matrix;
    }

    public static void buildCancerPosFrequencies(
            final PositionFrequencies posFrequencies, final Matrix posFreqCounts, final Map<String,Integer> sampleIndexMap,
            final Map<String,List<SampleData>> refCancerSampleData, final String filename)
    {
        if(posFreqCounts == null)
            return;

        try
        {
            BufferedWriter writer = createBufferedWriter(filename, false);

            int maxSampleCount = posFrequencies.getMaxSampleCount();
            int bucketCount = posFrequencies.getBucketCount();

            final Map<String,double[]> cancerPosCounts = Maps.newHashMap();

            for(Map.Entry<String,List<SampleData>> entry : refCancerSampleData.entrySet())
            {
                final String cancerType = entry.getKey();

                if(cancerType.equals(CANCER_TYPE_OTHER))
                    continue;

                final double[] cancerCounts = new double[bucketCount];

                for(final SampleData sample : entry.getValue())
                {
                    final double[] sampleCounts = posFreqCounts.getCol(sampleIndexMap.get(sample.Id));

                    if(sampleCounts == null)
                        continue;

                    double sampleTotal = sumVector(sampleCounts);

                    double reductionFactor = sampleTotal > maxSampleCount ? maxSampleCount / sampleTotal : 1.0;

                    for(int b = 0; b < sampleCounts.length; ++b)
                    {
                        cancerCounts[b] += reductionFactor * sampleCounts[b];
                    }
                }

                cancerPosCounts.put(cancerType, cancerCounts);
            }

            final List<String> cancerTypes = cancerPosCounts.keySet().stream().collect(Collectors.toList());
            writer.write(cancerTypes.get(0));
            for(int i = 1; i < cancerTypes.size(); ++i)
            {
                writer.write(String.format(",%s", cancerTypes.get(i)));
            }

            writer.newLine();

            for(int b = 0; b < bucketCount; ++b)
            {
                writer.write(String.format("%.1f", cancerPosCounts.get(cancerTypes.get(0))[b]));

                for(int i = 1; i < cancerTypes.size(); ++i)
                {
                    writer.write(String.format(",%.1f", cancerPosCounts.get(cancerTypes.get(i))[b]));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample pos data output: {}", e.toString());
        }
    }

}
