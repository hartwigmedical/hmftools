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
import static com.hartwig.hmftools.cup.somatics.CopyNumberProfile.normaliseGenPosCountsByCopyNumber;
import static com.hartwig.hmftools.cup.somatics.SomaticsCommon.DEC_3_FORMAT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;
import com.hartwig.hmftools.common.utils.Matrix;
import com.hartwig.hmftools.common.utils.VectorUtils;
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

        final Matrix matrix = new Matrix(1, posFrequencies.getCounts().length);

        matrix.setRow(0, posFrequencies.getCounts());
        samplePosFreqIndex.put(sampleId, 0);

        return matrix;
    }

    public static Matrix buildCancerMatrix(
            final Matrix samplePosFreqCounts, final Map<String,Integer> sampleIndexMap,
            final List<String> cancerTypes, final Map<String,List<SampleData>> refCancerSampleData, int maxSampleCount)
    {
        Matrix cancerGenPosMatrix = new Matrix(samplePosFreqCounts.Rows, cancerTypes.size());
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
