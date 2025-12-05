package com.hartwig.hmftools.common.segmentation;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.segmentation.copynumber.FixedPenalty;
import com.hartwig.hmftools.common.segmentation.copynumber.GammaPenaltyCalculator;
import com.hartwig.hmftools.common.segmentation.copynumber.PenaltyCalculator;
import com.hartwig.hmftools.common.segmentation.copynumber.PiecewiseConstantFit;
import com.hartwig.hmftools.common.segmentation.copynumber.Segmenter;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public abstract class PerArmSegmenter<T extends GenomePosition>
{
    private static final Logger SG_LOGGER = LogManager.getLogger(PerArmSegmenter.class);
    private final ListMultimap<ChrArm, T> ArmToRatios = ArrayListMultimap.create();
    private final Map<ChrArm, DataForSegmentation> mDataByArm = new HashMap<>();
    private final PenaltyCalculator mPenaltyCalculator;

    protected PerArmSegmenter(ListMultimap<Chromosome, T> ratios, ChrArmLocator chrArmLocator, double gamma)
    {
        ratios.keySet().forEach(chromosome ->
        {
            List<T> ratiosForChromosome = ratios.get(chromosome);
            ratiosForChromosome.forEach(cobaltRatio ->
            {
                if(value(cobaltRatio) >= 0.0)
                {
                    ArmToRatios.put(chrArmLocator.map(cobaltRatio), cobaltRatio);
                }
            });
        });
        ArmToRatios.keySet().forEach(chrArm -> mDataByArm.put(chrArm, buildSegmentationData(ArmToRatios.get(chrArm))));
        int totalCount = mDataByArm.values().stream().mapToInt(DataForSegmentation::count).sum();
        int position = 0;
        if(totalCount < 100_000)
        {
            double[] allRatios = new double[totalCount];
            for(ChrArm chrArm : mDataByArm.keySet())
            {
                DataForSegmentation data = mDataByArm.get(chrArm);
                System.arraycopy(data.valuesForSegmentation(), 0, allRatios, position, data.count());
                position += data.count();
            }
            SG_LOGGER.info("Using uniform segmentation penalty, number of ratios: {}", allRatios.length);
            GammaPenaltyCalculator oneOffCalculation = new GammaPenaltyCalculator(gamma, true);
            final double penalty = oneOffCalculation.getPenalty(allRatios);
            mPenaltyCalculator = new FixedPenalty(penalty);
            SG_LOGGER.info("Uniform segmentation penalty: {}", penalty);
        }
        else
        {
            mPenaltyCalculator = new GammaPenaltyCalculator(gamma, true);
        }
    }

    public Map<ChrArm, ChromosomeArmSegments<T>> getSegmentation(ExecutorService executor) throws Exception
    {
        final HashMap<ChrArm, ChromosomeArmSegments<T>> result = new HashMap<>();
        List<Future<Pair<ChrArm, ChromosomeArmSegments<T>>>> futures = Lists.newArrayList();

        ArmToRatios.keySet().forEach(chr ->
        {
            Callable<Pair<ChrArm, ChromosomeArmSegments<T>>> task = () -> calculateSegmentsForChromosome(ArmToRatios.get(chr), chr);
            futures.add(executor.submit(task));
        });

        for(Future<Pair<ChrArm, ChromosomeArmSegments<T>>> future : futures)
        {
            Pair<ChrArm, ChromosomeArmSegments<T>> chrPCF = future.get();
            if(chrPCF != null)
            {
                result.put(chrPCF.getLeft(), chrPCF.getRight());
            }
        }
        return result;
    }

    public abstract double value(T ratio);

    public abstract DataForSegmentation buildSegmentationData(List<T> ratios);

    public abstract boolean isWindowed();

    private Pair<ChrArm, ChromosomeArmSegments<T>> calculateSegmentsForChromosome(List<T> ratios, ChrArm arm)
    {
        DataForSegmentation dataForArm = mDataByArm.get(arm);
        if(dataForArm.isEmpty())
        {
            return null;
        }
        Segmenter segmenter = new Segmenter(dataForArm.valuesForSegmentation(), mPenaltyCalculator);
        PiecewiseConstantFit fit = segmenter.pcf();
        ChromosomeArmSegments<T> segments;
        if(isWindowed())
        {
            segments = new WindowSegments<>(arm, ratios, dataForArm.rawValues(), fit);
        }
        else
        {
            segments = new AbsoluteSegments<>(arm, ratios, dataForArm.rawValues(), fit);
        }
        return Pair.of(arm, segments);
    }
}
