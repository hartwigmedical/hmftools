package com.hartwig.hmftools.cobalt.segmentation;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.segmentation.PiecewiseConstantFit;
import com.hartwig.hmftools.common.segmentation.Segmenter;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.FastMath;

class RatioSegmenter
{
    private final ListMultimap<ChrArm, CobaltRatio> mChromosomeToRatios = ArrayListMultimap.create();
    private final double mGamma;

    public RatioSegmenter(List<CobaltRatio> ratios, ChrArmLocator chrArmLocator, double gamma)
    {
        ratios.forEach(ratio -> mChromosomeToRatios.put(chrArmLocator.map(ratio), ratio));
        mGamma = gamma;
    }

    public Map<ChrArm, PiecewiseConstantFit> getSegmentation(ExecutorService executor)
            throws Exception
    {
        final HashMap<ChrArm, PiecewiseConstantFit> result = new HashMap<>();
        List<Future<Pair<ChrArm, PiecewiseConstantFit>>> futures = Lists.newArrayList();

        mChromosomeToRatios.keySet().forEach(chr ->
        {
            Callable<Pair<ChrArm, PiecewiseConstantFit>> task = () -> calculateSegmentsForChromosome(mChromosomeToRatios.get(chr), chr);
            futures.add(executor.submit(task));
        });

        for(Future<Pair<ChrArm, PiecewiseConstantFit>> future : futures)
        {
            Pair<ChrArm, PiecewiseConstantFit> chrPCF = future.get();
            result.put(chrPCF.getLeft(), chrPCF.getRight());
        }
        return result;
    }

    private Pair<ChrArm, PiecewiseConstantFit> calculateSegmentsForChromosome(List<CobaltRatio> ratios, ChrArm arm)
    {
        List<CobaltRatio> nonZeroRatios = ratios.stream().filter(cobaltRatio -> cobaltRatio.tumorGCRatio() >= 0.0).toList();
        double[] d = new double[nonZeroRatios.size()];
        for(int i = 0; i < nonZeroRatios.size(); i++)
        {
            CobaltRatio ratio = nonZeroRatios.get(i);
            d[i] = FastMath.log(2, ratio.tumorGCRatio());
        }
        Segmenter segmenter = new Segmenter(d, mGamma, true);
        PiecewiseConstantFit fit = segmenter.pcf();
        return Pair.of(arm, fit);
    }
}
