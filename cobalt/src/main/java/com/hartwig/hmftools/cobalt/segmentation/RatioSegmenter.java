package com.hartwig.hmftools.cobalt.segmentation;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.segmentation.PiecewiseConstantFit;
import com.hartwig.hmftools.common.segmentation.Segmenter;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.FastMath;

public abstract class RatioSegmenter
{
    private final ListMultimap<ChrArm, CobaltRatio> ArmToRatios = ArrayListMultimap.create();
    private final ConcurrentHashMap<ChrArm, Integer> ArmToIndexOfStart = new ConcurrentHashMap<>();
    private final double mGamma;

    public static void writeTumorSegments(ListMultimap<Chromosome, CobaltRatio> ratios,
            double gamma,
            RefGenomeVersion genomeVersion,
            ExecutorService executor,
            String outputPath) throws Exception
    {
        writeSegments(ratios, gamma, genomeVersion, executor, outputPath, true);
    }

    public static void writeReferenceSegments(ListMultimap<Chromosome, CobaltRatio> ratios,
            double gamma,
            RefGenomeVersion genomeVersion,
            ExecutorService executor,
            String outputPath) throws Exception
    {
        writeSegments(ratios, gamma, genomeVersion, executor, outputPath, false);
    }

    private static void writeSegments(
            ListMultimap<Chromosome, CobaltRatio> ratios,
            double gamma,
            RefGenomeVersion genomeVersion,
            ExecutorService executor,
            String outputPath,
            boolean isForTumor) throws Exception
    {
        ChrArmLocator locator = ChrArmLocator.defaultLocator(genomeVersion);
        RatioSegmenter segmenter =
                isForTumor ? new TumorRatioSegmenter(ratios, locator, gamma) : new ReferenceRatioSegmenter(ratios, locator, gamma);
        Map<ChrArm, PiecewiseConstantFit> segmentsByChrArm = segmenter.getSegmentation(executor);
        SegmentsFile.write(segmentsByChrArm, segmenter.startPositions(), outputPath);
    }

    RatioSegmenter(ListMultimap<Chromosome, CobaltRatio> ratios, ChrArmLocator chrArmLocator, double gamma)
    {
        ratios.keySet().forEach(chromosome ->
        {
            List<CobaltRatio> ratiosForChromosome = ratios.get(chromosome);
            ratiosForChromosome.forEach(cobaltRatio ->
            {
                if(value(cobaltRatio) > 0.0)
                {
                    ArmToRatios.put(chrArmLocator.map(cobaltRatio), cobaltRatio);
                }
            });

        });
        mGamma = gamma;
    }

    Map<ChrArm, PiecewiseConstantFit> getSegmentation(ExecutorService executor) throws Exception
    {
        final HashMap<ChrArm, PiecewiseConstantFit> result = new HashMap<>();
        List<Future<Pair<ChrArm, PiecewiseConstantFit>>> futures = Lists.newArrayList();

        ArmToRatios.keySet().forEach(chr ->
        {
            Callable<Pair<ChrArm, PiecewiseConstantFit>> task = () -> calculateSegmentsForChromosome(ArmToRatios.get(chr), chr);
            futures.add(executor.submit(task));
        });

        for(Future<Pair<ChrArm, PiecewiseConstantFit>> future : futures)
        {
            Pair<ChrArm, PiecewiseConstantFit> chrPCF = future.get();
            if(chrPCF != null)
            {
                result.put(chrPCF.getLeft(), chrPCF.getRight());
            }
        }
        return result;
    }

    Map<ChrArm, Integer> startPositions()
    {
        return ArmToIndexOfStart;
    }

    abstract double value(CobaltRatio ratio);

    private Pair<ChrArm, PiecewiseConstantFit> calculateSegmentsForChromosome(List<CobaltRatio> ratios, ChrArm arm)
    {
        if(ratios.isEmpty())
        {
            return null;
        }
        ArmToIndexOfStart.put(arm, ratios.get(0).position());
        double[] d = new double[ratios.size()];
        for(int i = 0; i < ratios.size(); i++)
        {
            CobaltRatio ratio = ratios.get(i);
            d[i] = FastMath.log(2, value(ratio));
        }
        Segmenter segmenter = new Segmenter(d, mGamma, true);
        PiecewiseConstantFit fit = segmenter.pcf();
        return Pair.of(arm, fit);
    }
}

class TumorRatioSegmenter extends RatioSegmenter
{
    TumorRatioSegmenter(final ListMultimap<Chromosome, CobaltRatio> ratios, final ChrArmLocator chrArmLocator, final double gamma)
    {
        super(ratios, chrArmLocator, gamma);
    }

    @Override
    double value(final CobaltRatio ratio)
    {
        return ratio.tumorGCRatio();
    }
}

class ReferenceRatioSegmenter extends RatioSegmenter
{
    ReferenceRatioSegmenter(final ListMultimap<Chromosome, CobaltRatio> ratios, final ChrArmLocator chrArmLocator, final double gamma)
    {
        super(ratios, chrArmLocator, gamma);
    }

    @Override
    double value(final CobaltRatio ratio)
    {
        return ratio.referenceGCDiploidRatio();
    }
}
