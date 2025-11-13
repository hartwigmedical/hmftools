package com.hartwig.hmftools.cobalt.segmentation;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.IOException;
import java.io.Writer;
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
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.segmentation.PiecewiseConstantFit;
import com.hartwig.hmftools.common.segmentation.Segmenter;
import com.hartwig.hmftools.common.utils.pcf.CobaltSegment;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.FastMath;

public abstract class RatioSegmenter
{
    private final ListMultimap<ChrArm, CobaltRatio> ArmToRatios = ArrayListMultimap.create();
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
        Map<ChrArm, CobaltSegments> segmentsByChrArm = segmenter.getSegmentation(executor);
        write(segmentsByChrArm, outputPath);
    }

    RatioSegmenter(ListMultimap<Chromosome, CobaltRatio> ratios, ChrArmLocator chrArmLocator, double gamma)
    {
        ratios.keySet().forEach(chromosome ->
        {
            List<CobaltRatio> ratiosForChromosome = ratios.get(chromosome);
            ratiosForChromosome.forEach(cobaltRatio ->
            {
                if(value(cobaltRatio) >= 0.0)
                {
                    ArmToRatios.put(chrArmLocator.map(cobaltRatio), cobaltRatio);
                }
            });

        });
        mGamma = gamma;
    }

    private static void write(Map<ChrArm, CobaltSegments> data, String filename) throws IOException
    {
        List<ChrArm> armsInOrder = data.keySet().stream().sorted().toList();
        try(Writer writer = createBufferedWriter(filename))
        {
            String header = "Chromosome\tStart\tEnd\tMeanRatio\n";
            writer.write(header);
            for(ChrArm chrArm : armsInOrder)
            {
                CobaltSegments pcf = data.get(chrArm);
                String chromosome = chrArm.chromosome().shortName();
                for(CobaltSegment interval : pcf.Segments)
                {
                    int start = interval.start();
                    int end = interval.end();
                    writer.append(chromosome);
                    writer.append('\t');
                    writer.append(String.valueOf(start));
                    writer.append('\t');
                    writer.append(String.valueOf(end));
                    writer.append('\t');
                    writer.append(CobaltRatioFile.FORMAT.format(interval.MeanRatio));
                    writer.append('\n');
                }
            }
        }
    }

    private Map<ChrArm, CobaltSegments> getSegmentation(ExecutorService executor) throws Exception
    {
        final HashMap<ChrArm, CobaltSegments> result = new HashMap<>();
        List<Future<Pair<ChrArm, CobaltSegments>>> futures = Lists.newArrayList();

        ArmToRatios.keySet().forEach(chr ->
        {
            Callable<Pair<ChrArm, CobaltSegments>> task = () -> calculateSegmentsForChromosome(ArmToRatios.get(chr), chr);
            futures.add(executor.submit(task));
        });

        for(Future<Pair<ChrArm, CobaltSegments>> future : futures)
        {
            Pair<ChrArm, CobaltSegments> chrPCF = future.get();
            if(chrPCF != null)
            {
                result.put(chrPCF.getLeft(), chrPCF.getRight());
            }
        }
        return result;
    }

    abstract double value(CobaltRatio ratio);

    private Pair<ChrArm, CobaltSegments> calculateSegmentsForChromosome(List<CobaltRatio> ratios, ChrArm arm)
    {
        if(ratios.isEmpty())
        {
            return null;
        }
        double[] valuesForSegmentation = new double[ratios.size()];
        double[] rawValues = new double[ratios.size()];
        for(int i = 0; i < ratios.size(); i++)
        {
            CobaltRatio ratio = ratios.get(i);
            final double v = value(ratio);
            rawValues[i] = v;
            // Our R script that called copynumber put 0.001 as a floor for the ratios and converted
            // them to log_2 values. Note that negative ratio values have already been filtered out.
            if(v < 0.001)
            {
                valuesForSegmentation[i] = -9.965784;
            }
            else {
                valuesForSegmentation[i] = (float)FastMath.log(2, v);
            }
        }
        Segmenter segmenter = new Segmenter(valuesForSegmentation, mGamma, true);
        PiecewiseConstantFit fit = segmenter.pcf();
        CobaltSegments segments = new CobaltSegments(arm, ratios, rawValues, fit);
        return Pair.of(arm, segments);
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
