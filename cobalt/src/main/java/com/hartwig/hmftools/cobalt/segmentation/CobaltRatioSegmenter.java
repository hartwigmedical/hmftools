package com.hartwig.hmftools.cobalt.segmentation;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;
import com.hartwig.hmftools.common.segmentation.ChromosomeArmSegments;
import com.hartwig.hmftools.common.segmentation.DataForSegmentation;
import com.hartwig.hmftools.common.segmentation.PerArmSegmenter;
import com.hartwig.hmftools.common.segmentation.SegmentsFile;

import org.apache.commons.math3.util.FastMath;

public abstract class CobaltRatioSegmenter extends PerArmSegmenter<CobaltRatio>
{
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
        PerArmSegmenter<CobaltRatio> segmenter =
                isForTumor ? new TumorRatioSegmenter(ratios, locator, gamma) : new ReferenceRatioSegmenter(ratios, locator, gamma);
        Map<ChrArm, ChromosomeArmSegments<CobaltRatio>> segmentsByChrArm = segmenter.getSegmentation(executor);
        SegmentsFile.write(genomeVersion, segmentsByChrArm, outputPath);
    }

    CobaltRatioSegmenter(final ListMultimap<Chromosome, CobaltRatio> ratios, final ChrArmLocator chrArmLocator, final double gamma)
    {
        super(ratios, chrArmLocator, gamma);
    }

    @Override
    public DataForSegmentation buildSegmentationData(final List<CobaltRatio> ratios)
    {
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
            else
            {
                valuesForSegmentation[i] = (float) FastMath.log(2, v);
            }
        }

        return new DataForSegmentation(valuesForSegmentation, rawValues);
    }

    @Override
    public final boolean isWindowed()
    {
        return true;
    }
}

