package com.hartwig.hmftools.amber;

import java.util.List;
import java.util.concurrent.ExecutorService;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;
import com.hartwig.hmftools.common.segmentation.DataForSegmentation;
import com.hartwig.hmftools.common.segmentation.PerArmSegmenter;
import com.hartwig.hmftools.common.segmentation.SegmentsFile;

public class BAFSegmenter extends PerArmSegmenter<AmberBAF>
{
    public static void writeSegments(final List<AmberBAF> amberBAFS,
            RefGenomeVersion genomeVersion,
            ExecutorService executor,
            String outputPath) throws Exception
    {
        ArrayListMultimap<Chromosome, AmberBAF> bafMap = ArrayListMultimap.create();
        amberBAFS.forEach(amberBAF -> bafMap.put(amberBAF.chr(), amberBAF));
        BAFSegmenter segmenter = new BAFSegmenter(bafMap, ChrArmLocator.defaultLocator(genomeVersion));
        var segmentation = segmenter.getSegmentation(executor);
        SegmentsFile.write(genomeVersion, segmentation, outputPath);
    }

    public BAFSegmenter(ListMultimap<Chromosome, AmberBAF> bafs, final ChrArmLocator chrArmLocator)
    {
        super(bafs, chrArmLocator, 100.0);
    }

    @Override
    public double value(final AmberBAF ratio)
    {
        return ratio.tumorModifiedBAF();
    }

    @Override
    public DataForSegmentation buildSegmentationData(final List<AmberBAF> ratios)
    {
        double[] valuesForSegmentation = new double[ratios.size()];
        double[] rawValues = new double[ratios.size()];
        for(int i = 0; i < ratios.size(); i++)
        {
            AmberBAF baf = ratios.get(i);
            final double v = value(baf);
            rawValues[i] = v;
            valuesForSegmentation[i] = v;
        }
        return new DataForSegmentation(valuesForSegmentation, rawValues);
    }

    @Override
    public boolean isWindowed()
    {
        return false;
    }
}
