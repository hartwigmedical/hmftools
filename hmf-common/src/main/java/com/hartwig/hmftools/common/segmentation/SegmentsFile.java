package com.hartwig.hmftools.common.segmentation;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.IOException;
import java.io.Writer;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.pcf.PcfSegment;

public class SegmentsFile
{
    public static <T extends GenomePosition> void write(
            RefGenomeVersion genomeVersion,
            Map<ChrArm, ChromosomeArmSegments<T>> data,
            String filename) throws IOException
    {
        List<ChrArm> armsInOrder = data.keySet().stream().sorted().toList();
        try(Writer writer = createBufferedWriter(filename))
        {
            String header = "Chromosome\tStart\tEnd\tMeanRatio\n";
            writer.write(header);
            for(ChrArm chrArm : armsInOrder)
            {
                ChromosomeArmSegments<?> pcf = data.get(chrArm);
                String chromosome = genomeVersion.versionedChromosome(chrArm.chromosome());
                for(PcfSegment interval : pcf.Segments)
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
}
