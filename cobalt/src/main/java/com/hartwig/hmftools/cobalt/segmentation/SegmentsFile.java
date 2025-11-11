package com.hartwig.hmftools.cobalt.segmentation;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.IOException;
import java.io.Writer;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.cobalt.CobaltConstants;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.segmentation.PiecewiseConstantFit;

class SegmentsFile
{
    static void write(Map<ChrArm, PiecewiseConstantFit> data,
            Map<ChrArm, Integer> armStartPositions,
            String filename) throws IOException
    {
        int windowSize = CobaltConstants.WINDOW_SIZE;
        List<ChrArm> armsInOrder = data.keySet().stream().sorted().toList();
        try(Writer writer = createBufferedWriter(filename))
        {
            String header = "Chromosome\tStart\tEnd\n";
            writer.write(header);
            for(ChrArm chrArm : armsInOrder)
            {
                PiecewiseConstantFit pcf = data.get(chrArm);
                Integer positionOfFirstWindow = armStartPositions.get(chrArm);
                String chromosome = chrArm.chromosome().shortName();
                PCFCoordinates coordinates = new PCFCoordinates(pcf, windowSize, positionOfFirstWindow, chromosome);
                for(ChrBaseRegion interval : coordinates.intervals())
                {
                    String chr = interval.chromosome();
                    int start = interval.start();
                    int end = interval.end();
                    writer.append(chr);
                    writer.append('\t');
                    writer.append(String.valueOf(start));
                    writer.append('\t');
                    writer.append(String.valueOf(end));
                    writer.append('\n');
                }
            }
        }
    }
}
