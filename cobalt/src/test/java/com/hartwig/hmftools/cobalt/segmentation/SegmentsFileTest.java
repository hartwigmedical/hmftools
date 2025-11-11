package com.hartwig.hmftools.cobalt.segmentation;

import static com.hartwig.hmftools.cobalt.segmentation.Arm.P;
import static com.hartwig.hmftools.cobalt.segmentation.Arm.Q;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.segmentation.PiecewiseConstantFit;

import org.apache.commons.io.FileUtils;
import org.junit.Assert;
import org.junit.Test;

public class SegmentsFileTest
{
    @Test
    public void writeSegmentationFileTest() throws Exception
    {
        File tempDir = Files.createTempDirectory("sft").toFile();
        File outputFile = new File(tempDir, "blah.pcf");
        Assert.assertFalse(outputFile.exists());
        double[] means = { 20.0, 30.0 };
        double[] moreMeans = { 20.0, 30.0, 40 };

        int[] lengths_1_1 = { 100, 200, 200 };
        int[] positions_1_1 = { 1, 2000, 3000 };
        PiecewiseConstantFit pcf_1_1 = new PiecewiseConstantFit(lengths_1_1, positions_1_1, moreMeans);

        int[] lengths_2_1 = { 100, 200 };
        int[] positions_2_1 = { 1, 2000 };
        PiecewiseConstantFit pcf_2_1 = new PiecewiseConstantFit(lengths_2_1, positions_2_1, means);

        int[] lengths_3_1 = { 10, 10 };
        int[] positions_3_1 = { 1, 21 };
        PiecewiseConstantFit pcf_3_1 = new PiecewiseConstantFit(lengths_3_1, positions_3_1, means);

        Map<ChrArm, PiecewiseConstantFit> data = new HashMap<>();
        data.put(new ChrArm(_2, P), pcf_2_1); // 3 rows
        data.put(new ChrArm(_1, P), pcf_1_1); // 2
        data.put(new ChrArm(_3, Q), pcf_3_1); // 2

        Map<ChrArm, Integer> startPositions = new HashMap<>();
        startPositions.put(new ChrArm(_2, P), 0);
        startPositions.put(new ChrArm(_1, P), 0);
        startPositions.put(new ChrArm(_3, Q), 0);

        SegmentsFile.write(data, startPositions, outputFile.getAbsolutePath());

        List<String> written = FileUtils.readLines(outputFile, StandardCharsets.UTF_8);
        assertEquals(7 + 1, written.size());
        assertEquals("Chromosome\tStart\tEnd", written.get(0));
        assertEquals("1\t1000\t100999", written.get(1));
        assertEquals("1\t2000000\t2199999", written.get(2));
        assertEquals("3\t21000\t30999", written.get(7));
    }
}
