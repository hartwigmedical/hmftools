package com.hartwig.hmftools.amber;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TumorContaminationModelTest
{

    private static final double EPSILON = 1e-10;

    private static final String CONTAMINATION_HIGH = Resources.getResource("amber/contamination.high").getPath();
    private static final String CONTAMINATION_LOW = Resources.getResource("amber/contamination.low").getPath();
    private static final String CONTAMINATION_NONE = Resources.getResource("amber/contamination.none").getPath();

    @Test
    public void testHighContamination() throws IOException
    {
        final Map<Integer, Long> contaminationMap = fromFile(CONTAMINATION_HIGH);
        final TumorContaminationModel model = new TumorContaminationModel();
        assertEquals(249329, TumorContaminationModel.reads(2, contaminationMap));
        assertEquals(248936, TumorContaminationModel.reads(3, contaminationMap));

        final double contamination = model.contamination(107, contaminationMap);
        assertEquals(1, contamination, EPSILON);
    }

    @Test
    public void testLowContamination() throws IOException
    {
        final Map<Integer, Long> contaminationMap = fromFile(CONTAMINATION_LOW);
        final TumorContaminationModel model = new TumorContaminationModel();
        assertEquals(73058, TumorContaminationModel.reads(2, contaminationMap));
        assertEquals(31170, TumorContaminationModel.reads(3, contaminationMap));

        final double contamination = model.contamination(99, contaminationMap);
        assertEquals(0.02, contamination, EPSILON);
    }

    @Test
    public void testNoContamination() throws IOException
    {
        final Map<Integer, Long> contaminationMap = fromFile(CONTAMINATION_NONE);
        final TumorContaminationModel model = new TumorContaminationModel();
        assertEquals(2860, TumorContaminationModel.reads(2, contaminationMap));
        assertEquals(374, TumorContaminationModel.reads(3, contaminationMap));

        final double contamination = model.contamination(107, contaminationMap);
        assertEquals(0, contamination, EPSILON);
    }

    @Test
    public void testMedianReturnsZeroIfListEmpty()
    {
        assertEquals(0, TumorContaminationModel.medianDepth(Collections.emptyList()));
    }

    @NotNull
    private static Map<Integer, Long> fromFile(@NotNull final String fileName) throws IOException
    {
        Map<Integer, Long> result = Maps.newHashMap();

        List<String> lines = Files.readAllLines(new File(fileName).toPath());

        for(String line : lines)
        {
            String[] split = line.trim().split("\\s+");
            result.put(Integer.valueOf(split[1]), Long.valueOf(split[0]));
        }

        return result;
    }
}
