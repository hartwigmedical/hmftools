package com.hartwig.hmftools.amber.contamination;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.amber.AmberBase;
import com.hartwig.hmftools.common.amber.BaseDepthData;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

public class TumorContaminationModelTest
{
    private static final double EPSILON = 1e-10;

    private static final String CONTAMINATION_HIGH = Resources.getResource("amber/contamination.high").getPath();
    private static final String CONTAMINATION_LOW = Resources.getResource("amber/contamination.low").getPath();
    private static final String CONTAMINATION_NONE = Resources.getResource("amber/contamination.none").getPath();

    @Test
    public void testTumorContaminationVaf()
    {
        BaseDepthData tumorBdd = new BaseDepthData(
                AmberBase.T,
                AmberBase.G,
                100,
                0,
                95,
                5);
        TumorContamination tc = new TumorContamination("chr1", 1000, null, tumorBdd);
        Assert.assertEquals(0.05, tc.tumorVaf(), 0.001);
    }

    @Test
    public void testHighContamination() throws IOException
    {
        Map<Integer, Long> contaminationMap = fromFile(CONTAMINATION_HIGH);
        assertEquals(249329, calcAltReadCountAboveThreshold(2, contaminationMap));
        assertEquals(248936, calcAltReadCountAboveThreshold(3, contaminationMap));

        double contamination = calcContamination(107, contaminationMap);
        assertEquals(1, contamination, EPSILON);
    }

    @Test
    public void testLowContamination() throws IOException
    {
        Map<Integer, Long> contaminationMap = fromFile(CONTAMINATION_LOW);
        assertEquals(73058, calcAltReadCountAboveThreshold(2, contaminationMap));
        assertEquals(31170, calcAltReadCountAboveThreshold(3, contaminationMap));

        double contamination = calcContamination(99, contaminationMap);
        assertEquals(0.02, contamination, EPSILON);
    }

    /* would need to translate the alt frequencies into actual TumorContamination records to test this condition
    @Test
    public void testNoContamination() throws IOException
    {
        Map<Integer,Long> contaminationMap = fromFile(CONTAMINATION_NONE);
        assertEquals(2860, calcAltReadCountAboveThreshold(2, contaminationMap));
        assertEquals(374, calcAltReadCountAboveThreshold(3, contaminationMap));

        double contamination = calcContamination(107, contaminationMap);
        assertEquals(0, contamination, EPSILON);
    }
    */

    public static double calcContamination(int medianTumorReadDepth, final Map<Integer, Long> altSupportFrequencies)
    {
        long twoPlusReadCount = calcAltReadCountAboveThreshold(2, altSupportFrequencies);
        TumorContaminationModel model = new TumorContaminationModel();
        return model.calcContamination(medianTumorReadDepth, twoPlusReadCount, altSupportFrequencies);
    }

    static long calcAltReadCountAboveThreshold(int minAltSupport, final Map<Integer, Long> altSupportMap)
    {
        return altSupportMap.entrySet().stream().filter(x -> x.getKey() >= minAltSupport).mapToLong(Map.Entry::getValue).sum();
    }


    @NotNull
    private static Map<Integer, Long> fromFile(final String fileName) throws IOException
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
