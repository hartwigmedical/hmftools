package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.assembly.types.Junction.validateJunctionMap;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.assembly.types.Junction;

import org.junit.Test;

public class JunctionGroupsTest
{
    @Test
    public void testJunctionGroupMerging()
    {
        List<Junction> existingJunctions = Lists.newArrayList();
        List<Junction> newJunctions = Lists.newArrayList();

        Map<String,List<Junction>> existingMap = Maps.newHashMap();
        existingMap.put(CHR_1, existingJunctions);

        Map<String,List<Junction>> newMap = Maps.newHashMap();
        newMap.put(CHR_1, newJunctions);

        existingJunctions.add(new Junction(CHR_1, 100, FORWARD));
        existingJunctions.add(new Junction(CHR_1, 110, FORWARD));
        existingJunctions.add(new Junction(CHR_1, 120, FORWARD));
        existingJunctions.add(new Junction(CHR_1, 120, REVERSE));
        existingJunctions.add(new Junction(CHR_1, 130, REVERSE));
        existingJunctions.add(new Junction(CHR_1, 140, REVERSE));

        newJunctions.add(new Junction(CHR_1, 95, FORWARD));
        newJunctions.add(new Junction(CHR_1, 115, FORWARD));
        newJunctions.add(new Junction(CHR_1, 120, FORWARD)); // dup
        newJunctions.add(new Junction(CHR_1, 120, REVERSE)); // dup
        newJunctions.add(new Junction(CHR_1, 130, FORWARD));
        newJunctions.add(new Junction(CHR_1, 135, REVERSE));
        newJunctions.add(new Junction(CHR_1, 140, REVERSE)); // dup
        newJunctions.add(new Junction(CHR_1, 145, REVERSE));

        Junction.mergeJunctions(existingMap, newMap);

        assertEquals(11, existingJunctions.size());
        assertTrue(validateJunctionMap(existingMap));
    }
}
