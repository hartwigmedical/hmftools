package com.hartwig.hmftools.esvee;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
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

        existingJunctions.add(new Junction(CHR_1, 100, POS_ORIENT));
        existingJunctions.add(new Junction(CHR_1, 110, POS_ORIENT));
        existingJunctions.add(new Junction(CHR_1, 120, POS_ORIENT));
        existingJunctions.add(new Junction(CHR_1, 120, NEG_ORIENT));
        existingJunctions.add(new Junction(CHR_1, 130, NEG_ORIENT));
        existingJunctions.add(new Junction(CHR_1, 140, NEG_ORIENT));

        newJunctions.add(new Junction(CHR_1, 95, POS_ORIENT));
        newJunctions.add(new Junction(CHR_1, 115, POS_ORIENT));
        newJunctions.add(new Junction(CHR_1, 120, POS_ORIENT)); // dup
        newJunctions.add(new Junction(CHR_1, 120, NEG_ORIENT)); // dup
        newJunctions.add(new Junction(CHR_1, 130, POS_ORIENT));
        newJunctions.add(new Junction(CHR_1, 135, NEG_ORIENT));
        newJunctions.add(new Junction(CHR_1, 140, NEG_ORIENT)); // dup
        newJunctions.add(new Junction(CHR_1, 145, NEG_ORIENT));

        Junction.mergeJunctions(existingMap, newMap);

        assertEquals(11, existingJunctions.size());
        assertTrue(validateJunctionMap(existingMap));
    }
}
