package com.hartwig.hmftools.knowledgebasegenerator.vicc;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.junit.Test;

public class ViccHotspotExtractorTestApplicationTest {

    @Test
    public void canMergeIntoExistingMap() {
        Map<String, List<String>> map1 = Maps.newHashMap();
        map1.put("key", Lists.newArrayList("1", "2", "3", "4"));
        Map<String, List<String>> map2 = Maps.newHashMap();
        map2.put("key", Lists.newArrayList("4", "5", "6"));

        Map<String, List<String>> map3 = Maps.newHashMap();
        ViccHotspotExtractorTestApplication.mergeIntoExistingMap(map3, map1);
        ViccHotspotExtractorTestApplication.mergeIntoExistingMap(map3, map2);

        assertEquals(7, map3.get("key").size());
    }
}