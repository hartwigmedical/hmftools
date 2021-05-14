package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;

import org.junit.Test;

public class VirusBlackListModel2Test {

    @Test
    public void canMatchVirusToId() {
        Map<Integer, String> virusBlacklistMap = Maps.newHashMap();
        virusBlacklistMap.put(1, "virus1");
        VirusBlackListModel2 VirusBlackListModel = new VirusBlackListModel2(virusBlacklistMap);

        assertTrue(VirusBlackListModel.checkVirusForBlacklisting(1));
        assertFalse(VirusBlackListModel.checkVirusForBlacklisting(4));
    }
}