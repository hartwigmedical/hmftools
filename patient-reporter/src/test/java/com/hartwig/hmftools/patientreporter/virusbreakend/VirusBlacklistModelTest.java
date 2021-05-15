package com.hartwig.hmftools.patientreporter.virusbreakend;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;

import org.junit.Test;

public class VirusBlacklistModelTest {

    @Test
    public void canMatchVirusToId() {
        Map<Integer, String> virusBlacklistMap = Maps.newHashMap();
        virusBlacklistMap.put(1, "virus1");
        VirusBlacklistModel VirusBlackListModel = new VirusBlacklistModel(virusBlacklistMap);

        assertTrue(VirusBlackListModel.checkVirusForBlacklisting(1));
        assertFalse(VirusBlackListModel.checkVirusForBlacklisting(4));
    }
}