package com.hartwig.hmftools.common.lims.cohort;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.lims.LimsTestUtil;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class LimsCohortModelTest {

    @Test
    public void canQueryLimsCohortModel() {
        LimsCohortConfig cohortConfigData = LimsTestUtil.createAllDisabledCohortConfig("DRUP");
        Map<String, LimsCohortConfig> cohortMap = Maps.newHashMap();
        cohortMap.put("DRUP", cohortConfigData);

        LimsCohortModel model = ImmutableLimsCohortModel.builder().limsCohortMap(cohortMap).build();
        assertEquals(cohortConfigData, model.queryCohortData("DRUP", "DRUP01", Strings.EMPTY));
    }

    @Test(expected = IllegalStateException.class)
    public void canQueryLimsCohortModelException() {
        LimsCohortConfig cohortConfigData = LimsTestUtil.createAllDisabledCohortConfig("DRUP");
        Map<String, LimsCohortConfig> cohortMap = Maps.newHashMap();
        cohortMap.put("DRUP", cohortConfigData);

        LimsCohortModel model = ImmutableLimsCohortModel.builder().limsCohortMap(cohortMap).build();
        assertNull(model.queryCohortData(null, "DRUP01", Strings.EMPTY));
        assertNull(model.queryCohortData("CPCT", "DRUP01", Strings.EMPTY));
        assertNull(model.queryCohortData("DRUP", "CPCT01", Strings.EMPTY));
    }
}