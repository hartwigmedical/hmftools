package com.hartwig.hmftools.common.lims.cohort;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.lims.LimsTestUtil;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class LimsCohortModelTest {

    @Test
    public void canQueryLimsCohortModel() {
        LimsCohortConfig cohortConfig = LimsTestUtil.createAllDisabledCohortConfig("DRUP");
        Map<String, LimsCohortConfig> cohortMap = Maps.newHashMap();
        cohortMap.put("DRUP", cohortConfig);

        LimsCohortModel model = ImmutableLimsCohortModel.builder().limsCohortMap(cohortMap).build();
        assertEquals(cohortConfig, model.queryCohortData("DRUP", "DRUP01"));
    }

    @Test
    public void fallbackToSampleIdWhenCohortIsEmpty() {
        LimsCohortConfig cohortConfig = LimsTestUtil.createAllDisabledCohortConfig("DRUP");
        Map<String, LimsCohortConfig> cohortMap = Maps.newHashMap();
        cohortMap.put("DRUP", cohortConfig);

        LimsCohortModel model = ImmutableLimsCohortModel.builder().limsCohortMap(cohortMap).build();
        assertEquals(cohortConfig, model.queryCohortData(null, "DRUP01"));
        assertEquals(cohortConfig, model.queryCohortData(Strings.EMPTY, "DRUP01"));
    }

    @Test(expected = IllegalStateException.class)
    public void crashWhenSampleIdIsNotAsExpected() {
        Map<String, LimsCohortConfig> cohortMap = Maps.newHashMap();
        cohortMap.put("DRUP", LimsTestUtil.createAllDisabledCohortConfig("DRUP"));

        LimsCohortModel model = ImmutableLimsCohortModel.builder().limsCohortMap(cohortMap).build();
        model.queryCohortData("DRUP", "CPCT01");
    }

    @Test(expected = IllegalStateException.class)
    public void crashWhenCohortConfigIsNotPresent() {
        Map<String, LimsCohortConfig> cohortMap = Maps.newHashMap();
        cohortMap.put("DRUP", LimsTestUtil.createAllDisabledCohortConfig("DRUP"));

        LimsCohortModel model = ImmutableLimsCohortModel.builder().limsCohortMap(cohortMap).build();
        model.queryCohortData("NON-DRUP", "NON-DRUP01");
    }
}