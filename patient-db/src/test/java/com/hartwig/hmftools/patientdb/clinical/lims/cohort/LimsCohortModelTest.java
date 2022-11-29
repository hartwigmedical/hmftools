package com.hartwig.hmftools.patientdb.clinical.lims.cohort;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class LimsCohortModelTest {

    @Test
    public void canQueryLimsCohortModel() {
        LimsCohortConfig cohortConfig = LimsCohortTestFactory.createAllDisabledCohortConfig("DRUP");
        Map<String, LimsCohortConfig> cohortMap = Maps.newHashMap();
        cohortMap.put("DRUP", cohortConfig);

        LimsCohortModel model = ImmutableLimsCohortModel.builder().limsCohortMap(cohortMap).build();
        assertEquals(cohortConfig, model.queryCohortData("DRUP", "DRUP01"));
    }

    @Test
    public void fallbackToSampleIdWhenCohortIsEmpty() {
        LimsCohortConfig cohortConfig = LimsCohortTestFactory.createAllDisabledCohortConfig("DRUP");
        Map<String, LimsCohortConfig> cohortMap = Maps.newHashMap();
        cohortMap.put("DRUP", cohortConfig);

        LimsCohortModel model = ImmutableLimsCohortModel.builder().limsCohortMap(cohortMap).build();
        assertEquals(cohortConfig, model.queryCohortData(null, "DRUP01"));
        assertEquals(cohortConfig, model.queryCohortData(Strings.EMPTY, "DRUP01"));
    }

    @Test(expected = IllegalStateException.class)
    public void crashWhenSampleIdIsNotAsExpected() {
        Map<String, LimsCohortConfig> cohortMap = Maps.newHashMap();
        cohortMap.put("DRUP", LimsCohortTestFactory.createAllDisabledCohortConfig("DRUP"));

        LimsCohortModel model = ImmutableLimsCohortModel.builder().limsCohortMap(cohortMap).build();
        model.queryCohortData("DRUP", "CPCT01");
    }

    @Test
    public void nullWhenCohortConfigIsNotPresent() {
        Map<String, LimsCohortConfig> cohortMap = Maps.newHashMap();
        cohortMap.put("DRUP", LimsCohortTestFactory.createAllDisabledCohortConfig("DRUP"));

        LimsCohortModel model = ImmutableLimsCohortModel.builder().limsCohortMap(cohortMap).build();
        assertNull(model.queryCohortData("NON-DRUP", "NON-DRUP01"));
    }
}