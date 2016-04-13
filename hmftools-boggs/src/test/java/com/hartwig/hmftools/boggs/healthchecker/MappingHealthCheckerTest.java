package com.hartwig.hmftools.boggs.healthchecker;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.boggs.PatientData;
import com.hartwig.hmftools.boggs.SampleData;
import com.hartwig.hmftools.boggs.flagstatreader.FlagStatData;
import com.hartwig.hmftools.boggs.flagstatreader.FlagStatTestFactory;
import com.hartwig.hmftools.boggs.healthcheck.HealthChecker;
import com.hartwig.hmftools.boggs.healthcheck.MappingHealthChecker;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class MappingHealthCheckerTest {

    @Test
    public void verifyMappingHealthChecker() {
        HealthChecker checker = new MappingHealthChecker();

        PatientData patient = new PatientData(dummyData(), dummyData());

        assertTrue(checker.isHealthy(patient));
    }

    @NotNull
    private static SampleData dummyData() {
        FlagStatData testData = FlagStatTestFactory.createTestData();
        return new SampleData("DUMMY", Lists.newArrayList(testData), testData, testData);

    }
}
