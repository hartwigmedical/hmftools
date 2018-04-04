package com.hartwig.hmftools.patientdb.readers.cpct;

import static org.junit.Assert.assertNotNull;

import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.curators.TestCuratorFactory;
import com.hartwig.hmftools.patientdb.data.BaselineData;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BaselineReaderTest {

    @Test
    public void canReadEmptyBaselineInfo() {
        final Map<Integer, String> hospitals = Maps.newHashMap();
        final BaselineReader baselineReader = new BaselineReader(TestCuratorFactory.tumorLocationCurator(), hospitals);

        final EcrfPatient cpctPatient = emptyPatient();
        final BaselineData baselineData = baselineReader.read(cpctPatient);
        assertNotNull(baselineData);
    }

    @NotNull
    private static EcrfPatient emptyPatient() {
        return new EcrfPatient("patient", Maps.newHashMap(), Lists.newArrayList());
    }
}
