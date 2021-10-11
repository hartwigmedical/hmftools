package com.hartwig.hmftools.patientdb.clinical.readers.cpct;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.patientdb.clinical.consents.ConsentConfigFactory;
import com.hartwig.hmftools.patientdb.clinical.curators.TestCuratorFactory;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BaselineData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;

import org.junit.Test;

public class BaselineReaderTest {

    private static final String INFORMED_CONSENTS_TSV = Resources.getResource("consents/informed_consents.tsv").getPath();

    @Test
    public void canReadEmptyBaselineInfo() throws IOException {
        Map<Integer, String> hospitals = Maps.newHashMap();
        BaselineReader baselineReader =
                new BaselineReader(TestCuratorFactory.primaryTumorCurator(), hospitals);

        EcrfPatient cpctPatient = new EcrfPatient("patient", Maps.newHashMap(), Lists.newArrayList());
        BaselineData baselineData = baselineReader.read(cpctPatient, ConsentConfigFactory.read(INFORMED_CONSENTS_TSV), "CPCT");
        assertNotNull(baselineData);
    }
}
