package com.hartwig.hmftools.common.ecrf;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.formstatus.ImmutableFormStatusModel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CpctEcrfModelTest {

    private static final Logger LOGGER = LogManager.getLogger(CpctEcrfModel.class);

    private static final String BASE_RESOURCE_DIR = Resources.getResource("ecrf").getPath();
    private static final String TEST_ECRF = BASE_RESOURCE_DIR + File.separator + "example" + File.separator + "ecrf.xml";

    @Test
    public void loadDataFromRealXML() throws IOException, XMLStreamException {
        if (Files.exists(new File(TEST_ECRF).toPath())) {
            final CpctEcrfModel model = CpctEcrfModel.loadFromXML(TEST_ECRF, new ImmutableFormStatusModel(Maps.newHashMap()));
            assertTrue(hasField(model, "BASELINE.CARCINOMA.CARCINOMA.PTUMLOC"));
            assertFalse(hasField(model, "Does Not Exist"));
            assertNotNull(model.datamodel().studyEvents().get("SE.BASELINE"));
            assertNotNull(model.datamodel().forms().get("FRM.CARCINOMA"));
            assertNotNull(model.datamodel().itemGroups().get("GRP.CARCINOMA.CARCINOMA"));
            assertNotNull(model.datamodel().items().get("FLD.CARCINOMA.PTUMLOC"));
            assertNull(model.datamodel().studyEvents().get("Does Not Exist"));
            assertNull(model.datamodel().forms().get("Does Not Exist"));
            assertNull(model.datamodel().itemGroups().get("Does Not Exist"));
            assertNull(model.datamodel().items().get("Does Not Exist"));

            assertTrue(hasPatient(model, "CPCT02252500"));
            assertFalse(hasPatient(model, "Does Not Exist"));
        } else {
            LOGGER.warn("Could not run CpctEcrfModelTest since test xml is not present: " + TEST_ECRF);
        }
    }

    private static boolean hasField(@NotNull final CpctEcrfModel model, @NotNull final String fieldId) {
        return Lists.newArrayList(model.findFieldsById(Lists.newArrayList(fieldId))).size() > 0;
    }

    private static boolean hasPatient(@NotNull final CpctEcrfModel model, @NotNull final String patientId) {
        final List<EcrfPatient> patients = Lists.newArrayList(model.findPatientsById(Lists.newArrayList(patientId)));
        return !patients.get(0).fields().isEmpty();
    }
}
