package com.hartwig.hmftools.ecrfanalyser;

import java.io.File;
import java.io.IOException;
import java.util.List;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Ignore;
import org.junit.Test;

public class EcrfAnalysisApplicationTest {

    private static final String TEST_ECRF = Resources.getResource("example").getPath() + File.separator + "ecrf.xml";
    private static final List<String> PATIENTS = Lists.newArrayList("CPCT02252500");

    @Test
    @Ignore
    public void runApplication() throws IOException, XMLStreamException {
        new EcrfAnalysisApplication(TEST_ECRF, null, PATIENTS, null, false, true).run();
    }

    @Test
    @Ignore
    public void runWithFieldSelectionAndRowsTransposed() throws IOException, XMLStreamException {
        final List<String> fields = Lists.newArrayList("AFTERBIOPT.TRTAFTER.TRTAFTER.SYSREGPOST",
                "BASELINE.CARCINOMA.CARCINOMA.PTUMLOC");
        new EcrfAnalysisApplication(TEST_ECRF, null, PATIENTS, fields, true, true).run();
    }
}