package com.hartwig.hmftools.ecrfanalyser;

import java.io.IOException;
import java.util.List;

import javax.xml.stream.XMLStreamException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.junit.Test;

public class EcrfAnalysisApplicationTest {

    private static final String TEST_ECRF = Resources.getResource("example/ecrf.xml").getPath();
    private static final String CSV_OUT = "/Users/kduyvesteyn/hmf/tmp/ecrf.csv";

    private static final List<String> PATIENTS = Lists.newArrayList("CPCT02252500");

    @Test
    public void tryIt() throws IOException, XMLStreamException {
        EcrfAnalysisApplication app = new EcrfAnalysisApplication(TEST_ECRF, CSV_OUT, PATIENTS);

        app.run();
    }
}