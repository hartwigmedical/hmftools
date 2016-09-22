package com.hartwig.hmftools.ecrfanalyser;

import static org.junit.Assert.assertTrue;

import java.io.FileNotFoundException;

import javax.xml.stream.XMLStreamException;

import com.google.common.io.Resources;

import org.junit.Test;

public class EcrfAnalysisApplicationTest {

    private static final String TEST_ECRF = Resources.getResource("example/ecrf.xml").getPath();

    @Test
    public void tryIt() throws FileNotFoundException, XMLStreamException {
        EcrfAnalysisApplication app = new EcrfAnalysisApplication();
        app.run(TEST_ECRF);
        assertTrue(true);
    }
}