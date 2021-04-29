package com.hartwig.hmftools.patientreporter.cuppa;

import static org.junit.Assert.*;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class MolecularTissueOriginFactoryTest {

    private static final String MOLECULAR_TISSUE_ORIGIN_TXT = Resources.getResource("test_run/cuppa/sample.cuppa.conclusion.txt").getPath();

    @Test
    public void canReadMolecularTissueOriginTsv() throws IOException {
        String molecularTissueOrigin = MolecularTissueOriginFactory.readMolecularTissueOriginResult(MOLECULAR_TISSUE_ORIGIN_TXT);

        assertEquals("results inconclusive", molecularTissueOrigin);

    }

}