package com.hartwig.hmftools.common.cuppa;

import static org.junit.Assert.*;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOriginFactory;

import org.junit.Test;

public class MolecularTissueOriginFactoryTest {

    private static final String MOLECULAR_TISSUE_ORIGIN_TXT = Resources.getResource("cuppa/sample.cuppa.conclusion.txt").getPath();

    @Test
    public void canReadMolecularTissueOriginTsv() throws IOException {
        String molecularTissueOrigin = MolecularTissueOriginFactory.readMolecularTissueOriginResult(MOLECULAR_TISSUE_ORIGIN_TXT);

        assertEquals("results inconclusive", molecularTissueOrigin);

    }

}