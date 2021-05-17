package com.hartwig.hmftools.common.cuppa;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class MolecularTissueOriginFileTest {

    private static final String MOLECULAR_TISSUE_ORIGIN_TXT = Resources.getResource("cuppa/sample.cuppa.conclusion.txt").getPath();

    @Test
    public void canReadMolecularTissueOriginTxt() throws IOException {
        String molecularTissueOrigin = MolecularTissueOriginFile.read(MOLECULAR_TISSUE_ORIGIN_TXT);

        assertEquals("results inconclusive", molecularTissueOrigin);
    }
}