package com.hartwig.hmftools.common.cuppa;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class MolecularTissueOriginFileTest {

    private static final String MOLECULAR_TISSUE_ORIGIN_TXT = Resources.getResource("cuppa/sample.cuppa.conclusion.txt").getPath();

    @Test
    public void canReadMolecularTissueOriginTxt() throws IOException {
        MolecularTissueOrginData molecularTissueOrigin = MolecularTissueOriginFile.read(MOLECULAR_TISSUE_ORIGIN_TXT);

        assertEquals("results inconclusive", molecularTissueOrigin.conclusion());
    }

    @Test
    public void canDetermineMolecularTissueOriginDta() {
        MolecularTissueOrginData molecularTissueOrigin = MolecularTissueOriginFile.extractPedictionDataOrigin("Lower GI tract (likelihood=80.4%)");

        assertEquals("Lower GI tract (likelihood=80.4%)", molecularTissueOrigin.conclusion());
        assertEquals("Lower GI tract", molecularTissueOrigin.predictedOrigin());
        assertEquals("likelihood=80.4%", molecularTissueOrigin.predictionLikelihood());

    }
}