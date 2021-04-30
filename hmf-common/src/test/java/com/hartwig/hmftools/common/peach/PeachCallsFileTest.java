package com.hartwig.hmftools.common.peach;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class PeachCallsFileTest {

    private static final String FILE = Resources.getResource("peach/sample_calls.txt").getPath();

    private static final String GENE_1 = "GENE";
    private static final String CHROMOSOME_1 = "1";
    private static final String POSITIONV37_1 = "10";
    private static final String POSITIONV38_1 = "12";
    private static final String REFV37_1 = "A";
    private static final String REFV38_1 = "A";
    private static final String ALLELE1_1 = "A";
    private static final String ALLELE2_1 = "A";
    private static final String RSID_1 = "rs12";
    private static final String VARIANTANNOTATIONV37_1 = "REF_CALL";
    private static final String FILTERV37_1 = "NO_CALL";
    private static final String VARIANTANNOTATIONV38_1 = "REF_CALL";
    private static final String FILTERV38_1 = "NO_CALL";
    private static final String PANEL_VERSION_1 = "PGx_min_DPYD_v0.3";
    private static final String REPO_VERSION_1 = "1.0";

    private static final String GENE_2 = "GENE";
    private static final String CHROMOSOME_2 = "1";
    private static final String POSITIONV37_2 = "100";
    private static final String POSITIONV38_2 = "110";
    private static final String REFV37_2 = "T";
    private static final String REFV38_2 = "T";
    private static final String ALLELE1_2 = "T";
    private static final String ALLELE2_2 = "T";
    private static final String RSID_2 = "rs54";
    private static final String VARIANTANNOTATIONV37_2 = "REF_CALL";
    private static final String FILTERV37_2 = "NO_CALL";
    private static final String VARIANTANNOTATIONV38_2 = "REF_CALL";
    private static final String FILTERV38_2 = "NO_CALL";
    private static final String PANEL_VERSION_2 = "PGx_min_DPYD_v0.3";
    private static final String REPO_VERSION_2 = "1.0";

    @Test
    public void loadPeachCallsFile() throws IOException {
        List<PeachCalls> peachCalls = PeachCallsFile.read(FILE);

        assertEquals(6, peachCalls.size());

        PeachCalls peachCalls1 = peachCalls.get(0);
        assertEquals(GENE_1, peachCalls1.gene());
        assertEquals(CHROMOSOME_1, peachCalls1.chromosome());
        assertEquals(POSITIONV37_1, peachCalls1.positionV37());
        assertEquals(POSITIONV38_1, peachCalls1.positionV38());
        assertEquals(REFV37_1, peachCalls1.refV37());
        assertEquals(REFV38_1, peachCalls1.refV38());
        assertEquals(ALLELE1_1, peachCalls1.allele1());
        assertEquals(ALLELE2_1, peachCalls1.allele2());
        assertEquals(RSID_1, peachCalls1.rsid());
        assertEquals(VARIANTANNOTATIONV37_1, peachCalls1.variantAnnotationV37());
        assertEquals(FILTERV37_1, peachCalls1.filterV37());
        assertEquals(VARIANTANNOTATIONV38_1, peachCalls1.variantAnnotationV38());
        assertEquals(FILTERV38_1, peachCalls1.filterV38());
        assertEquals(PANEL_VERSION_1, peachCalls1.panelVersion());
        assertEquals(REPO_VERSION_1, peachCalls1.repoVersion());

        PeachCalls peachCalls2 = peachCalls.get(1);
        assertEquals(GENE_2, peachCalls2.gene());
        assertEquals(CHROMOSOME_2, peachCalls2.chromosome());
        assertEquals(POSITIONV37_2, peachCalls2.positionV37());
        assertEquals(POSITIONV38_2, peachCalls2.positionV38());
        assertEquals(REFV37_2, peachCalls2.refV37());
        assertEquals(REFV38_2, peachCalls2.refV38());
        assertEquals(ALLELE1_2, peachCalls2.allele1());
        assertEquals(ALLELE2_2, peachCalls2.allele2());
        assertEquals(RSID_2, peachCalls2.rsid());
        assertEquals(VARIANTANNOTATIONV37_2, peachCalls2.variantAnnotationV37());
        assertEquals(FILTERV37_2, peachCalls2.filterV37());
        assertEquals(VARIANTANNOTATIONV38_2, peachCalls2.variantAnnotationV38());
        assertEquals(FILTERV38_2, peachCalls2.filterV38());
        assertEquals(PANEL_VERSION_2, peachCalls2.panelVersion());
        assertEquals(REPO_VERSION_2, peachCalls2.repoVersion());
    }
}