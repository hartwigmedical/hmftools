package com.hartwig.hmftools.finding.clinicalrelevantgenecopynumber;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.nio.file.Path;

import com.google.common.io.Resources;

import org.junit.Test;

public class ClinicalRelevantGeneCopyNumberFileTest
{

    private static final String RELEVANT_CLINICAL_GENE_COPY_NUMBERS_TSV =
            Resources.getResource("clinicalrelevantgenecopynumbers/clinical_relevant_gene_copy_numbers.tsv").getPath();

    @Test
    public void testCanReadClinicalRelevantGeneCopyNumbersTsv() throws IOException
    {
        ClinicalRelevantGeneCopyNumberModel clinicalCopyNumberModel =
                ClinicalRelevantGeneCopyNumberFile.buildFromTsv(Path.of(RELEVANT_CLINICAL_GENE_COPY_NUMBERS_TSV));
        assertNotNull(clinicalCopyNumberModel);
    }
}