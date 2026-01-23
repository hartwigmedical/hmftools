package com.hartwig.hmftools.finding.clinicaltranscript;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.IOException;
import java.nio.file.Path;

import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;

import org.junit.Test;

public class ClinicalTranscriptFileTest
{

    private static final String CLINICAL_TRANSCRIPT_TSV =
            ClinicalTranscriptsModelTest.class.getClassLoader().getResource("clinicaltranscript/clinical_transcripts.tsv").getPath();

    @Test
    public void canReadClinicalTranscriptsTsvHG37() throws IOException
    {
        ClinicalTranscriptsModel clinicalTranscriptsModel =
                ClinicalTranscriptFile.buildFromTsv(OrangeRefGenomeVersion.V37, Path.of(CLINICAL_TRANSCRIPT_TSV));

        assertEquals("NM_789", clinicalTranscriptsModel.findCanonicalTranscriptForGene("BRCA1"));
        assertEquals("NM_123", clinicalTranscriptsModel.findCanonicalTranscriptForGene("BRCA2"));
    }

    @Test
    public void canReadClinicalTranscriptsTsvHG38() throws IOException
    {
        ClinicalTranscriptsModel clinicalTranscriptsModel =
                ClinicalTranscriptFile.buildFromTsv(OrangeRefGenomeVersion.V38, Path.of(CLINICAL_TRANSCRIPT_TSV));

        assertNull(clinicalTranscriptsModel.findCanonicalTranscriptForGene("BRCA1"));
        assertNull(clinicalTranscriptsModel.findCanonicalTranscriptForGene("BRCA2"));
    }
}