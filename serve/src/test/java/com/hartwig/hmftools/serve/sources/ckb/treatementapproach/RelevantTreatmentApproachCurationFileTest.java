package com.hartwig.hmftools.serve.sources.ckb.treatementapproach;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;

import org.junit.Test;

public class RelevantTreatmentApproachCurationFileTest {

    private static final String TEST_CKB_DRUGCLASS_CURATION_TSV =
            Resources.getResource("ckb_curation/ckb_treatment_approach_curation.tsv").getPath();

    @Test
    public void canReadCkbDrugClassCurationTsv() throws IOException {
        Map<RelevantTreatmentApprochCurationEntryKey, RelevantTreatmentApprochCurationEntry> treatmentApproachMap =
                RelevantTreatmentApproachCurationFile.read(TEST_CKB_DRUGCLASS_CURATION_TSV);
        assertEquals(4, treatmentApproachMap.size());
    }
}