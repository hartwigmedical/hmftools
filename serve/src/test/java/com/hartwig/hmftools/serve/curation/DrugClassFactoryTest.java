package com.hartwig.hmftools.serve.curation;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.serve.treatementapproach.curation.RelevantTreatmentApproachFactory;
import com.hartwig.hmftools.serve.treatementapproach.curation.RelevantTreatmentApproachKey;
import com.hartwig.hmftools.serve.treatementapproach.curation.RelevantTreatmentApproch;

import org.junit.Test;

public class DrugClassFactoryTest {

    private static final String TEST_CKB_DRUGCLASS_CURATION_TSV =
            Resources.getResource("ckb_curation/ckb_drugclass_curation.tsv").getPath();

    @Test
    public void canReadCkbDrugClassCurationTsv() throws IOException {
        Map<RelevantTreatmentApproachKey, RelevantTreatmentApproch> drugClassesMap = RelevantTreatmentApproachFactory.read(TEST_CKB_DRUGCLASS_CURATION_TSV);
        assertEquals(3, drugClassesMap.size());
    }
}