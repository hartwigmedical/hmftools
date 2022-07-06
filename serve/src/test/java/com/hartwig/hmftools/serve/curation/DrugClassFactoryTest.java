package com.hartwig.hmftools.serve.curation;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Map;

import com.google.common.io.Resources;

import org.junit.Test;

public class DrugClassFactoryTest {

    private static final String TEST_CKB_DRUGCLASS_CURATION_TSV =
            Resources.getResource("ckb_curation/ckb_drugclass_curation.tsv").getPath();

    @Test
    public void canReadCkbDrugClassCurationTsv() throws IOException {
        Map<DrugClassKey, DrugClasses> drugClassesMap = DrugClassFactory.read(TEST_CKB_DRUGCLASS_CURATION_TSV);
        assertEquals(3, drugClassesMap.size());
    }
}