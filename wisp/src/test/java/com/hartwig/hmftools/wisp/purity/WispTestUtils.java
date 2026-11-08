package com.hartwig.hmftools.wisp.purity;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.wisp.ResultsWriter;
import com.hartwig.hmftools.wisp.SampleData;
import com.hartwig.hmftools.wisp.WispConfig;

public class WispTestUtils
{
    public static final WispConfig TEST_CONFIG = new WispConfig(new MockRefGenome());

    public static final String TEST_PATIENT_ID = "PATIENT_001";
    public static final String TEST_TUMOR_ID = "TUMOR_ID_001";
    public static final String TEST_REFERENCE_ID = "REF_ID_001";
    public static final String TEST_SAMPLE_ID = "SAMPLE_ID_001";

    public static final SampleData TEST_SAMPLE = new SampleData(
            TEST_PATIENT_ID, TEST_TUMOR_ID, TEST_REFERENCE_ID, Lists.newArrayList(TEST_SAMPLE_ID),
            "", false, "");

    public static final ResultsWriter TEST_RESULTS_WRITER = new ResultsWriter(TEST_CONFIG);


}
