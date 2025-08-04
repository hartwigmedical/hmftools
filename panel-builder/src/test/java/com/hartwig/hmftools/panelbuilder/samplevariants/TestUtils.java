package com.hartwig.hmftools.panelbuilder.samplevariants;

import com.hartwig.hmftools.common.test.MockRefGenome;

// TODO

public final class TestUtils
{
    public static final String TEST_SAMPLE_ID = "SAMPLE_ID";
    public static final String TEST_REF_ID = "REF_ID";
    public static final double DEFAULT_QUAL = 1000;

    //                                             1        10        20        30        40        50
    //                                             12345678901234567890123456789012345678901234567890123456789
    public static final String REF_BASES_CHR_1 = "XAACCGGTTACGTAAACCCGGGTTTAACCGGTTACGTAAACCCGGGTTTAAACCCGGGTTT";
    public static final String REF_BASES_CHR_2 = "XAAACCCGGGTTTAACCGGTTACGTAAACCCGGGTTTAACCGGTTACGTACCGGTTACGTT";

    public static final MockRefGenome MOCK_REF_GENOME = new MockRefGenome();

    public static final ProbeConfig TEST_CONFIG = new ProbeConfig(
            1, 0, 0, 0);

}
