package com.hartwig.hmftools.protect.variants.germline;

import static org.junit.Assert.*;

import org.junit.Test;

public class ConditionReportingVariantTest {

    @Test
    public void canConvertConditionReportingVariant() {
        assertEquals(ConditionReportingVariant.BIALLELIC_ONLY, ConditionReportingVariant.fromConditionString("Biallelic"));
        assertEquals(ConditionReportingVariant.ALL, ConditionReportingVariant.fromConditionString("Monoallelic"));
        assertEquals(ConditionReportingVariant.UNKNOWN, ConditionReportingVariant.fromConditionString("test"));
    }

}