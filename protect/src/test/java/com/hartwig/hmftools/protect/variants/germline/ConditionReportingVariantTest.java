package com.hartwig.hmftools.protect.variants.germline;

import static org.junit.Assert.*;

import org.junit.Test;

public class ConditionReportingVariantTest {

    @Test
    public void canConvertConditionReportingVariant() {
        assertEquals(ConditionReportingVariant.BIALLELIC, ConditionReportingVariant.fromConditionString("Biallelic"));
        assertEquals(ConditionReportingVariant.MONOALLELIC, ConditionReportingVariant.fromConditionString("Monoallelic"));
        assertEquals(ConditionReportingVariant.UNKNOWN, ConditionReportingVariant.fromConditionString("test"));
    }

}