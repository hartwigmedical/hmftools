package com.hartwig.hmftools.orange.util;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.orange.OrangeReportTestFactory;

import org.junit.Test;

public class OrangeReportModifierTest {

    @Test
    public void canLimitJsonOutput() {
        assertNotNull(OrangeReportModifier.limitAllListsToMaxOne(OrangeReportTestFactory.createProperTestReport()));
    }
}