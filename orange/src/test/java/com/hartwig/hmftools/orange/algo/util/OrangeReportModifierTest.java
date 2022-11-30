package com.hartwig.hmftools.orange.algo.util;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.orange.TestOrangeReportFactory;

import org.junit.Test;

public class OrangeReportModifierTest {

    @Test
    public void canLimitAllListsToOne() {
        assertNotNull(OrangeReportModifier.limitAllListsToMaxOne(TestOrangeReportFactory.createMinimalTestReport()));

        assertNotNull(OrangeReportModifier.limitAllListsToMaxOne(TestOrangeReportFactory.createProperTestReport()));
    }
}