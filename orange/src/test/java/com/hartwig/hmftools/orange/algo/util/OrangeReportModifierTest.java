package com.hartwig.hmftools.orange.algo.util;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.orange.OrangeReportTestFactory;

import org.junit.Test;

public class OrangeReportModifierTest {

    @Test
    public void canLimitAllListsToOne() {
        assertNotNull(OrangeReportModifier.limitAllListsToMaxOne(OrangeReportTestFactory.createMinimalTestReport()));

        assertNotNull(OrangeReportModifier.limitAllListsToMaxOne(OrangeReportTestFactory.createProperTestReport()));
    }
}