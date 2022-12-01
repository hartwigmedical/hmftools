package com.hartwig.hmftools.orange.algo.interpretation;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.orange.TestOrangeReportFactory;

import org.junit.Test;

public class OrangeReportLimiterTest {

    @Test
    public void canLimitAllListsToOne() {
        assertNotNull(OrangeReportLimiter.limitAllListsToMaxOne(TestOrangeReportFactory.createMinimalTestReport()));

        assertNotNull(OrangeReportLimiter.limitAllListsToMaxOne(TestOrangeReportFactory.createProperTestReport()));
    }
}