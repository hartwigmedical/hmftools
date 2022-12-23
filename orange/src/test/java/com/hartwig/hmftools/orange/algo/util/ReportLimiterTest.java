package com.hartwig.hmftools.orange.algo.util;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.orange.TestOrangeReportFactory;

import org.junit.Test;

public class ReportLimiterTest {

    @Test
    public void canLimitAllListsToOne() {
        assertNotNull(ReportLimiter.limitAllListsToMaxOne(TestOrangeReportFactory.createMinimalTestReport()));

        assertNotNull(ReportLimiter.limitAllListsToMaxOne(TestOrangeReportFactory.createProperTestReport()));
    }
}