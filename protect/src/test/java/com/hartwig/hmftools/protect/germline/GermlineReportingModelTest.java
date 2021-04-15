package com.hartwig.hmftools.protect.germline;

import static com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel.NO_REPORTING;
import static com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION;
import static com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;

import org.junit.Test;

public class GermlineReportingModelTest {

    @Test
    public void modelWorksAsExpected() {
        String notifyGene = "Notify";
        String reportGene = "Report";

        GermlineReportingEntry germlineReportingTrue = ImmutableGermlineReportingEntry.builder()
                .gene(notifyGene)
                .notifyClinicalGeneticist(true)
                .exclusiveHgvsProteinFilter(null)
                .build();

        GermlineReportingEntry germlineReportingFalse = ImmutableGermlineReportingEntry.builder()
                .gene(reportGene)
                .notifyClinicalGeneticist(false)
                .exclusiveHgvsProteinFilter(null)
                .build();

        GermlineReportingModel victim = new GermlineReportingModel(Lists.newArrayList(germlineReportingTrue, germlineReportingFalse));
        assertNotNull(victim.entryForGene(notifyGene));
        assertNotNull(victim.entryForGene(reportGene));
        assertNull(victim.entryForGene("Other"));

        assertTrue(victim.notifyAboutGene(notifyGene, REPORT_WITH_NOTIFICATION));
        assertFalse(victim.notifyAboutGene(notifyGene, REPORT_WITHOUT_NOTIFICATION));
        assertFalse(victim.notifyAboutGene(notifyGene, NO_REPORTING));

        assertFalse(victim.notifyAboutGene(reportGene, REPORT_WITH_NOTIFICATION));
        assertFalse(victim.notifyAboutGene(reportGene, REPORT_WITHOUT_NOTIFICATION));
        assertFalse(victim.notifyAboutGene(reportGene, NO_REPORTING));

        assertFalse(victim.notifyAboutGene("Other", REPORT_WITH_NOTIFICATION));
    }
}