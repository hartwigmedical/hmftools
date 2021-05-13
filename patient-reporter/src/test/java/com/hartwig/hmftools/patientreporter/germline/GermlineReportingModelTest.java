package com.hartwig.hmftools.patientreporter.germline;

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
                .notifyClinicalGeneticist(GermlineCondition.ALWAYS)
                .conditionFilter(null)
                .build();

        GermlineReportingEntry germlineReportingFalse = ImmutableGermlineReportingEntry.builder()
                .gene(reportGene)
                .notifyClinicalGeneticist(GermlineCondition.NEVER)
                .conditionFilter(null)
                .build();

        GermlineReportingModel model = new GermlineReportingModel(Lists.newArrayList(germlineReportingTrue, germlineReportingFalse));
        assertNotNull(model.entryForGene(notifyGene));
        assertNotNull(model.entryForGene(reportGene));
        assertNull(model.entryForGene("Other"));

        assertTrue(model.notifyAboutGene(notifyGene, REPORT_WITH_NOTIFICATION));
        assertFalse(model.notifyAboutGene(notifyGene, REPORT_WITHOUT_NOTIFICATION));
        assertFalse(model.notifyAboutGene(notifyGene, NO_REPORTING));

        assertFalse(model.notifyAboutGene(reportGene, REPORT_WITHOUT_NOTIFICATION));
        assertFalse(model.notifyAboutGene(reportGene, REPORT_WITHOUT_NOTIFICATION));
        assertFalse(model.notifyAboutGene(reportGene, NO_REPORTING));

        assertFalse(model.notifyAboutGene("Other", REPORT_WITH_NOTIFICATION));
    }
}