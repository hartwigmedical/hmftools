package com.hartwig.hmftools.protect.variants.germline;

import static com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel.NO_REPORTING;
import static com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION;
import static com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;

import org.junit.Test;

public class GermlineReportingModelTest {

    @Test
    public void testBehaviour() {
        GermlineReportingEntry germlineReportingTrue = ImmutableGermlineReportingEntry.builder()
                .gene("Notify")
                .notifyClinicalGeneticist(true)
                .reportBiallelicOnly(true)
                .reportableHgvsProtein("p.Arg876Cys")
                .build();
        GermlineReportingEntry germlineReportingFalse = ImmutableGermlineReportingEntry.builder()
                .gene("Report")
                .notifyClinicalGeneticist(false)
                .reportBiallelicOnly(false)
                .reportableHgvsProtein(null)
                .build();

        GermlineReportingModel victim = new GermlineReportingModel(Lists.newArrayList(germlineReportingTrue, germlineReportingFalse));
        assertTrue(victim.notifyAboutGene(REPORT_WITH_NOTIFICATION, "Notify"));
        assertFalse(victim.notifyAboutGene(REPORT_WITHOUT_NOTIFICATION, "Notify"));
        assertFalse(victim.notifyAboutGene(NO_REPORTING, "Notify"));

        assertTrue(victim.reportableGermlineGenes().contains("Report"));
        assertTrue(victim.reportableGermlineGenes().contains("Notify"));

        assertTrue(victim.monoallelicGenesReportable().contains("Report"));

        assertTrue(victim.reportableSpecificVariants().containsKey("Notify"));
        assertTrue(victim.reportableSpecificVariants().containsValue("p.Arg876Cys"));

        assertFalse(victim.notifiableGenes(REPORT_WITH_NOTIFICATION).contains("Report"));
        assertTrue(victim.notifiableGenes(REPORT_WITH_NOTIFICATION).contains("Notify"));
        assertFalse(victim.notifiableGenes(REPORT_WITHOUT_NOTIFICATION).contains("Notify"));

        assertFalse(victim.notifyAboutGene(REPORT_WITH_NOTIFICATION, "DoesNotExist"));
    }
}