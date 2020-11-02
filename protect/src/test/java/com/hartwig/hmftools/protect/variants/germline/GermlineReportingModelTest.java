package com.hartwig.hmftools.protect.variants.germline;

import static com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel.NO_REPORTING;
import static com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION;
import static com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineReportingModelTest {

    @Test
    public void testBehaviour() {
        Map<String, GermlineReporting> notifyMap = Maps.newHashMap();
        GermlineReporting germlineReportingTrue = ImmutableGermlineReporting.builder()
                .notifyClinicalGeneticus(true)
                .condition(ConditionReportingVariant.BIALLELIC)
                .variant("p.Arg876Cys")
                .build();
        GermlineReporting germlineReportingFalse = ImmutableGermlineReporting.builder()
                .notifyClinicalGeneticus(false)
                .condition(ConditionReportingVariant.MONOALLELIC)
                .variant("")
                .build();
        notifyMap.put("Report", germlineReportingFalse);
        notifyMap.put("Notify", germlineReportingTrue);

        GermlineReportingModel victim = new GermlineReportingModel(notifyMap);
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