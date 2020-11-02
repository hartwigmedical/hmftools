package com.hartwig.hmftools.patientreporter.cfreport.data;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.actionability.WithEvent;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.protect.variants.ImmutableReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.ReportableVariantSource;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class EventFilterTest {

    @Test
    public void filterWorks() {
        String somaticEventString = "somatic";
        String germlineEventString = "germline";

        ReportableVariant somaticVariant = variant(ReportableVariantSource.SOMATIC, somaticEventString);
        ReportableVariant germlineVariant = variant(ReportableVariantSource.GERMLINE, germlineEventString);

        // Extra space is added because of protein impact which is missing in this test.
        WithEvent somaticEvent = new DummyEvent(somaticEventString + " ");
        WithEvent germlineEvent = new DummyEvent(germlineEventString + " ");

        assertEquals(2,
                EventFilter.removeEvidenceOnFilteredGermlineVariants(Lists.newArrayList(somaticEvent, germlineEvent),
                        Lists.newArrayList(somaticVariant, germlineVariant),
                        LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION).size());

        assertEquals(1,
                EventFilter.removeEvidenceOnFilteredGermlineVariants(Lists.newArrayList(somaticEvent, germlineEvent),
                        Lists.newArrayList(somaticVariant, germlineVariant),
                        LimsGermlineReportingLevel.NO_REPORTING).size());

        assertEquals(2,
                EventFilter.removeEvidenceOnFilteredGermlineVariants(Lists.newArrayList(somaticEvent, somaticEvent),
                        Lists.newArrayList(somaticVariant, germlineVariant),
                        LimsGermlineReportingLevel.NO_REPORTING).size());

        assertEquals(2,
                EventFilter.removeEvidenceOnFilteredGermlineVariants(Lists.newArrayList(somaticEvent, germlineEvent),
                        Lists.newArrayList(somaticVariant),
                        LimsGermlineReportingLevel.NO_REPORTING).size());
    }

    private static ReportableVariant variant(@NotNull ReportableVariantSource source, @NotNull String event) {
        return ImmutableReportableVariant.builder()
                .source(source)
                .gene(event)
                .chromosome(Strings.EMPTY)
                .position(1)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .alleleReadCount(0)
                .totalReadCount(0)
                .alleleCopyNumber(0D)
                .totalCopyNumber(0D)
                .hotspot(Hotspot.NON_HOTSPOT)
                .driverLikelihood(1D)
                .clonalLikelihood(1D)
                .biallelic(false)
                .build();
    }

    private static class DummyEvent implements WithEvent {

        @NotNull
        private final String event;

        public DummyEvent(@NotNull final String event) {
            this.event = event;
        }

        @NotNull
        @Override
        public String event() {
            return event;
        }
    }

}