package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.protect.EvidenceType;
import com.hartwig.hmftools.common.protect.ImmutableKnowledgebaseSource;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantSource;
import com.hartwig.hmftools.common.variant.ReportableVariantTestFactory;
import com.hartwig.serve.datamodel.EvidenceDirection;
import com.hartwig.serve.datamodel.EvidenceLevel;
import com.hartwig.serve.datamodel.ImmutableTreatment;
import com.hartwig.serve.datamodel.Knowledgebase;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ConsentFilterFunctionsTest {

    @Test
    public void canFilterVariantsForGermlineConsent() {
        ReportableVariant somaticVariant = createTestReportableVariantBuilder().source(ReportableVariantSource.SOMATIC).build();
        ReportableVariant germlineVariant = createTestReportableVariantBuilder().source(ReportableVariantSource.GERMLINE).build();

        Map<ReportableVariant, Boolean> notifyGermlineVariants = Maps.newHashMap();
        notifyGermlineVariants.put(somaticVariant, false);
        notifyGermlineVariants.put(germlineVariant, true);

        List<ReportableVariantWithNotify> variantsWithNotify =
                ConsentFilterFunctions.filterVariants(Lists.newArrayList(somaticVariant, germlineVariant),
                        notifyGermlineVariants,
                        LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION);
        assertEquals(2, variantsWithNotify.size());
        assertEquals(1, variantsWithNotify.stream().filter(x -> x.variant().source() == ReportableVariantSource.GERMLINE).count());
        assertEquals(1, variantsWithNotify.stream().filter(x -> x.notifyVariant()).count());

        Map<ReportableVariant, Boolean> noNotifyGermlineVariants = Maps.newHashMap();
        noNotifyGermlineVariants.put(somaticVariant, false);
        noNotifyGermlineVariants.put(germlineVariant, false);

        List<ReportableVariantWithNotify> variantsWithoutNotify =
                ConsentFilterFunctions.filterVariants(Lists.newArrayList(somaticVariant, germlineVariant),
                        noNotifyGermlineVariants,
                        LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION);
        assertEquals(2, variantsWithoutNotify.size());
        assertEquals(0, variantsWithoutNotify.stream().filter(x -> x.variant().source() == ReportableVariantSource.GERMLINE).count());
        assertEquals(0, variantsWithoutNotify.stream().filter(x -> x.notifyVariant()).count());

        List<ReportableVariantWithNotify> noGermlineReporting =
                ConsentFilterFunctions.filterVariants(Lists.newArrayList(somaticVariant, germlineVariant),
                        notifyGermlineVariants,
                        LimsGermlineReportingLevel.NO_REPORTING);
        assertEquals(1, noGermlineReporting.size());
        assertEquals(0, noGermlineReporting.stream().filter(x -> x.variant().source() == ReportableVariantSource.GERMLINE).count());
        assertEquals(0, noGermlineReporting.stream().filter(x -> x.notifyVariant()).count());
    }

    @Test
    public void canFilterEvidenceForGermlineConsent() {
        ProtectEvidence evidence = ProtectTestFactory.builder()
                .event("HR deficiency")
                .germline(true)
                .reported(true)
                .treatment(ImmutableTreatment.builder()
                        .treament("TryMe")
                        .sourceRelevantTreatmentApproaches(Sets.newHashSet())
                        .relevantTreatmentApproaches(Sets.newHashSet())
                        .build())
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(ImmutableKnowledgebaseSource.builder()
                        .name(Knowledgebase.ICLUSION)
                        .sourceEvent(Strings.EMPTY)
                        .sourceUrls(Sets.newHashSet())
                        .evidenceType(EvidenceType.AMPLIFICATION)
                        .build()))
                .build();

        List<ProtectEvidence> withNotify = ConsentFilterFunctions.filterEvidenceForGermlineConsent(Lists.newArrayList(evidence),
                LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION);
        assertEquals(1, withNotify.size());
        assertEquals(0, withNotify.stream().filter(x -> x.germline()).count());

        assertEquals(0,
                ConsentFilterFunctions.filterEvidenceForGermlineConsent(Lists.newArrayList(evidence),
                        LimsGermlineReportingLevel.NO_REPORTING).size());
    }

    @NotNull
    private static ImmutableReportableVariant.Builder createTestReportableVariantBuilder() {
        return ImmutableReportableVariant.builder().from(ReportableVariantTestFactory.create());
    }
}