package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.protect.purple.ImmutableReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantSource;

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

        List<ReportableVariant> variantsWithNotify =
                ConsentFilterFunctions.filterVariants(Lists.newArrayList(somaticVariant, germlineVariant),
                        notifyGermlineVariants,
                        LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION);
        assertEquals(2, variantsWithNotify.size());
        assertEquals(1, variantsWithNotify.stream().filter(x -> x.source() == ReportableVariantSource.GERMLINE).count());

        Map<ReportableVariant, Boolean> noNotifyGermlineVariants = Maps.newHashMap();
        noNotifyGermlineVariants.put(somaticVariant, false);
        noNotifyGermlineVariants.put(germlineVariant, false);

        List<ReportableVariant> variantsWithoutNotify =
                ConsentFilterFunctions.filterVariants(Lists.newArrayList(somaticVariant, germlineVariant),
                        noNotifyGermlineVariants,
                        LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION);
        assertEquals(2, variantsWithoutNotify.size());
        assertEquals(0, variantsWithoutNotify.stream().filter(x -> x.source() == ReportableVariantSource.GERMLINE).count());

        List<ReportableVariant> noGermlineReporting =
                ConsentFilterFunctions.filterVariants(Lists.newArrayList(somaticVariant, germlineVariant),
                        notifyGermlineVariants,
                        LimsGermlineReportingLevel.NO_REPORTING);
        assertEquals(1, noGermlineReporting.size());
    }

    @Test
    public void canFilterEvidenceForGermlineConsent() {
        ProtectEvidence evidence = ImmutableProtectEvidence.builder()
                .genomicEvent("HR deficiency signature")
                .germline(true)
                .reported(true)
                .treatment("TryMe")
                .onLabel(true)
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .sources(Sets.newHashSet(Knowledgebase.ICLUSION))
                .build();

        assertEquals(1,
                ConsentFilterFunctions.filterEvidenceForGermlineConsent(Lists.newArrayList(evidence),
                        LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION).size());

        assertEquals(0,
                ConsentFilterFunctions.filterEvidenceForGermlineConsent(Lists.newArrayList(evidence),
                        LimsGermlineReportingLevel.NO_REPORTING).size());

        assertEquals(0,
                ConsentFilterFunctions.filterEvidenceForGermlineConsent(Lists.newArrayList(evidence),
                        LimsGermlineReportingLevel.NO_REPORTING).size());

        assertEquals(0,
                ConsentFilterFunctions.filterEvidenceForGermlineConsent(Lists.newArrayList(evidence),
                        LimsGermlineReportingLevel.NO_REPORTING).size());
    }

    @NotNull
    private static ImmutableReportableVariant.Builder createTestReportableVariantBuilder() {
        return ImmutableReportableVariant.builder()
                .type(VariantType.SNP)
                .gene(Strings.EMPTY)
                .genotypeStatus(GenotypeStatus.UNKNOWN)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .canonicalTranscript(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .totalCopyNumber(0)
                .alleleCopyNumber(0D)
                .hotspot(Hotspot.HOTSPOT)
                .clonalLikelihood(1D)
                .driverLikelihood(0D)
                .biallelic(false);
    }
}