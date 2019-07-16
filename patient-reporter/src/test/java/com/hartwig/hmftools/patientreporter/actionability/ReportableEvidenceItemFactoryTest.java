package com.hartwig.hmftools.patientreporter.actionability;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;
import com.hartwig.hmftools.patientreporter.variants.somatic.SomaticVariantAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ReportableEvidenceItemFactoryTest {
    private static final Logger LOGGER = LogManager.getLogger(ReportableEvidenceItemFactoryTest.class);

    @Test
    public void higherLevelWorks() {
        EvidenceItem item1 = builder().level(EvidenceLevel.LEVEL_A).build();
        EvidenceItem item2 = builder().level(EvidenceLevel.LEVEL_B).build();

        assertTrue(ReportableEvidenceItemFactory.hasHigherEvidence(item1, item2));
        assertFalse(ReportableEvidenceItemFactory.hasHigherEvidence(item2, item1));
        assertFalse(ReportableEvidenceItemFactory.hasHigherEvidence(item1, item1));
    }

    @Test
    public void selectHighDriverVariantForEvidence() {
        EnrichedSomaticVariant variant1 =
                PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene("TP53").build();

        Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant = Maps.newHashMap();

        evidencePerVariant.put(variant1, Lists.newArrayList(builder().level(EvidenceLevel.LEVEL_A).isOnLabel(true).build()));

        List<DriverCatalog> catalog = Lists.newArrayList(ImmutableDriverCatalog.builder()
                .gene("TP53")
                .driverLikelihood(0.9)
                .category(DriverCategory.TSG)
                .driver(DriverType.MUTATION)
                .likelihoodMethod(LikelihoodMethod.HOTSPOT)
                .dndsLikelihood(0.4)
                .missense(10)
                .nonsense(20)
                .splice(10)
                .inframe(40)
                .frameshift(50)
                .biallelic(true)
                .build());

        List<EvidenceItem> evidenceForHighDrivers =
                ReportableEvidenceItemFactory.reportableFlatListDriversOnly(evidencePerVariant, catalog);

        assertEquals(1, evidenceForHighDrivers.size());
    }

    @Test
    public void filterNoneHighDriverVariantsOutForEvidence() {
        EnrichedSomaticVariant variant1 =
                PatientReporterTestFactory.createTestEnrichedSomaticVariantBuilder().gene("KRAS").build();

        Map<EnrichedSomaticVariant, List<EvidenceItem>> evidencePerVariant = Maps.newHashMap();

        evidencePerVariant.put(variant1, Lists.newArrayList(builder().level(EvidenceLevel.LEVEL_A).isOnLabel(true).build()));

        List<DriverCatalog> catalog = Lists.newArrayList(ImmutableDriverCatalog.builder()
                .gene("KRAS")
                .driverLikelihood(0.6)
                .category(DriverCategory.ONCO)
                .driver(DriverType.MUTATION)
                .likelihoodMethod(LikelihoodMethod.HOTSPOT)
                .dndsLikelihood(0.4)
                .missense(10)
                .nonsense(20)
                .splice(10)
                .inframe(40)
                .frameshift(50)
                .biallelic(true)
                .build());

        List<EvidenceItem> evidenceForHighDrivers =
                ReportableEvidenceItemFactory.reportableFlatListDriversOnly(evidencePerVariant, catalog);

        assertEquals(0, evidenceForHighDrivers.size());

    }

    @Test
    public void canSelectHighest() {
        EvidenceItem item1 = builder().level(EvidenceLevel.LEVEL_A).isOnLabel(true).build();
        EvidenceItem item2 = builder().level(EvidenceLevel.LEVEL_B).isOnLabel(false).build();
        EvidenceItem item3 = builder().level(EvidenceLevel.LEVEL_C).isOnLabel(false).build();

        EvidenceItem highestOffLabel = ReportableEvidenceItemFactory.highestOffLabel(Lists.newArrayList(item1, item2, item3));
        assertNotNull(highestOffLabel);
        assertEquals(EvidenceLevel.LEVEL_B, highestOffLabel.level());

        EvidenceItem highestOnLabel = ReportableEvidenceItemFactory.highestOnLabel(Lists.newArrayList(item1, item2, item3));
        assertNotNull(highestOnLabel);
        assertEquals(EvidenceLevel.LEVEL_A, highestOnLabel.level());

        EvidenceItem onLabelOnly = ReportableEvidenceItemFactory.highestOffLabel(Lists.newArrayList(item1));
        assertNull(onLabelOnly);
    }

    @Test
    public void hasKnownLevel() {
        EvidenceItem item1 = builder().level(EvidenceLevel.LEVEL_A).isOnLabel(true).build();
        EvidenceItem item2 = builder().level(EvidenceLevel.LEVEL_B).isOnLabel(true).build();
        EvidenceItem item3 = builder().level(EvidenceLevel.LEVEL_C).isOnLabel(false).build();
        EvidenceItem item4 = builder().level(EvidenceLevel.LEVEL_D).isOnLabel(false).build();
        EvidenceItem item5 = builder().level(EvidenceLevel.LEVEL_E).isOnLabel(false).build();

        assertTrue(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item1));
        assertTrue(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item2));
        assertFalse(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item3));
        assertFalse(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item4));
        assertFalse(ReportableEvidenceItemFactory.hasReportableEvidenceLevel(item5));
    }

    @NotNull
    private static ImmutableEvidenceItem.Builder builder() {
        return ImmutableEvidenceItem.builder()
                .event(Strings.EMPTY)
                .source(ActionabilitySource.CIVIC)
                .reference(Strings.EMPTY)
                .drug(Strings.EMPTY)
                .drugsType(Strings.EMPTY)
                .level(EvidenceLevel.LEVEL_A)
                .response(Strings.EMPTY)
                .isOnLabel(false)
                .cancerType(Strings.EMPTY)
                .scope(EvidenceScope.SPECIFIC);
    }
}