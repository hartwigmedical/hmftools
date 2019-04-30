package com.hartwig.hmftools.patientreporter.variants.germline;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientreporter.PatientReporterTestUtil;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FilterGermlineVariantsTest {

    @Test
    public void checkForGermlineGenesReportedONCO() throws IOException {
        List<GermlineVariant> germlineVariants = createTestGermlineVariantsONCOGene();
        Map<String, Boolean> germlineGenesReporting = PatientReporterTestUtil.testGermlineModel().germlineGenes();
        Map<String, DriverCategory> driverCategoryMapMatch = Maps.newHashMap();
        driverCategoryMapMatch.put("KIT", DriverCategory.ONCO);

        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();
        List<EnrichedSomaticVariant> variants = Lists.newArrayList();

        List<GermlineVariant> filteredGermlineVariantMatch = FilterGermlineVariants.filterGermlineVariantsForReporting(germlineVariants,
                germlineGenesReporting,
                driverCategoryMapMatch,
                geneCopyNumbers,
                variants,
                LimsGermlineReportingChoice.ACTIONABLE_ONLY);
        assertEquals(filteredGermlineVariantMatch.size(), 1);

        Map<String, DriverCategory> driverCategoryMapNotMatch = Maps.newHashMap();
        driverCategoryMapNotMatch.put("AAAA", DriverCategory.ONCO);

        List<GermlineVariant> filteredGermlineVariantNonMatch = FilterGermlineVariants.filterGermlineVariantsForReporting(germlineVariants,
                germlineGenesReporting,
                driverCategoryMapNotMatch,
                geneCopyNumbers,
                variants,
                LimsGermlineReportingChoice.ACTIONABLE_ONLY);
        assertEquals(filteredGermlineVariantNonMatch.size(), 0);
    }

    @Test
    public void checkForGermlineGenesReportedTSG() throws IOException {
        Map<String, Boolean> germlineGenesReporting = PatientReporterTestUtil.testGermlineModel().germlineGenes();
        Map<String, DriverCategory> driverCategoryMapMatch = Maps.newHashMap();
        driverCategoryMapMatch.put("BRCA2", DriverCategory.TSG);
        List<GermlineVariant> germlineVariantsMatch = createTestGermlineVariantsTSGGene(true);
        List<GeneCopyNumber> geneCopyNumbersMatch = Lists.newArrayList(createTestCopyNumberBuilder(1).build());
        List<EnrichedSomaticVariant> variantsMatch = Lists.newArrayList(createEnriched("BRCA2").build());
        List<GermlineVariant> filteredGermlineVariantMatch =
                FilterGermlineVariants.filterGermlineVariantsForReporting(germlineVariantsMatch,
                        germlineGenesReporting,
                        driverCategoryMapMatch,
                        geneCopyNumbersMatch,
                        variantsMatch,
                        LimsGermlineReportingChoice.ACTIONABLE_ONLY);
        assertEquals(filteredGermlineVariantMatch.size(), 1); // all three options matched

        List<GermlineVariant> germlineVariantsNonMatchBiallelic = createTestGermlineVariantsTSGGene(false);
        List<GeneCopyNumber> geneCopyNumbersNonMatchBiallelic = Lists.newArrayList(createTestCopyNumberBuilder(1).build());
        Map<String, DriverCategory> driverCategoryMapNonMatchBilallelic = Maps.newHashMap();
        driverCategoryMapNonMatchBilallelic.put("BRCA2", DriverCategory.TSG);
        List<EnrichedSomaticVariant> variantsNonMatchBiallelic = Lists.newArrayList(createEnriched("BRCA2").build());
        List<GermlineVariant> filteredGermlineVariantNonMatchBiallelic = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsNonMatchBiallelic,
                germlineGenesReporting,
                driverCategoryMapNonMatchBilallelic,
                geneCopyNumbersNonMatchBiallelic,
                variantsNonMatchBiallelic,
                LimsGermlineReportingChoice.ACTIONABLE_ONLY);
        assertEquals(filteredGermlineVariantNonMatchBiallelic.size(), 1); // option biallelic failed

        List<GermlineVariant> germlineVariantsNonMatchVariant = createTestGermlineVariantsTSGGene(true);
        List<GeneCopyNumber> geneCopyNumbersNonMatchVariant = Lists.newArrayList(createTestCopyNumberBuilder(1).build());
        Map<String, DriverCategory> driverCategoryMapNonMatchVariant = Maps.newHashMap();
        driverCategoryMapNonMatchVariant.put("BRCA2", DriverCategory.TSG);
        List<EnrichedSomaticVariant> variantsNonMatchVariant = Lists.newArrayList(createEnriched("ATM").build());
        List<GermlineVariant> filteredGermlineVariantNonMatchVariant = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsNonMatchVariant,
                germlineGenesReporting,
                driverCategoryMapNonMatchVariant,
                geneCopyNumbersNonMatchVariant,
                variantsNonMatchVariant,
                LimsGermlineReportingChoice.ACTIONABLE_ONLY);
        assertEquals(filteredGermlineVariantNonMatchVariant.size(), 1); // option variant failed

        List<GermlineVariant> germlineVariantsNonMatchCopy = createTestGermlineVariantsTSGGene(true);
        List<GeneCopyNumber> geneCopyNumbersNonMatchCopy = Lists.newArrayList(createTestCopyNumberBuilder(0).build());
        Map<String, DriverCategory> driverCategoryMapNonMatchCopy = Maps.newHashMap();
        driverCategoryMapNonMatchCopy.put("BRCA2", DriverCategory.TSG);
        List<EnrichedSomaticVariant> variantsNonMatchCopy = Lists.newArrayList(createEnriched("BRCA2").build());
        List<GermlineVariant> filteredGermlineVariantNonMatchCopy = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsNonMatchCopy,
                germlineGenesReporting,
                driverCategoryMapNonMatchCopy,
                geneCopyNumbersNonMatchCopy,
                variantsNonMatchCopy,
                LimsGermlineReportingChoice.ACTIONABLE_ONLY);
        assertEquals(filteredGermlineVariantNonMatchCopy.size(), 1); // option copy number failed

        List<GermlineVariant> germlineVariantsNonMatch = createTestGermlineVariantsTSGGene(false);
        List<GeneCopyNumber> geneCopyNumbersNonMatch = Lists.newArrayList(createTestCopyNumberBuilder(0).build());
        Map<String, DriverCategory> driverCategoryMapNonMatch = Maps.newHashMap();
        driverCategoryMapNonMatch.put("BRCA2", DriverCategory.TSG);
        List<EnrichedSomaticVariant> variantsNonMatch = Lists.newArrayList(createEnriched("ATM").build());
        List<GermlineVariant> filteredGermlineVariantNonMatch = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsNonMatch,
                germlineGenesReporting,
                driverCategoryMapNonMatch,
                geneCopyNumbersNonMatch,
                variantsNonMatch,
                LimsGermlineReportingChoice.ACTIONABLE_ONLY);
        assertEquals(filteredGermlineVariantNonMatch.size(), 0); // all option failed

        List<GermlineVariant> germlineVariantsOptionBaillelic = createTestGermlineVariantsTSGGene(true);
        List<GeneCopyNumber> geneCopyNumbersOptionBaillelic = Lists.newArrayList(createTestCopyNumberBuilder(0).build());
        Map<String, DriverCategory> driverCategoryMapOptionBaillelic = Maps.newHashMap();
        driverCategoryMapOptionBaillelic.put("BRCA2", DriverCategory.TSG);
        List<EnrichedSomaticVariant> variantsOptionBaillelic = Lists.newArrayList(createEnriched("ATM").build());
        List<GermlineVariant> filteredGermlineVariantOptionBaillelic = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsOptionBaillelic,
                germlineGenesReporting,
                driverCategoryMapOptionBaillelic,
                geneCopyNumbersOptionBaillelic,
                variantsOptionBaillelic,
                LimsGermlineReportingChoice.ACTIONABLE_ONLY);
        assertEquals(filteredGermlineVariantOptionBaillelic.size(), 1); // only match biallelic

        List<GermlineVariant> germlineVariantsOptionVariant = createTestGermlineVariantsTSGGene(false);
        List<GeneCopyNumber> geneCopyNumbersOptionVariant = Lists.newArrayList(createTestCopyNumberBuilder(0).build());
        Map<String, DriverCategory> driverCategoryMapOptionVaraint = Maps.newHashMap();
        driverCategoryMapOptionVaraint.put("BRCA2", DriverCategory.TSG);
        List<EnrichedSomaticVariant> variantsOptionVariant = Lists.newArrayList(createEnriched("BRCA2").build());
        List<GermlineVariant> filteredGermlineVariantOptionVariant = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsOptionVariant,
                germlineGenesReporting,
                driverCategoryMapOptionVaraint,
                geneCopyNumbersOptionVariant,
                variantsOptionVariant,
                LimsGermlineReportingChoice.ACTIONABLE_ONLY);
        assertEquals(filteredGermlineVariantOptionVariant.size(), 1); // only match variant

        List<GermlineVariant> germlineVariantsOptionCopyNumber = createTestGermlineVariantsTSGGene(false);
        List<GeneCopyNumber> geneCopyNumbersCopyNumber = Lists.newArrayList(createTestCopyNumberBuilder(1).build());
        Map<String, DriverCategory> driverCategoryMapOptionCopyNumber = Maps.newHashMap();
        driverCategoryMapOptionCopyNumber.put("BRCA2", DriverCategory.TSG);
        List<EnrichedSomaticVariant> variantsOptionCopyNumber = Lists.newArrayList(createEnriched("ATM").build());
        List<GermlineVariant> filteredGermlineVariantOptionCopyNumber = FilterGermlineVariants.filterGermlineVariantsForReporting(
                germlineVariantsOptionCopyNumber,
                germlineGenesReporting,
                driverCategoryMapOptionCopyNumber,
                geneCopyNumbersCopyNumber,
                variantsOptionCopyNumber,
                LimsGermlineReportingChoice.ACTIONABLE_ONLY);
        assertEquals(filteredGermlineVariantOptionCopyNumber.size(), 1); // only match copy number

    }

    @NotNull
    private static List<GermlineVariant> createTestGermlineVariantsONCOGene() {
        List<GermlineVariant> germlineVariants = Lists.newArrayList();

        germlineVariants.add(ImmutableGermlineVariant.builder()
                .passFilter(true)
                .gene("KIT")
                .hgvsCodingImpact(Strings.EMPTY)
                .hgvsProteinImpact(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .germlineStatus(Strings.EMPTY)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAllelePloidy(0)
                .biallelic(false)
                .build());

        return germlineVariants;
    }

    @NotNull
    private static List<GermlineVariant> createTestGermlineVariantsTSGGene(boolean biallelicFilter) {
        List<GermlineVariant> germlineVariants = Lists.newArrayList();

        germlineVariants.add(ImmutableGermlineVariant.builder()
                .passFilter(true)
                .gene("BRCA2")
                .hgvsCodingImpact(Strings.EMPTY)
                .hgvsProteinImpact(Strings.EMPTY)
                .totalReadCount(0)
                .alleleReadCount(0)
                .germlineStatus(Strings.EMPTY)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAllelePloidy(1D)
                .biallelic(biallelicFilter)
                .build());

        return germlineVariants;
    }

    @NotNull
    private static ImmutableGeneCopyNumber.Builder createTestCopyNumberBuilder(int minNumberFilter) {
        return ImmutableGeneCopyNumber.builder()
                .start(0)
                .end(0)
                .gene("BRCA2")
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .minRegionStart(0)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEnd(0)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegions(0)
                .germlineHet2HomRegions(0)
                .germlineHomRegions(0)
                .somaticRegions(0)
                .minCopyNumber(minNumberFilter)
                .maxCopyNumber(0)
                .transcriptID(Strings.EMPTY)
                .transcriptVersion(0)
                .minMinorAllelePloidy(0);
    }

    @NotNull
    private static ImmutableEnrichedSomaticVariant.Builder createEnriched(@NotNull String gene) {
        return ImmutableEnrichedSomaticVariant.builder()
                .from(create(gene).build())
                .trinucleotideContext(Strings.EMPTY)
                .highConfidenceRegion(false)
                .microhomology(Strings.EMPTY)
                .repeatSequence(Strings.EMPTY)
                .repeatCount(0)
                .clonality(Clonality.UNKNOWN);
    }

    @NotNull
    private static ImmutableSomaticVariantImpl.Builder create(@NotNull String gene) {
        return ImmutableSomaticVariantImpl.builder()
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .type(VariantType.UNDEFINED)
                .filter("PASS")
                .totalReadCount(0)
                .alleleReadCount(0)
                .gene(gene)
                .genesEffected(0)
                .worstEffect(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.NONE)
                .worstEffectTranscript(Strings.EMPTY)
                .canonicalEffect(Strings.EMPTY)
                .canonicalCodingEffect(CodingEffect.UNDEFINED)
                .canonicalHgvsCodingImpact(Strings.EMPTY)
                .canonicalHgvsProteinImpact(Strings.EMPTY)
                .hotspot(Hotspot.NON_HOTSPOT)
                .recovered(false)
                .adjustedCopyNumber(0)
                .adjustedVAF(0)
                .minorAllelePloidy(0)
                .germlineStatus(GermlineStatus.UNKNOWN)
                .ploidy(0)
                .mappability(0);
    }

}