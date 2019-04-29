package com.hartwig.hmftools.patientreporter.germline;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
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
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.patientreporter.PatientReporterTestUtil;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FilterGermlineVariantsTest {

    @Test
    public void checkForGermlineGenesReportedONCO() throws IOException{
        List<GermlineVariant> germlineVariants = createTestGermlineVariantsONCOGene();
        Map<String,Boolean> germlineGenesReporting = PatientReporterTestUtil.testGermlineModel().germlineGenes();
        Map<String, DriverCategory> driverCategoryMapMatch = Maps.newHashMap();
        driverCategoryMapMatch.put("BRCA2", DriverCategory.ONCO);

        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();
        String sampleId = "CPCT02990001T";
        List<EnrichedSomaticVariant> variants = Lists.newArrayList();

        List<GermlineVariant> filteredGermlineVariantMatch = FilterGermlineVariants.filteringReportedGermlineVariant(germlineVariants,
                germlineGenesReporting,
                driverCategoryMapMatch,
                geneCopyNumbers,
                sampleId,
                variants);
        assertEquals(filteredGermlineVariantMatch.size(), 1);

        Map<String, DriverCategory> driverCategoryMapNotMatch = Maps.newHashMap();
        driverCategoryMapNotMatch.put("BRCA1", DriverCategory.ONCO);

        List<GermlineVariant> filteredGermlineVariantNonMatch = FilterGermlineVariants.filteringReportedGermlineVariant(germlineVariants,
                germlineGenesReporting,
                driverCategoryMapNotMatch,
                geneCopyNumbers,
                sampleId,
                variants);
        assertEquals(filteredGermlineVariantNonMatch.size(), 0);
    }

    @Test
    public void checkForGermlineGenesReportedTSG() throws IOException {
        Map<String,Boolean> germlineGenesReporting = PatientReporterTestUtil.testGermlineModel().germlineGenes();
        Map<String, DriverCategory> driverCategoryMapMatch = Maps.newHashMap();
        driverCategoryMapMatch.put("ATM", DriverCategory.TSG);
        String sampleId = "CPCT02990001T";
        List<GermlineVariant> germlineVariantsMatch = createTestGermlineVariantsTSGGene(true);
        List<GeneCopyNumber> geneCopyNumbersMatch = Lists.newArrayList(createTestCopyNumberBuilder(1).build());
        List<EnrichedSomaticVariant> variantsMatch = Lists.newArrayList(createEnriched("ATM").build());
        List<GermlineVariant> filteredGermlineVariantMatch = FilterGermlineVariants.filteringReportedGermlineVariant(germlineVariantsMatch,
                germlineGenesReporting,
                driverCategoryMapMatch,
                geneCopyNumbersMatch,
                sampleId,
                variantsMatch);
        assertEquals(filteredGermlineVariantMatch.size(), 1);

        List<GermlineVariant> germlineVariantsNonMatch = createTestGermlineVariantsTSGGene(true);
        List<GeneCopyNumber> geneCopyNumbersNonMatch = Lists.newArrayList(createTestCopyNumberBuilder(1).build());
        Map<String, DriverCategory> driverCategoryMapNonMatch = Maps.newHashMap();
        driverCategoryMapNonMatch.put("BRCA1", DriverCategory.TSG);
        List<EnrichedSomaticVariant> variantsNonMatch = Lists.newArrayList(createEnriched("ATM").build());
        List<GermlineVariant> filteredGermlineVariantNonMatch = FilterGermlineVariants.filteringReportedGermlineVariant(germlineVariantsNonMatch,
                germlineGenesReporting,
                driverCategoryMapNonMatch,
                geneCopyNumbersNonMatch,
                sampleId,
                variantsNonMatch);
        assertEquals(filteredGermlineVariantNonMatch.size(), 1);
    }

    @Test
    public void filteringONCOGenesForGermlineVariantsCheckONCO() throws IOException {
        List<GermlineVariant> germlineVariants = createTestGermlineVariantsONCOGene();
        Map<String, Boolean> germlineGenesReporting = PatientReporterTestUtil.testGermlineModel().germlineGenes();
        Map<String, DriverCategory> driverCategoryMap = Maps.newHashMap();
        driverCategoryMap.put("BRCA2", DriverCategory.ONCO);

        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();
        String sampleId = "CPCT02990001T";
        List<EnrichedSomaticVariant> variants = Lists.newArrayList();

        List<GermlineVariant> filteredGermlineVariantONCO = FilterGermlineVariants.filteringReportedGermlineVariant(germlineVariants,
                germlineGenesReporting,
                driverCategoryMap,
                geneCopyNumbers,
                sampleId,
                variants);
        assertEquals(filteredGermlineVariantONCO.size(), 1);

        Map<String, DriverCategory> driverCategoryMapTSG = Maps.newHashMap();
        driverCategoryMapTSG.put("ATM", DriverCategory.TSG);

        List<GermlineVariant> filteredGermlineVariantTSG = FilterGermlineVariants.filteringReportedGermlineVariant(germlineVariants,
                germlineGenesReporting,
                driverCategoryMapTSG,
                geneCopyNumbers,
                sampleId,
                variants);
        assertEquals(filteredGermlineVariantTSG.size(), 0);
    }

    @Test
    public void filterTSGGenesForGermlineVariantsCheckTSG() throws IOException {
        Map<String, Boolean> germlineGenesReporting = PatientReporterTestUtil.testGermlineModel().germlineGenes();
        Map<String, DriverCategory> driverCategoryMap = Maps.newHashMap();
        driverCategoryMap.put("ATM", DriverCategory.TSG);
        String sampleId = "CPCT02990001T";
        List<GermlineVariant> germlineVariantsAll = createTestGermlineVariantsTSGGene(true);
        List<GeneCopyNumber> geneCopyNumbersAll = Lists.newArrayList(createTestCopyNumberBuilder(1).build());
        List<EnrichedSomaticVariant> variantsAll = Lists.newArrayList(createEnriched("ATM").build());
        List<GermlineVariant> filteredGermlineVariantAll = FilterGermlineVariants.filteringReportedGermlineVariant(germlineVariantsAll,
                germlineGenesReporting,
                driverCategoryMap,
                geneCopyNumbersAll,
                sampleId,
                variantsAll);
        assertEquals(filteredGermlineVariantAll.size(), 1);

        List<GermlineVariant> germlineVariantsWrongBiallelic = createTestGermlineVariantsTSGGene(false);
        List<GeneCopyNumber> geneCopyNumbersWrongBiallelic = Lists.newArrayList(createTestCopyNumberBuilder(1).build());
        List<EnrichedSomaticVariant> variantsWrongBiallelic = Lists.newArrayList(createEnriched("ATM").build());
        List<GermlineVariant> filteredGermlineVariantWrongBiallelic = FilterGermlineVariants.filteringReportedGermlineVariant(
                germlineVariantsWrongBiallelic,
                germlineGenesReporting,
                driverCategoryMap,
                geneCopyNumbersWrongBiallelic,
                sampleId,
                variantsWrongBiallelic);
        assertEquals(filteredGermlineVariantWrongBiallelic.size(), 1);

        List<GermlineVariant> germlineVariantsWrongSomaticGene = createTestGermlineVariantsTSGGene(true);
        List<GeneCopyNumber> geneCopyNumbersWrongSomaticGene = Lists.newArrayList(createTestCopyNumberBuilder(1).build());
        List<EnrichedSomaticVariant> variantsWrongSomaticGene = Lists.newArrayList(createEnriched("ALK").build());
        List<GermlineVariant> filteredGermlineVariantWrongSomaticGene = FilterGermlineVariants.filteringReportedGermlineVariant(
                germlineVariantsWrongSomaticGene,
                germlineGenesReporting,
                driverCategoryMap,
                geneCopyNumbersWrongSomaticGene,
                sampleId,
                variantsWrongSomaticGene);
        assertEquals(filteredGermlineVariantWrongSomaticGene.size(), 1);

        List<GermlineVariant> germlineVariantsWrongGeneCopyNumber = createTestGermlineVariantsTSGGene(true);
        List<GeneCopyNumber> geneCopyNumbersWrongGeneCopyNumber = Lists.newArrayList(createTestCopyNumberBuilder(2).build());
        List<EnrichedSomaticVariant> variantsWrongGeneCopyNumber = Lists.newArrayList(createEnriched("ATM").build());
        List<GermlineVariant> filteredGermlineVariantWrongGeneCopyNumber = FilterGermlineVariants.filteringReportedGermlineVariant(
                germlineVariantsWrongGeneCopyNumber,
                germlineGenesReporting,
                driverCategoryMap,
                geneCopyNumbersWrongGeneCopyNumber,
                sampleId,
                variantsWrongGeneCopyNumber);
        assertEquals(filteredGermlineVariantWrongGeneCopyNumber.size(), 1);
    }

    @NotNull
    public static ImmutableGeneCopyNumber.Builder createTestCopyNumberBuilder(int maxNumberFilter) {
        return ImmutableGeneCopyNumber.builder()
                .start(1)
                .end(2)
                .gene("ATM")
                .chromosome("1")
                .chromosomeBand(Strings.EMPTY)
                .minRegionStart(0)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEnd(0)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegions(1)
                .germlineHet2HomRegions(0)
                .germlineHomRegions(0)
                .somaticRegions(1)
                .minCopyNumber(0)
                .maxCopyNumber(maxNumberFilter)
                .transcriptID("trans")
                .transcriptVersion(0)
                .minMinorAllelePloidy(0);
    }

    @NotNull
    public static ImmutableSomaticVariantImpl.Builder create(@NotNull String gene) {
        return ImmutableSomaticVariantImpl.builder()
                .chromosome(Strings.EMPTY)
                .position(0L)
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
                .adjustedCopyNumber(0D)
                .adjustedVAF(0D)
                .minorAllelePloidy(0D)
                .germlineStatus(GermlineStatus.UNKNOWN)
                .ploidy(0)
                .mappability(0D);
    }

    @NotNull
    public static ImmutableEnrichedSomaticVariant.Builder createEnriched(@NotNull String gene) {
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
    private static ImmutableEnrichedSomaticVariant.Builder builder() {
        return SomaticVariantTestBuilderFactory.createEnriched().filter("PASS");
    }

    @NotNull
    private static List<GermlineVariant> createTestGermlineVariantsONCOGene() {
        List<GermlineVariant> germlineVariants = Lists.newArrayList();

        int totalReads = 112;
        int altReads = 67;
        double adjustedCopyNumber = 3D;

        germlineVariants.add(ImmutableGermlineVariant.builder()
                .passFilter(true)
                .gene("BRCA2")
                .hgvsCodingImpact("c.5946delT")
                .hgvsProteinImpact("p.Ser1982fs")
                .totalReadCount(totalReads)
                .alleleReadCount(altReads)
                .germlineStatus("HET")
                .adjustedCopyNumber(adjustedCopyNumber)
                .adjustedVAF(12)
                .minorAllelePloidy(1D)
                .biallelic(false)
                .build());

        return germlineVariants;
    }

    @NotNull
    private static List<GermlineVariant> createTestGermlineVariantsTSGGene(boolean biallelicFilter) {
        List<GermlineVariant> germlineVariants = Lists.newArrayList();

        int totalReads = 112;
        int altReads = 67;
        double adjustedCopyNumber = 3D;

        germlineVariants.add(ImmutableGermlineVariant.builder()
                .passFilter(true)
                .gene("ATM")
                .hgvsCodingImpact("c.5946delT")
                .hgvsProteinImpact("p.Ser1982fs")
                .totalReadCount(totalReads)
                .alleleReadCount(altReads)
                .germlineStatus("HET")
                .adjustedCopyNumber(adjustedCopyNumber)
                .adjustedVAF(12)
                .minorAllelePloidy(1D)
                .biallelic(biallelicFilter)
                .build());

        return germlineVariants;
    }

}