package com.hartwig.hmftools.orange.algo.wildtype;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.ImmutableHomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.orange.algo.linx.GeneDisruption;
import com.hartwig.hmftools.orange.algo.linx.ImmutableGeneDisruption;
import com.hartwig.hmftools.orange.algo.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.orange.algo.purple.GainLossTestFactory;
import com.hartwig.hmftools.orange.algo.purple.PurpleGainLoss;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariantTestFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class WildTypeFactoryTest {

    @Test
    public void canDetermineWildTypeSomatic() {
        PurpleVariant variantSomatic =
                PurpleVariantTestFactory.builder().gene("BRCA2").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList(variantSomatic);
        List<PurpleVariant> reportableGermlineVariants = null;

        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<GeneDisruption> reportableGeneDisruptions = Lists.newArrayList();
        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BRCA2"));

        List<WildTypeGene> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableGeneDisruptions,
                driverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeGermline() {
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        PurpleVariant variantGermline =
                PurpleVariantTestFactory.builder().gene("BRCA1").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableGermlineVariants = Lists.newArrayList(variantGermline);

        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<GeneDisruption> geneDisruptions = Lists.newArrayList();

        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BRCA1"));

        List<WildTypeGene> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                geneDisruptions,
                driverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeCNV() {
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        PurpleGainLoss reportableAmp = GainLossTestFactory.createGainLoss("KRAS", CopyNumberInterpretation.FULL_GAIN);
        PurpleGainLoss reportableDel = GainLossTestFactory.createGainLoss("APC", CopyNumberInterpretation.FULL_LOSS);
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList(reportableAmp, reportableDel);
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<GeneDisruption> reportableGeneDisruptions = Lists.newArrayList();

        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("APC", "KRAS"));

        List<WildTypeGene> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableGeneDisruptions,
                driverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeFusion5prime() {
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        LinxFusion reportedFusionMatch = create("BAG4", "EGFR");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<GeneDisruption> reportableGeneDisruptions = Lists.newArrayList();

        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BAG4"));

        List<WildTypeGene> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableGeneDisruptions,
                driverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeFusion3prime() {
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        LinxFusion reportedFusionMatch = create("EGFR", "BAG4");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<GeneDisruption> reportableGeneDisruptions = Lists.newArrayList();

        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BAG4"));

        List<WildTypeGene> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableGeneDisruptions,
                driverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeHomozygousDisruption() {
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        HomozygousDisruption homozygousDisruption = create("NRAS");
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList(homozygousDisruption);
        List<GeneDisruption> reportableGeneDisruptions = Lists.newArrayList();

        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("NRAS"));

        List<WildTypeGene> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableGeneDisruptions,
                driverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeGeneDisruption() {
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList();

        GeneDisruption geneDisruption = createDisruption("MYC");
        List<GeneDisruption> geneDisruptions = Lists.newArrayList(geneDisruption);

        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("MYC"));

        List<WildTypeGene> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                geneDisruptions,
                driverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildType() {
        PurpleVariant variantSomatic =
                PurpleVariantTestFactory.builder().gene("BRCA2").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList(variantSomatic);

        PurpleVariant variantGermline =
                PurpleVariantTestFactory.builder().gene("BRCA1").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableGermlineVariants = Lists.newArrayList(variantGermline);

        PurpleGainLoss reportableAmp = GainLossTestFactory.createGainLoss("APC", CopyNumberInterpretation.FULL_GAIN);
        PurpleGainLoss reportableDel = GainLossTestFactory.createGainLoss("KRAS", CopyNumberInterpretation.FULL_LOSS);
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList(reportableAmp, reportableDel);

        LinxFusion reportedFusionMatch = create("BAG4", "FGFR1");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);

        HomozygousDisruption homozygousDisruption = create("NRAS");
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList(homozygousDisruption);

        GeneDisruption geneDisruption = createDisruption("MYC");
        List<GeneDisruption> reportableGeneDisruptions = Lists.newArrayList(geneDisruption);

        List<DriverGene> driverGenes =
                createDriverMap(Lists.newArrayList("BRCA1", "BRCA2", "APC", "KRAS", "BAG4", "FGFR1", "NRAS", "EGFR", "MYC"));

        List<WildTypeGene> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableGeneDisruptions,
                driverGenes);
        assertEquals(1, wildTypes.size());
    }

    @Test
    public void canFilterWildType() {
        PurpleVariant variantSomatic =
                PurpleVariantTestFactory.builder().gene("BRCA2").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList(variantSomatic);

        PurpleVariant variantGermline =
                PurpleVariantTestFactory.builder().gene("BRCA1").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableGermlineVariants = Lists.newArrayList(variantGermline);

        PurpleGainLoss reportableAmp = GainLossTestFactory.createGainLoss("APC", CopyNumberInterpretation.FULL_GAIN);
        PurpleGainLoss reportableDel = GainLossTestFactory.createGainLoss("KRAS", CopyNumberInterpretation.FULL_LOSS);
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList(reportableAmp, reportableDel);

        LinxFusion reportedFusionMatch = create("BAG4", "FGFR1");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);

        HomozygousDisruption homozygousDisruption = create("NRAS");
        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList(homozygousDisruption);

        GeneDisruption geneDisruption = createDisruption("MYC");
        List<GeneDisruption> reportableGeneDisruptions = Lists.newArrayList(geneDisruption);

        List<DriverGene> driverGenes =
                createDriverMap(Lists.newArrayList("BRCA1", "BRCA2", "APC", "KRAS", "BAG4", "FGFR1", "NRAS", "EGFR", "MYC"));

        List<WildTypeGene> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableGeneDisruptions,
                driverGenes);

        Set<PurpleQCStatus> purpleQCStatusSetPASS = Sets.newHashSet();
        purpleQCStatusSetPASS.add(PurpleQCStatus.PASS);
        assertEquals(1, WildTypeFactory.filterQCWildTypes(purpleQCStatusSetPASS, wildTypes).size());

        Set<PurpleQCStatus> purpleQCStatusSetWarnDeleted = Sets.newHashSet();
        purpleQCStatusSetWarnDeleted.add(PurpleQCStatus.WARN_DELETED_GENES);
        assertEquals(1, WildTypeFactory.filterQCWildTypes(purpleQCStatusSetWarnDeleted, wildTypes).size());

        Set<PurpleQCStatus> purpleQCStatusSetFailPurity = Sets.newHashSet();
        purpleQCStatusSetFailPurity.add(PurpleQCStatus.FAIL_NO_TUMOR);
        assertEquals(0, WildTypeFactory.filterQCWildTypes(purpleQCStatusSetFailPurity, wildTypes).size());

        Set<PurpleQCStatus> purpleQCStatusSetWarnPurity = Sets.newHashSet();
        purpleQCStatusSetWarnPurity.add(PurpleQCStatus.WARN_LOW_PURITY);
        assertEquals(0, WildTypeFactory.filterQCWildTypes(purpleQCStatusSetWarnPurity, wildTypes).size());
    }

    @NotNull
    private static List<DriverGene> createDriverMap(@NotNull List<String> genes) {
        List<DriverGene> driverGeneList = Lists.newArrayList();
        for (String gene : genes) {
            driverGeneList.add(DriverGeneTestFactory.builder().gene(gene).build());
        }
        return driverGeneList;
    }

    @NotNull
    private static HomozygousDisruption create(@NotNull String gene) {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript("123")
                .isCanonical(true)
                .build();
    }

    @NotNull
    public static GeneDisruption createDisruption(@NotNull String gene) {
        return createTestReportableGeneDisruptionBuilder().gene(gene).isCanonical(true).build();
    }

    @NotNull
    private static ImmutableGeneDisruption.Builder createTestReportableGeneDisruptionBuilder() {
        return ImmutableGeneDisruption.builder()
                .location(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .range(Strings.EMPTY)
                .type(Strings.EMPTY)
                .junctionCopyNumber(2.012)
                .undisruptedCopyNumber(0.0)
                .firstAffectedExon(5)
                .clusterId(2)
                .transcriptId(Strings.EMPTY);
    }

    @NotNull
    private static LinxFusion create(@NotNull String geneStart, @NotNull String geneEnd) {
        return LinxTestFactory.fusionBuilder().geneStart(geneStart).geneEnd(geneEnd).reported(true).build();
    }
}