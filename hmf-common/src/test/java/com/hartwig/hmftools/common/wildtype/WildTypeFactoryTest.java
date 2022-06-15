package com.hartwig.hmftools.common.wildtype;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.linx.ImmutableReportableGeneDisruption;
import com.hartwig.hmftools.common.linx.ImmutableReportableHomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.linx.ReportableGeneDisruption;
import com.hartwig.hmftools.common.linx.ReportableGeneDisruptionFactory;
import com.hartwig.hmftools.common.linx.ReportableGeneDisruptionFactoryTest;
import com.hartwig.hmftools.common.linx.ReportableHomozygousDisruption;
import com.hartwig.hmftools.common.purple.interpretation.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.purple.interpretation.GainLossTestFactory;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.ImmutableReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantTestFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class WildTypeFactoryTest {

    @Test
    public void canDetermineWildTypeSomatic() {
        List<ReportableVariant> reportableGermlineVariant = Lists.newArrayList();

        ReportableVariant variantSomatic = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("BRCA2")
                .chromosome("1")
                .position(56412)
                .ref("A")
                .alt("C")
                .build();
        List<ReportableVariant> reportableSomaticVariant = Lists.newArrayList(variantSomatic);

        List<GainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> geneDisruptions = Lists.newArrayList();
        List<DriverGene> listDriverGenes =
                createDriverMap(Lists.newArrayList("BRCA2"));

        List<WildType> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions,
                listDriverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeGermline() {
        ReportableVariant variantGermline = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("BRCA1")
                .chromosome("1")
                .position(56412)
                .ref("A")
                .alt("C")
                .build();
        List<ReportableVariant> reportableGermlineVariant = Lists.newArrayList(variantGermline);

        List<ReportableVariant> reportableSomaticVariant = Lists.newArrayList();
        List<GainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> geneDisruptions = Lists.newArrayList();

        List<DriverGene> listDriverGenes =
                createDriverMap(Lists.newArrayList("BRCA1"));

        List<WildType> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions,
                listDriverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeCNV() {
        List<ReportableVariant> reportableGermlineVariant = Lists.newArrayList();
        List<ReportableVariant> reportableSomaticVariant = Lists.newArrayList();
        GainLoss reportableAmp = GainLossTestFactory.createGainLoss("KRAS", CopyNumberInterpretation.FULL_GAIN);
        GainLoss reportableDel = GainLossTestFactory.createGainLoss("APC", CopyNumberInterpretation.FULL_LOSS);
        List<GainLoss> reportableSomaticGainsLosses = Lists.newArrayList(reportableAmp, reportableDel);
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> geneDisruptions = Lists.newArrayList();

        List<DriverGene> listDriverGenes =
                createDriverMap(Lists.newArrayList("APC", "KRAS"));

        List<WildType> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions,
                listDriverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeFusion5prime() {
        List<ReportableVariant> reportableGermlineVariant = Lists.newArrayList();
        List<ReportableVariant> reportableSomaticVariant = Lists.newArrayList();
        List<GainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        LinxFusion reportedFusionMatch = create("BAG4", "EGFR");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> geneDisruptions = Lists.newArrayList();

        List<DriverGene> listDriverGenes =
                createDriverMap(Lists.newArrayList("BAG4"));

        List<WildType> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions,
                listDriverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeFusion3prime() {
        List<ReportableVariant> reportableGermlineVariant = Lists.newArrayList();
        List<ReportableVariant> reportableSomaticVariant = Lists.newArrayList();
        List<GainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        LinxFusion reportedFusionMatch = create("EGFR", "BAG4");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<ReportableGeneDisruption> geneDisruptions = Lists.newArrayList();

        List<DriverGene> listDriverGenes =
                createDriverMap(Lists.newArrayList("BAG4"));

        List<WildType> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions,
                listDriverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeHomozygousDisruption() {
        List<ReportableVariant> reportableGermlineVariant = Lists.newArrayList();
        List<ReportableVariant> reportableSomaticVariant = Lists.newArrayList();
        List<GainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        ReportableHomozygousDisruption homozygousDisruption = create("NRAS");
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList(homozygousDisruption);
        List<ReportableGeneDisruption> geneDisruptions = Lists.newArrayList();

        List<DriverGene> listDriverGenes =
                createDriverMap(Lists.newArrayList("NRAS"));

        List<WildType> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions,
                listDriverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeGeneDisruption() {
        List<ReportableVariant> reportableGermlineVariant = Lists.newArrayList();
        List<ReportableVariant> reportableSomaticVariant = Lists.newArrayList();
        List<GainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();

        ReportableGeneDisruption geneDisruption = createDisruption("MYC");
        List<ReportableGeneDisruption> geneDisruptions = Lists.newArrayList(geneDisruption);

        List<DriverGene> listDriverGenes =
                createDriverMap(Lists.newArrayList("MYC"));

        List<WildType> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions,
                listDriverGenes);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildType() {

        ReportableVariant variantGermline = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("BRCA1")
                .chromosome("1")
                .position(56412)
                .ref("A")
                .alt("C")
                .build();
        List<ReportableVariant> reportableGermlineVariant = Lists.newArrayList(variantGermline);

        ReportableVariant variantSomatic = ImmutableReportableVariant.builder()
                .from(ReportableVariantTestFactory.create())
                .gene("BRCA2")
                .chromosome("1")
                .position(56412)
                .ref("A")
                .alt("C")
                .build();
        List<ReportableVariant> reportableSomaticVariant = Lists.newArrayList(variantSomatic);

        GainLoss reportableAmp = GainLossTestFactory.createGainLoss("APC", CopyNumberInterpretation.FULL_GAIN);
        GainLoss reportableDel = GainLossTestFactory.createGainLoss("KRAS", CopyNumberInterpretation.FULL_LOSS);
        List<GainLoss> reportableSomaticGainsLosses = Lists.newArrayList(reportableAmp, reportableDel);

        LinxFusion reportedFusionMatch = create("BAG4", "FGFR1");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);

        ReportableHomozygousDisruption homozygousDisruption = create("NRAS");
        List<ReportableHomozygousDisruption> homozygousDisruptions = Lists.newArrayList(homozygousDisruption);

        ReportableGeneDisruption geneDisruption = createDisruption("MYC");
        List<ReportableGeneDisruption> geneDisruptions = Lists.newArrayList(geneDisruption);

        List<DriverGene> listDriverGenes =
                createDriverMap(Lists.newArrayList("BRCA1", "BRCA2", "APC", "KRAS", "BAG4", "FGFR1", "NRAS", "EGFR", "MYC"));

        List<WildType> wildTypes = WildTypeFactory.determineWildTypeGenes(reportableGermlineVariant,
                reportableSomaticVariant,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions, geneDisruptions,
                listDriverGenes);
        assertEquals(1, wildTypes.size());
    }

    @NotNull
    private static List<DriverGene> createDriverMap(@NotNull List<String> genes) {
        List<DriverGene> driverGeneList = Lists.newArrayList();
        for (String gene : genes) {
            driverGeneList.add(createDriverGene(gene));
        }
        return driverGeneList;
    }

    public static DriverGene createDriverGene(final String name) {
        return ImmutableDriverGene.builder()
                .gene(name)
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(true)
                .reportAmplification(false)
                .reportSomaticHotspot(false)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .likelihoodType(TSG)
                .reportGermlineDisruption(true)
                .build();
    }

    @NotNull
    private static ReportableHomozygousDisruption create(@NotNull String gene) {
        return ImmutableReportableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript("123")
                .isCanonical(true)
                .build();
    }

    @NotNull
    public static ReportableGeneDisruption createDisruption(@NotNull String gene) {
        return createTestReportableGeneDisruptionBuilder().gene(gene).isCanonical(true).build();
    }

    @Override
    public boolean equals(final Object obj) {
        return super.equals(obj);
    }

    private static ImmutableReportableGeneDisruption.Builder createTestReportableGeneDisruptionBuilder() {
        return ImmutableReportableGeneDisruption.builder()
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
        return linxFusionBuilder(geneStart, geneEnd).build();
    }

    @NotNull
    private static ImmutableLinxFusion.Builder linxFusionBuilder(@NotNull String geneStart, @NotNull String geneEnd) {
        return ImmutableLinxFusion.builder()
                .from(LinxTestFactory.createMinimalTestFusion())
                .geneStart(geneStart)
                .geneEnd(geneEnd)
                .reported(true);
    }
}