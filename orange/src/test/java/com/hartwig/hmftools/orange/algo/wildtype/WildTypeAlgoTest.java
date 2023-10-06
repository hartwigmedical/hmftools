package com.hartwig.hmftools.orange.algo.wildtype;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.wildtype.WildTypeGene;
import com.hartwig.hmftools.orange.algo.linx.LinxOrangeTestFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleGainLossFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class WildTypeAlgoTest
{
    @Test
    public void canDetermineWildTypeSomatic()
    {
        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BRCA2"));

        PurpleVariant variantSomatic =
                TestPurpleVariantFactory.builder().gene("BRCA2").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList(variantSomatic);
        List<PurpleVariant> reportableGermlineVariants = null;

        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeGermline()
    {
        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BRCA1"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        PurpleVariant variantGermline =
                TestPurpleVariantFactory.builder().gene("BRCA1").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableGermlineVariants = Lists.newArrayList(variantGermline);

        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeCNV()
    {
        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("APC", "KRAS"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        PurpleGainLoss reportableAmp = TestPurpleGainLossFactory.createGainLoss("KRAS", CopyNumberInterpretation.FULL_GAIN);
        PurpleGainLoss reportableDel = TestPurpleGainLossFactory.createGainLoss("APC", CopyNumberInterpretation.FULL_LOSS);
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList(reportableAmp, reportableDel);
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeFusion5prime()
    {
        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BAG4"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        LinxFusion reportedFusionMatch = createFusion("BAG4", "EGFR");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeFusion3prime()
    {
        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BAG4"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        LinxFusion reportedFusionMatch = createFusion("EGFR", "BAG4");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeHomozygousDisruption()
    {
        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("NRAS"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        LinxHomozygousDisruption homozygousDisruption = createHomDisruption("NRAS");
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList(homozygousDisruption);
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);

        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeGeneDisruption()
    {
        List<DriverGene> driverGenes = createDriverMap(Lists.newArrayList("MYC"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();

        LinxBreakend breakend = createBreakend("MYC");
        List<LinxBreakend> reportableBreakends = Lists.newArrayList(breakend);

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildType()
    {
        List<DriverGene> driverGenes =
                createDriverMap(Lists.newArrayList("BRCA1", "BRCA2", "APC", "KRAS", "BAG4", "FGFR1", "NRAS", "EGFR", "MYC"));

        PurpleVariant variantSomatic =
                TestPurpleVariantFactory.builder().gene("BRCA2").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList(variantSomatic);

        PurpleVariant variantGermline =
                TestPurpleVariantFactory.builder().gene("BRCA1").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableGermlineVariants = Lists.newArrayList(variantGermline);

        PurpleGainLoss reportableAmp = TestPurpleGainLossFactory.createGainLoss("APC", CopyNumberInterpretation.FULL_GAIN);
        PurpleGainLoss reportableDel = TestPurpleGainLossFactory.createGainLoss("KRAS", CopyNumberInterpretation.FULL_LOSS);
        List<PurpleGainLoss> reportableSomaticGainsLosses = Lists.newArrayList(reportableAmp, reportableDel);

        LinxFusion reportedFusionMatch = createFusion("BAG4", "FGFR1");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);

        LinxHomozygousDisruption homozygousDisruption = createHomDisruption("NRAS");
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList(homozygousDisruption);

        LinxBreakend breakend = createBreakend("MYC");
        List<LinxBreakend> reportableBreakends = Lists.newArrayList(breakend);

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsLosses,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(1, wildTypes.size());
    }

    @Test
    public void canAssessQualityForWildTypeCalling()
    {
        Set<PurpleQCStatus> purpleQCStatusSetPASS = Sets.newHashSet();
        purpleQCStatusSetPASS.add(PurpleQCStatus.PASS);
        assertTrue(WildTypeAlgo.wildTypeCallingAllowed(purpleQCStatusSetPASS));

        Set<PurpleQCStatus> purpleQCStatusSetWarnDeleted = Sets.newHashSet();
        purpleQCStatusSetWarnDeleted.add(PurpleQCStatus.WARN_DELETED_GENES);
        assertTrue(WildTypeAlgo.wildTypeCallingAllowed(purpleQCStatusSetWarnDeleted));

        Set<PurpleQCStatus> purpleQCStatusSetFailPurity = Sets.newHashSet();
        purpleQCStatusSetFailPurity.add(PurpleQCStatus.FAIL_NO_TUMOR);
        assertFalse(WildTypeAlgo.wildTypeCallingAllowed(purpleQCStatusSetFailPurity));

        Set<PurpleQCStatus> purpleQCStatusSetWarnPurity = Sets.newHashSet();
        purpleQCStatusSetWarnPurity.add(PurpleQCStatus.WARN_LOW_PURITY);
        assertFalse(WildTypeAlgo.wildTypeCallingAllowed(purpleQCStatusSetWarnPurity));
    }

    @NotNull
    private static List<DriverGene> createDriverMap(@NotNull List<String> genes)
    {
        List<DriverGene> driverGeneList = Lists.newArrayList();
        for(String gene : genes)
        {
            driverGeneList.add(DriverGeneTestFactory.builder().gene(gene).build());
        }
        return driverGeneList;
    }

    @NotNull
    private static LinxHomozygousDisruption createHomDisruption(@NotNull String gene)
    {
        return LinxOrangeTestFactory.homozygousDisruptionBuilder().gene(gene).build();
    }

    @NotNull
    public static LinxBreakend createBreakend(@NotNull String gene)
    {
        return LinxOrangeTestFactory.breakendBuilder().gene(gene).isCanonical(true).build();
    }

    @NotNull
    private static LinxFusion createFusion(@NotNull String geneStart, @NotNull String geneEnd)
    {
        return LinxOrangeTestFactory.fusionBuilder().geneStart(geneStart).geneEnd(geneEnd).reported(true).build();
    }
}