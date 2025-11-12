package com.hartwig.hmftools.orange.algo.wildtype;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.wildtype.WildTypeGene;
import com.hartwig.hmftools.orange.algo.linx.LinxOrangeTestFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleGainDeletionFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class WildTypeAlgoTest
{
    @Test
    public void canDetermineWildTypeSomatic()
    {
        Map<String,DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BRCA2"));

        PurpleVariant variantSomatic =
                TestPurpleVariantFactory.builder().gene("BRCA2").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList(variantSomatic);
        List<PurpleVariant> reportableGermlineVariants = null;

        List<PurpleGainDeletion> reportableSomaticGainsDels = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsDels,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeGermline()
    {
        Map<String,DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BRCA1"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        PurpleVariant variantGermline =
                TestPurpleVariantFactory.builder().gene("BRCA1").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableGermlineVariants = Lists.newArrayList(variantGermline);

        List<PurpleGainDeletion> reportableSomaticGainsDels = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsDels,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeCNV()
    {
        Map<String,DriverGene> driverGenes = createDriverMap(Lists.newArrayList("APC", "KRAS"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        PurpleGainDeletion reportableAmp = TestPurpleGainDeletionFactory.createGainDel("KRAS", CopyNumberInterpretation.FULL_GAIN);
        PurpleGainDeletion reportableDel = TestPurpleGainDeletionFactory.createGainDel("APC", CopyNumberInterpretation.FULL_DEL);
        List<PurpleGainDeletion> reportableSomaticGainsDels = Lists.newArrayList(reportableAmp, reportableDel);
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsDels,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeFusion5prime()
    {
        Map<String,DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BAG4"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainDeletion> reportableSomaticGainsDels = Lists.newArrayList();
        LinxFusion reportedFusionMatch = createFusion("BAG4", "EGFR");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsDels,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeFusion3prime()
    {
        Map<String,DriverGene> driverGenes = createDriverMap(Lists.newArrayList("BAG4"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainDeletion> reportableSomaticGainsDels = Lists.newArrayList();
        LinxFusion reportedFusionMatch = createFusion("EGFR", "BAG4");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsDels,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeHomozygousDisruption()
    {
        Map<String,DriverGene> driverGenes = createDriverMap(Lists.newArrayList("NRAS"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainDeletion> reportableSomaticGainsDels = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        LinxHomozygousDisruption homozygousDisruption = createHomDisruption("NRAS");
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList(homozygousDisruption);
        List<LinxBreakend> reportableBreakends = Lists.newArrayList();

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsDels,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);

        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildTypeGeneDisruption()
    {
        Map<String,DriverGene> driverGenes = createDriverMap(Lists.newArrayList("MYC"));

        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList();
        List<PurpleVariant> reportableGermlineVariants = null;
        List<PurpleGainDeletion> reportableSomaticGainsDels = Lists.newArrayList();
        List<LinxFusion> reportableFusions = Lists.newArrayList();
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList();

        LinxBreakend breakend = createBreakend("MYC");
        List<LinxBreakend> reportableBreakends = Lists.newArrayList(breakend);

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsDels,
                reportableFusions,
                homozygousDisruptions,
                reportableBreakends);
        assertEquals(0, wildTypes.size());
    }

    @Test
    public void canDetermineWildType()
    {
        Map<String,DriverGene> driverGenes =
                createDriverMap(Lists.newArrayList("BRCA1", "BRCA2", "APC", "KRAS", "BAG4", "FGFR1", "NRAS", "EGFR", "MYC"));

        PurpleVariant variantSomatic =
                TestPurpleVariantFactory.builder().gene("BRCA2").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableSomaticVariants = Lists.newArrayList(variantSomatic);

        PurpleVariant variantGermline =
                TestPurpleVariantFactory.builder().gene("BRCA1").chromosome("1").position(56412).ref("A").alt("C").build();
        List<PurpleVariant> reportableGermlineVariants = Lists.newArrayList(variantGermline);

        PurpleGainDeletion reportableAmp = TestPurpleGainDeletionFactory.createGainDel("APC", CopyNumberInterpretation.FULL_GAIN);
        PurpleGainDeletion reportableDel = TestPurpleGainDeletionFactory.createGainDel("KRAS", CopyNumberInterpretation.FULL_DEL);
        List<PurpleGainDeletion> reportableSomaticGainsDels = Lists.newArrayList(reportableAmp, reportableDel);

        LinxFusion reportedFusionMatch = createFusion("BAG4", "FGFR1");
        List<LinxFusion> reportableFusions = Lists.newArrayList(reportedFusionMatch);

        LinxHomozygousDisruption homozygousDisruption = createHomDisruption("NRAS");
        List<LinxHomozygousDisruption> homozygousDisruptions = Lists.newArrayList(homozygousDisruption);

        LinxBreakend breakend = createBreakend("MYC");
        List<LinxBreakend> reportableBreakends = Lists.newArrayList(breakend);

        List<WildTypeGene> wildTypes = WildTypeAlgo.determineWildTypeGenes(driverGenes,
                reportableSomaticVariants,
                reportableGermlineVariants,
                reportableSomaticGainsDels,
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

    private static Map<String,DriverGene> createDriverMap(final List<String> genes)
    {
        Map<String,DriverGene> driverGeneList = Maps.newHashMap();

        for(String gene : genes)
        {
            driverGeneList.put(gene, DriverGeneTestFactory.builder().gene(gene).build());
        }
        return driverGeneList;
    }

    private static LinxHomozygousDisruption createHomDisruption(final String gene)
    {
        return LinxOrangeTestFactory.homozygousDisruptionBuilder().gene(gene).build();
    }

    public static LinxBreakend createBreakend(final String gene)
    {
        return LinxOrangeTestFactory.breakendBuilder().gene(gene).isCanonical(true).build();
    }

    private static LinxFusion createFusion(final String geneStart, final String geneEnd)
    {
        return LinxOrangeTestFactory.fusionBuilder().geneStart(geneStart).geneEnd(geneEnd).reported(true).build();
    }
}