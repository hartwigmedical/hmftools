package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.ANY;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_HET_DELETION_THRESHOLD;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_2;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.common.variant.HotspotType.HOTSPOT_FLAG;
import static com.hartwig.hmftools.common.variant.HotspotType.NEAR_HOTSPOT_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser.VAR_IMPACT;
import static com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser.toVcfData;
import static com.hartwig.hmftools.purple.drivers.OncoDriversTest.createGeneCopyNumber;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.driver.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.test.VariantContextFromString;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineReportedTest
{
    private final GermlineDrivers mGermlineDrivers = new GermlineDrivers(Maps.newHashMap());

    private final static Map<String,GeneCopyNumber> mGeneCopyNumberMap = Maps.newHashMap();

    private final static String CLIN_SIG = "UNKNOWN";

    static
    {
        GeneCopyNumber mGeneCopyNumber = createGeneCopyNumber(GENE_NAME_1);
        GeneCopyNumber mGeneCopyNumber2 = createGeneCopyNumber(GENE_NAME_2);
        mGeneCopyNumberMap.put(mGeneCopyNumber.GeneName, mGeneCopyNumber);
        mGeneCopyNumberMap.put(mGeneCopyNumber2.GeneName, mGeneCopyNumber2);
    }

    private static Map<String,DriverGene> buildDriverGeneMap(final DriverGene driverGene)
    {
        Map<String,DriverGene> map = Maps.newHashMap();
        map.put(driverGene.gene(), driverGene);
        return map;
    }

    private void setDrivers(final DriverGene driverGene)
    {
        mGermlineDrivers.driverGeneMap().clear();
        mGermlineDrivers.driverGeneMap().put(driverGene.gene(), driverGene);
    }

    @Test
    public void testReportHotspot()
    {
        setDrivers(createDriverGene(GENE_NAME_1, DriverGeneGermlineReporting.NONE, ANY));

        GermlineVariant var = createGermlineVariant(
                PASS_FILTER, GENE_NAME_1, false, true, CLIN_SIG, CodingEffect.NONE, 0.5);

        var.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);

        List<DriverCatalog> driverCatalog = mGermlineDrivers.findDrivers(List.of(var), mGeneCopyNumberMap, Collections.emptySet());
        assertEquals(1, driverCatalog.size());
        assertEquals(ReportedStatus.REPORTED, driverCatalog.get(0).reportedStatus());
    }

    @Test
    public void testReportHotspotWhenMultipleHits()
    {
        DriverGene driverGene = createDriverGene(GENE_NAME_1, DriverGeneGermlineReporting.WILDTYPE_LOST, DriverGeneGermlineReporting.WILDTYPE_LOST);
        setDrivers(driverGene);

        GermlineVariant var1 = createGermlineVariant(PASS_FILTER, GENE_NAME_1, false, true, CLIN_SIG, CodingEffect.NONE, 0.5);
        GermlineVariant var2 = createGermlineVariant(PASS_FILTER, GENE_NAME_1, false, false, "Pathogenic", CodingEffect.NONE, 0.5);

        var1.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);
        var2.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);

        List<DriverCatalog> driverCatalog = mGermlineDrivers.findDrivers(List.of(var1, var2), mGeneCopyNumberMap, Collections.emptySet());
        assertEquals(1, driverCatalog.size());
        assertEquals(ReportedStatus.REPORTED, driverCatalog.get(0).reportedStatus());
    }

    @Test
    public void testReportHotspotWhenMultipleHitsWithSameLPS()
    {
        DriverGene driverGene = createDriverGene(GENE_NAME_1, DriverGeneGermlineReporting.WILDTYPE_LOST, DriverGeneGermlineReporting.WILDTYPE_LOST);
        setDrivers(driverGene);

        GermlineVariant var1 = createGermline(
                PASS_FILTER, GENE_NAME_1, false, true, CLIN_SIG, CodingEffect.NONE, 0.5, 4);

        GermlineVariant var2 = createGermline(
                PASS_FILTER, GENE_NAME_1, false, false, "Pathogenic", CodingEffect.NONE, 0.5, 4) ;

        List<DriverCatalog> driverCatalog = mGermlineDrivers.findDrivers(List.of(var1, var2), mGeneCopyNumberMap, Collections.emptySet());
        assertEquals(1, driverCatalog.size());
        assertEquals(ReportedStatus.NOT_REPORTED, driverCatalog.get(0).reportedStatus());
    }

    @Test
    public void testReportUnknownFrameshift()
    {
        DriverGene driverGene = createDriverGene(GENE_NAME_1, ANY, DriverGeneGermlineReporting.NONE);

        GermlineVariant var = createGermlineVariant(
                PASS_FILTER, GENE_NAME_1, false, false, CLIN_SIG, CodingEffect.NONSENSE_OR_FRAMESHIFT, 0.5);

        GermlineReportedEnrichment germlineReportedEnrichment = new GermlineReportedEnrichment(buildDriverGeneMap(driverGene));

        germlineReportedEnrichment.processVariant(var);
        germlineReportedEnrichment.flush();

        assertTrue(var.reported());
    }

    @Test
    public void testReportPathogenic()
    {
        DriverGene driverGene = createDriverGene(GENE_NAME_1, ANY, DriverGeneGermlineReporting.NONE);
        setDrivers(driverGene);

        GermlineVariant var = createGermlineVariant(
                PASS_FILTER, GENE_NAME_1, false, false, "Pathogenic", CodingEffect.NONE, 0.5);

        GermlineReportedEnrichment germlineReportedEnrichment = new GermlineReportedEnrichment(buildDriverGeneMap(driverGene));

        germlineReportedEnrichment.processVariant(var);
        germlineReportedEnrichment.flush();

        assertTrue(var.reported());
    }

    @Test
    public void testReportBiallelicAsMultipleHit()
    {
        DriverGene driverGene = createDriverGene(GENE_NAME_1, ANY, DriverGeneGermlineReporting.WILDTYPE_LOST);

        GermlineVariant var = createGermlineVariant(
                PASS_FILTER, GENE_NAME_1, true, false, "Pathogenic", CodingEffect.NONE, 0.5);

        GermlineReportedEnrichment germlineReportedEnrichment = new GermlineReportedEnrichment(buildDriverGeneMap(driverGene));

        germlineReportedEnrichment.processVariant(var);
        germlineReportedEnrichment.flush();
    }

    @Test
    public void testIgnoreSinglePathogenicWhenNoMultipleHits()
    {
        DriverGene driverGene = createDriverGene(GENE_NAME_1, DriverGeneGermlineReporting.WILDTYPE_LOST, DriverGeneGermlineReporting.NONE);

        GermlineVariant var = createGermlineVariant(
                PASS_FILTER, GENE_NAME_1, false, false, "Pathogenic", CodingEffect.NONE, 0.5);

        GermlineReportedEnrichment germlineReportedEnrichment = new GermlineReportedEnrichment(buildDriverGeneMap(driverGene));

        germlineReportedEnrichment.processVariant(var);
        germlineReportedEnrichment.flush();
    }

    @Test
    public void testSomaticCountTowardsMultipleHits()
    {
        DriverGene driverGene = createDriverGene(GENE_NAME_1, DriverGeneGermlineReporting.WILDTYPE_LOST, DriverGeneGermlineReporting.NONE);
        setDrivers(driverGene);

        GermlineVariant var = createGermlineVariant(
                PASS_FILTER, GENE_NAME_1, false, false, "Pathogenic", CodingEffect.NONE, 0.5);

        List<DriverCatalog> driverCatalog = mGermlineDrivers.findDrivers(List.of(var), mGeneCopyNumberMap, Sets.newHashSet(GENE_NAME_1));
        assertEquals(1, driverCatalog.size());
        assertEquals(ReportedStatus.REPORTED, driverCatalog.get(0).reportedStatus());
    }

    @Test
    public void testVariantNotLost()
    {
        DriverGene driverGene = createDriverGene(GENE_NAME_1, DriverGeneGermlineReporting.VARIANT_NOT_LOST, DriverGeneGermlineReporting.NONE);
        DriverGene driverGene2 = createDriverGene(GENE_NAME_2, DriverGeneGermlineReporting.VARIANT_NOT_LOST, DriverGeneGermlineReporting.NONE);

        mGermlineDrivers.driverGeneMap().clear();
        mGermlineDrivers.driverGeneMap().put(driverGene.gene(), driverGene);
        mGermlineDrivers.driverGeneMap().put(driverGene2.gene(), driverGene2);

        GermlineVariant var1 = createGermlineVariant(
                PASS_FILTER, GENE_NAME_1, false, false, "Pathogenic", CodingEffect.MISSENSE, 0.4);

        GermlineVariant var2 = createGermlineVariant(
                PASS_FILTER, GENE_NAME_2, false, false, "Pathogenic", CodingEffect.MISSENSE, 0.6);

        var1.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);
        var2.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);

        List<DriverCatalog> driverCatalogs = mGermlineDrivers.findDrivers(List.of(var1, var2), mGeneCopyNumberMap, Collections.emptySet());
        assertEquals(2, driverCatalogs.size());

        DriverCatalog driverCatalog = driverCatalogs.stream().filter(x -> x.gene().equals(GENE_NAME_1)).findFirst().orElse(null);

        assertEquals(ReportedStatus.NOT_REPORTED, driverCatalog.reportedStatus());

        driverCatalog = driverCatalogs.stream().filter(x -> x.gene().equals(GENE_NAME_2)).findFirst().orElse(null);
        assertEquals(ReportedStatus.REPORTED, driverCatalog.reportedStatus());
    }

    private static GermlineVariant createGermlineVariant(
            final String filter, String gene, boolean biallelic, boolean isHotspot,
            final String clinSig, final CodingEffect codingEffect, double variantCopyNumber)
    {
        final String hotspotFlag = isHotspot ? HOTSPOT_FLAG : NEAR_HOTSPOT_FLAG;

        VariantImpact impact = new VariantImpact(
                gene, "ENST00000393562", "UTR_variant", codingEffect,
                "c.-275T>G", "", false, "", codingEffect, 1);

        final String line2 =
                "11\t1000\tCOSM123;COSM456\tG\tA\t100\t" + filter + "\t" + PURPLE_BIALLELIC_FLAG + "=" + biallelic + ";" + hotspotFlag
                        + ";" + VAR_IMPACT + "=" + impactToString(impact) + ";CLNSIG=" + clinSig + ";"
                        + PURPLE_VARIANT_CN + "=" + variantCopyNumber + ";\"\tGT:AD:DP\t0/1:73,17:91";

        GermlineVariant variant = new GermlineVariant(VariantContextFromString.decode(line2));
        return variant;
    }

    private static GermlineVariant createGermline(
            final String filter, String gene, boolean biallelic, boolean isHotspot,
            final String clinSig, CodingEffect codingEffect, double variantCopyNumber, int localPhaseSet)
    {
        final String hotspotFlag = isHotspot ? HOTSPOT_FLAG : NEAR_HOTSPOT_FLAG;
        VariantImpact impact = new VariantImpact(
                gene, "ENST00000393562", "UTR_variant", codingEffect,
                "c.-275T>G", "", false, "", codingEffect, 1);

        String line2 =
                "11\t1000\tCOSM123;COSM456\tG\tA\t100\t" + filter + "\t" + PURPLE_BIALLELIC_FLAG + "=" + biallelic + ";" + hotspotFlag
                        + ";" + VAR_IMPACT + "=" + impactToString(impact) + ";CLNSIG=" + clinSig + ";"
                        + PURPLE_VARIANT_CN + "=" + variantCopyNumber + ";" + LOCAL_PHASE_SET + "=" + localPhaseSet
                        + "\tGT:AD:DP\t0/1:73,17:91";

        GermlineVariant variant = new GermlineVariant(VariantContextFromString.decode(line2));
        variant.context().getCommonInfo().putAttribute(REPORTED_FLAG, true);
        return variant;
    }

    private static String impactToString(final VariantImpact impact)
    {
        StringJoiner sj = new StringJoiner(",");
        toVcfData(impact).forEach(x -> sj.add(x));
        return sj.toString();
    }

    @NotNull
    private static DriverGene createDriverGene(
            final String gene, final DriverGeneGermlineReporting variantReporting, final DriverGeneGermlineReporting hotspotReporting)
    {
        return ImmutableDriverGene.builder()
                .gene(gene)
                .reportDisruption(false)
                .reportDeletion(false)
                .reportNonsenseAndFrameshift(false)
                .reportMissenseAndInframe(false)
                .reportGermlineHotspot(hotspotReporting)
                .reportGermlineVariant(variantReporting)
                .likelihoodType(DriverCategory.TSG)
                .reportAmplification(false)
                .amplificationRatio(DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO)
                .reportHetDeletion(true)
                .hetDeletionThreshold(DEFAULT_DRIVER_HET_DELETION_THRESHOLD)
                .reportSomaticHotspot(false)
                .reportSplice(false)
                .reportGermlineDisruption(ANY)
                .reportGermlineDeletion(ANY)
                .reportPGX(false)
                .build();
    }
}
