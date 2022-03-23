package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_VARIANT_CN_INFO;
import static com.hartwig.hmftools.common.variant.Hotspot.HOTSPOT_FLAG;
import static com.hartwig.hmftools.common.variant.Hotspot.NEAR_HOTSPOT_FLAG;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.sage.SageMetaData;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.test.VariantContextFromString;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineReportedEnrichmentTest
{
    @Test
    public void testReportHotspot()
    {
        DriverGene driverGene = createDriverGene("KD53", DriverGeneGermlineReporting.NONE, DriverGeneGermlineReporting.ANY);

        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Collections.emptySet());

        VariantContext var = createGermline("PASS", "KD53", false, true, "UNKNOWN", CodingEffect.NONE, 0.5);
        victim.processVariant(var);
        victim.flush();

        assertTrue(new VariantContextDecorator(var).reported());
    }

    @Test
    public void testDoNotReportFailedHotspot()
    {
        DriverGene driverGene = createDriverGene("KD53", DriverGeneGermlineReporting.NONE, DriverGeneGermlineReporting.ANY);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Collections.emptySet());

        VariantContext var1 = createGermline("PASS", "KD53", false, true, "UNKNOWN", CodingEffect.NONE, 0.5);
        VariantContext var2 = createGermline("FILTERED", "KD53", false, true, "UNKNOWN", CodingEffect.NONE, 0.5);
        victim.processVariant(var1);
        victim.processVariant(var2);
        victim.flush();

        assertTrue(new VariantContextDecorator(var1).reported());
        assertFalse(new VariantContextDecorator(var2).reported());
    }

    @Test
    public void testDoNotReportHotspotWhenNone()
    {
        DriverGene driverGene = createDriverGene("KD53", DriverGeneGermlineReporting.NONE, DriverGeneGermlineReporting.NONE);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Collections.emptySet());

        VariantContext var = createGermline("PASS", "KD53", false, true, "UNKNOWN", CodingEffect.NONE, 0.5);
        victim.processVariant(var);
        victim.flush();

        assertFalse(new VariantContextDecorator(var).reported());
    }

    @Test
    public void testDoNotReportSingleHotspotWhenNoMultipleHits()
    {
        DriverGene driverGene =
                createDriverGene("KD53", DriverGeneGermlineReporting.WILDTYPE_LOST, DriverGeneGermlineReporting.WILDTYPE_LOST);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Collections.emptySet());

        VariantContext var1 = createGermline("PASS", "KD53", false, true, "UNKNOWN", CodingEffect.NONE, 0.5);
        VariantContext var2 = createGermline("FILTERED", "KD53", false, false, "Pathogenic", CodingEffect.NONE, 0.5);
        victim.processVariant(var1);
        victim.processVariant(var2);
        victim.flush();

        assertFalse(new VariantContextDecorator(var1).reported());
        assertFalse(new VariantContextDecorator(var2).reported());
    }

    @Test
    public void testDoNotReportSingleHotspotWhenMultipleHitsViaUnreportedVariant()
    {
        DriverGene driverGene = createDriverGene("KD53", DriverGeneGermlineReporting.NONE, DriverGeneGermlineReporting.WILDTYPE_LOST);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Collections.emptySet());

        VariantContext var1 = createGermline("PASS", "KD53", false, true, "UNKNOWN", CodingEffect.NONE, 0.5);
        VariantContext var2 = createGermline("PASS", "KD53", false, false, "Pathogenic", CodingEffect.NONE, 0.5);
        victim.processVariant(var1);
        victim.processVariant(var2);
        victim.flush();

        assertFalse(new VariantContextDecorator(var1).reported());
        assertFalse(new VariantContextDecorator(var2).reported());
    }

    @Test
    public void testReportHotspotWhenMultipleHits()
    {
        DriverGene driverGene =
                createDriverGene("KD53", DriverGeneGermlineReporting.WILDTYPE_LOST, DriverGeneGermlineReporting.WILDTYPE_LOST);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Collections.emptySet());

        VariantContext var1 = createGermline("PASS", "KD53", false, true, "UNKNOWN", CodingEffect.NONE, 0.5);
        VariantContext var2 = createGermline("PASS", "KD53", false, false, "Pathogenic", CodingEffect.NONE, 0.5);
        victim.processVariant(var1);
        victim.processVariant(var2);
        victim.flush();

        assertTrue(new VariantContextDecorator(var1).reported());
        assertTrue(new VariantContextDecorator(var2).reported());
    }

    @Test
    public void testReportHotspotWhenMultipleHitsWithSameLPS()
    {
        DriverGene driverGene =
                createDriverGene("KD53", DriverGeneGermlineReporting.WILDTYPE_LOST, DriverGeneGermlineReporting.WILDTYPE_LOST);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Collections.emptySet());

        VariantContext var1 = createGermline("PASS", "KD53", false, true, "UNKNOWN", CodingEffect.NONE, 0.5, 4);
        VariantContext var2 =createGermline("PASS", "KD53", false, false, "Pathogenic", CodingEffect.NONE, 0.5, 4) ;
        victim.processVariant(var1);
        victim.processVariant(var2);
        victim.flush();

        assertFalse(new VariantContextDecorator(var1).reported());
        assertFalse(new VariantContextDecorator(var2).reported());
    }

    @Test
    public void testReportUnknownFrameshift()
    {
        DriverGene driverGene = createDriverGene("KD53", DriverGeneGermlineReporting.ANY, DriverGeneGermlineReporting.NONE);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Collections.emptySet());

        VariantContext var = createGermline("PASS", "KD53", false, false, "UNKNOWN", CodingEffect.NONSENSE_OR_FRAMESHIFT, 0.5);
        victim.processVariant(var);
        victim.flush();

        assertTrue(new VariantContextDecorator(var).reported());
    }

    @Test
    public void testReportPathogenic()
    {
        DriverGene driverGene = createDriverGene("KD53", DriverGeneGermlineReporting.ANY, DriverGeneGermlineReporting.NONE);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Collections.emptySet());

        VariantContext var = createGermline("PASS", "KD53", false, false, "Pathogenic", CodingEffect.NONE, 0.5);
        victim.processVariant(var);
        victim.flush();

        assertTrue(new VariantContextDecorator(var).reported());
    }

    @Test
    public void testReportBiallelicAsMultipleHit()
    {
        DriverGene driverGene = createDriverGene("KD53", DriverGeneGermlineReporting.ANY, DriverGeneGermlineReporting.WILDTYPE_LOST);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Collections.emptySet());

        VariantContext var = createGermline("PASS", "KD53", true, false, "Pathogenic", CodingEffect.NONE, 0.5);
        victim.processVariant(var);
        victim.flush();

        assertTrue(new VariantContextDecorator(var).reported());
    }

    @Test
    public void testIgnoreSinglePathogenicWhenNoMultipleHits()
    {
        DriverGene driverGene = createDriverGene("KD53", DriverGeneGermlineReporting.WILDTYPE_LOST, DriverGeneGermlineReporting.NONE);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Collections.emptySet());
        VariantContext var = createGermline("PASS", "KD53", false, false, "Pathogenic", CodingEffect.NONE, 0.5);
        victim.processVariant(var);
        victim.flush();

        assertFalse(new VariantContextDecorator(var).reported());
    }

    @Test
    public void testSomaticCountTowardsMultipleHits()
    {
        DriverGene driverGene = createDriverGene("KD53", DriverGeneGermlineReporting.WILDTYPE_LOST, DriverGeneGermlineReporting.NONE);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Sets.newHashSet("KD53"));

        VariantContext var = createGermline("PASS", "KD53", false, false, "Pathogenic", CodingEffect.NONE, 0.5);
        victim.processVariant(var);
        victim.flush();

        assertTrue(new VariantContextDecorator(var).reported());
    }

    @Test
    public void testVariantNotLost()
    {
        DriverGene driverGene = createDriverGene("KD53", DriverGeneGermlineReporting.VARIANT_NOT_LOST, DriverGeneGermlineReporting.NONE);
        GermlineReportedEnrichment victim = new GermlineReportedEnrichment(Lists.newArrayList(driverGene), Sets.newHashSet("KD53"));

        VariantContext var1 = createGermline("PASS", "KD53", false, false, "Pathogenic", CodingEffect.NONE, 0.4);
        VariantContext var2 = createGermline("PASS", "KD53", false, false, "Pathogenic", CodingEffect.NONE, 0.6);
        victim.processVariant(var1);
        victim.processVariant(var2);
        victim.flush();

        assertFalse(new VariantContextDecorator(var1).reported());
        assertTrue(new VariantContextDecorator(var2).reported());
    }

    @NotNull
    private static VariantContext createGermline(@NotNull String filter, @NotNull String gene, boolean biallelic, boolean isHotspot,
            @NotNull String clinSig, @NotNull CodingEffect codingEffect, double variantCopyNumber)
    {
        final String hotspotFlag = isHotspot ? HOTSPOT_FLAG : NEAR_HOTSPOT_FLAG;
        final String line =
                "11\t1000\tCOSM123;COSM456\tG\tA\t100\t" + filter + "\t" + PURPLE_BIALLELIC_FLAG + "=" + biallelic + ";" + hotspotFlag
                        + ";SEC=" + gene + ",ENST00000393562,UTR_variant," + codingEffect.toString() + ",c.-275T>G,;CLNSIG=" + clinSig + ";"
                        + PURPLE_VARIANT_CN_INFO + "=" + variantCopyNumber + ";\"\tGT:AD:DP\t0/1:73,17:91";
        return VariantContextFromString.decode(line);
    }

    @NotNull
    private static VariantContext createGermline(@NotNull String filter, @NotNull String gene, boolean biallelic, boolean isHotspot,
            @NotNull String clinSig, @NotNull CodingEffect codingEffect, double variantCopyNumber, int localPhaseSet)
    {
        final String hotspotFlag = isHotspot ? HOTSPOT_FLAG : NEAR_HOTSPOT_FLAG;
        final String line =
                "11\t1000\tCOSM123;COSM456\tG\tA\t100\t" + filter + "\t" + PURPLE_BIALLELIC_FLAG + "=" + biallelic + ";" + hotspotFlag
                        + ";SEC=" + gene + ",ENST00000393562,UTR_variant," + codingEffect.toString() + ",c.-275T>G,;CLNSIG=" + clinSig + ";"
                        + PURPLE_VARIANT_CN_INFO + "=" + variantCopyNumber + ";" + SageMetaData.LOCAL_PHASE_SET + "=" + localPhaseSet
                        + "\tGT:AD:DP\t0/1:73,17:91";
        return VariantContextFromString.decode(line);
    }

    @NotNull
    private static DriverGene createDriverGene(@NotNull String gene, @NotNull DriverGeneGermlineReporting variantReporting,
            @NotNull DriverGeneGermlineReporting hotspotReporting)
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
                .reportSomaticHotspot(false)
                .reportSplice(false)
                .reportGermlineDisruption(false)
                .build();
    }
}
