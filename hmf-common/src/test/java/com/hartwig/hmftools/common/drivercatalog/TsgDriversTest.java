package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.dnds.ImmutableDndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.dnds.ImmutableDndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class TsgDriversTest {

    private DndsDriverGeneLikelihood geneLikelihood;
    private EnrichedSomaticVariant missense;
    private EnrichedSomaticVariant nonsense;
    private EnrichedSomaticVariant indel;
    private EnrichedSomaticVariant hotspot;
    private EnrichedSomaticVariant biallelic;

    @Before
    public void setup() {
        DndsDriverImpactLikelihood missenseLikelihood = createLikelihood(0.872024711427939, 0.0166790589295988, 9.89844130535209e-08);
        DndsDriverImpactLikelihood nonsenseLikelihood = createLikelihood(0.975537913457847, 0.00730132326080717, 7.40370091523547e-09);
        DndsDriverImpactLikelihood spliceLikelihood = createLikelihood(1, 0.00540540540540541, 0);
        DndsDriverImpactLikelihood indelLikelihood = createLikelihood(1, 0.0137214137214137, 0);
        geneLikelihood = ImmutableDndsDriverGeneLikelihood.builder()
                .gene("TP53")
                .missense(missenseLikelihood)
                .nonsense(nonsenseLikelihood)
                .splice(spliceLikelihood)
                .indel(indelLikelihood)
                .build();

        missense = create(VariantType.SNP, CodingEffect.MISSENSE, false, false);
        nonsense = create(VariantType.SNP, CodingEffect.NONSENSE_OR_FRAMESHIFT, false, false);
        indel = create(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, false, false);
        hotspot = create(VariantType.INDEL, CodingEffect.MISSENSE, true, false);
        biallelic = create(VariantType.MNP, CodingEffect.NONSENSE_OR_FRAMESHIFT, false, true);
    }

    @Test
    public void testHotspotFirst() {
        final DriverCatalog victim = TsgDrivers.geneDriver(5161, 10000, geneLikelihood, Lists.newArrayList(hotspot, biallelic, missense));
        assertEquals(DriverType.HOTSPOT, victim.driver());
        assertEquals(1, victim.driverLikelihood(), 0.01);
        assertEquals(1, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void testBiallelicSecond() {
        final DriverCatalog victim = TsgDrivers.geneDriver(5161, 10000, geneLikelihood, Lists.newArrayList(biallelic, missense));
        assertEquals(DriverType.BIALLELIC, victim.driver());
        assertEquals(1, victim.driverLikelihood(), 0.01);
        assertEquals(0.98, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void testSingleMissense() {
        final DriverCatalog victim = TsgDrivers.geneDriver(351610, 10000, geneLikelihood, Lists.newArrayList(missense));
        assertEquals(DriverType.DNDS, victim.driver());
        assertEquals(0.33, victim.driverLikelihood(), 0.01);
        assertEquals(0.87, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void testMultiMissense() {
        final DriverCatalog victim = TsgDrivers.geneDriver(351610, 10000, geneLikelihood, Lists.newArrayList(missense, missense));
        assertEquals(DriverType.DNDS, victim.driver());
        assertEquals(0.97, victim.driverLikelihood(), 0.01);
        assertEquals(0.87, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void testSingleNonsense() {
        final DriverCatalog victim = TsgDrivers.geneDriver(351610, 10000, geneLikelihood, Lists.newArrayList(nonsense));
        assertEquals(DriverType.DNDS, victim.driver());
        assertEquals(0.74, victim.driverLikelihood(), 0.01);
        assertEquals(0.98, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void testMixed() {
        final DriverCatalog victim = TsgDrivers.geneDriver(351610, 10000, geneLikelihood, Lists.newArrayList(missense, nonsense));
        assertEquals(DriverType.DNDS, victim.driver());
        assertEquals(1, victim.driverLikelihood(), 0.01);
        assertEquals(0.98, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void testIndel() {
        final DriverCatalog victim = TsgDrivers.geneDriver(351610, 10000, geneLikelihood, Lists.newArrayList(indel));
        assertEquals(DriverType.DNDS, victim.driver());
        assertEquals(1, victim.driverLikelihood(), 0.01);
        assertEquals(1, victim.dndsLikelihood(), 0.01);
    }

    @NotNull
    static DndsDriverImpactLikelihood createLikelihood(double dndsLikelihood, double pDriver, double pVariantNonDriver) {
        return ImmutableDndsDriverImpactLikelihood.builder()
                .dndsLikelihood(dndsLikelihood)
                .pDriver(pDriver)
                .pVariantNonDriverFactor(pVariantNonDriver)
                .build();
    }

    @NotNull
    private EnrichedSomaticVariant create(@NotNull final VariantType type, @NotNull final CodingEffect codingEffect, boolean hotspot,
            boolean biallelic) {
        return SomaticVariantTestBuilderFactory.createEnriched()
                .type(type)
                .canonicalCodingEffect(codingEffect)
                .hotspot(hotspot)
                .biallelic(biallelic)
                .build();
    }

}
