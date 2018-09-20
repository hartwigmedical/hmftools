package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.dnds.ImmutableDndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class OncoDriversTest {

    private static final double EPSILON = 0.0001;
    private static final double UNADJUSTED_LIKELIHOOD = 0.5;

    private EnrichedSomaticVariant frameshiftHotspot;
    private EnrichedSomaticVariant frameshiftNearHotspot;
    private EnrichedSomaticVariant inframe;
    private EnrichedSomaticVariant invalidInframe;
    private EnrichedSomaticVariant frameshift;
    private EnrichedSomaticVariant missense;
    private DndsDriverImpactLikelihood likelihood;

    private OncoDrivers oncoDrivers;

    @Before
    public void setup() {
        oncoDrivers = new OncoDrivers();

        likelihood = createLikelihood();
        frameshiftHotspot = create(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, 3, Hotspot.HOTSPOT);
        frameshiftNearHotspot = create(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, 3, Hotspot.NEAR_HOTSPOT);

        inframe = create(VariantType.INDEL, CodingEffect.MISSENSE, 3, Hotspot.NON_HOTSPOT);
        invalidInframe = create(VariantType.INDEL, CodingEffect.MISSENSE, OncoDrivers.MAX_REPEAT_COUNT + 1, Hotspot.NON_HOTSPOT);
        frameshift = create(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, 3, Hotspot.NON_HOTSPOT);
        missense = create(VariantType.SNP, CodingEffect.MISSENSE, 0, Hotspot.NON_HOTSPOT);
    }

    @Test
    public void testHIST2H3D() throws IOException {
        final Map<String, DndsDriverGeneLikelihood> dnds = DndsDriverGeneLikelihoodSupplier.oncoLikelihood();
        double value = OncoDrivers.missenseProbabilityDriverVariant(27742, dnds.get("HIST2H3D").missense());
        assertEquals(0.7042, value, EPSILON);
    }

    @Test
    public void testABL1() throws IOException {
        final Map<String, DndsDriverGeneLikelihood> dnds = DndsDriverGeneLikelihoodSupplier.oncoLikelihood();
        double value = OncoDrivers.missenseProbabilityDriverVariant(996698, dnds.get("ABL1").missense());
        assertEquals(0.0057, value, EPSILON);
    }

    @Test
    public void favourHotspot() {
        final DriverCatalog victim =
                oncoDrivers.geneDriver(10000, "Gene", likelihood, Lists.newArrayList(frameshiftHotspot, frameshiftNearHotspot, inframe, missense));
        assertEquals(DriverType.HOTSPOT, victim.driver());
        assertEquals(1, victim.driverLikelihood(), 0.01);
        assertEquals(UNADJUSTED_LIKELIHOOD, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void favourNearHotspot() {
        final DriverCatalog victim =
                oncoDrivers.geneDriver(10000, "Gene", likelihood, Lists.newArrayList(frameshiftNearHotspot, inframe, missense));
        assertEquals(DriverType.NEAR_HOTSPOT, victim.driver());
        assertEquals(1, victim.driverLikelihood(), 0.01);
        assertEquals(UNADJUSTED_LIKELIHOOD, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void favourValidInframe() {
        final DriverCatalog victim = oncoDrivers.geneDriver(10000, "Gene", likelihood, Lists.newArrayList(inframe, missense));
        assertEquals(DriverType.INFRAME, victim.driver());
        assertEquals(1, victim.driverLikelihood(), 0.01);
        assertEquals(UNADJUSTED_LIKELIHOOD, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void ignoreFrameshiftAndInvalidInframe() {
        final DriverCatalog victim = oncoDrivers.geneDriver(10000, "Gene", likelihood, Lists.newArrayList(invalidInframe, frameshift));
        assertEquals(DriverType.NONE, victim.driver());
        assertEquals(0, victim.driverLikelihood(), 0.01);
        assertEquals(0, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void singleMissense() {
        final DriverCatalog victim = oncoDrivers.geneDriver(10000, "Gene", likelihood, Lists.newArrayList(missense));
        assertEquals(DriverType.DNDS, victim.driver());
        assertEquals(0.68, victim.driverLikelihood(), 0.01);
        assertEquals(0.5, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void multiMissense() {
        final DriverCatalog victim = oncoDrivers.geneDriver(10000, "Gene", likelihood, Lists.newArrayList(missense, missense));
        assertEquals(DriverType.DNDS, victim.driver());
        assertEquals(0.68, victim.driverLikelihood(), 0.01);
        assertEquals(0.5, victim.dndsLikelihood(), 0.01);
    }

    @NotNull
    private DndsDriverImpactLikelihood createLikelihood() {
        return ImmutableDndsDriverImpactLikelihood.builder()
                .dndsLikelihood(UNADJUSTED_LIKELIHOOD)
                .pDriver(0.002)
                .pVariantNonDriverFactor(9.26e-08)
                .build();
    }

    @NotNull
    private EnrichedSomaticVariant create(@NotNull final VariantType type, @NotNull final CodingEffect codingEffect, int repeatCount,
            Hotspot hotspot) {
        return SomaticVariantTestBuilderFactory.createEnriched()
                .type(type)
                .canonicalCodingEffect(codingEffect)
                .repeatCount(repeatCount)
                .hotspot(hotspot)
                .build();
    }

}
