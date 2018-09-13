package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.dnds.DndsDriverLikelihood;
import com.hartwig.hmftools.common.dnds.DndsDriverLikelihoodSupplier;
import com.hartwig.hmftools.common.dnds.ImmutableDndsDriverLikelihood;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class OncoDriversTest {

    private static final double EPSILON = 0.0001;
    private static final double UNADJUSTED_LIKELIHOOD = 0.5;

    private EnrichedSomaticVariant frameshiftHotspot;

    private EnrichedSomaticVariant inframe;
    private EnrichedSomaticVariant invalidInframe;

    private EnrichedSomaticVariant frameshift;
    private EnrichedSomaticVariant missense;

    private DndsDriverLikelihood likelihood;

    @Before
    public void setup() {
        likelihood = createLikelihood();
        frameshiftHotspot = create(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, 3, true);

        inframe = create(VariantType.INDEL, CodingEffect.MISSENSE, 3, false);
        invalidInframe = create(VariantType.INDEL, CodingEffect.MISSENSE, OncoDrivers.MAX_REPEAT_COUNT + 1, false);
        frameshift = create(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, 3, false);
        missense = create(VariantType.SNP, CodingEffect.MISSENSE, 0, false);
    }

    @Test
    public void testHIST2H3D() throws IOException {
        final Map<String, DndsDriverLikelihood> dnds = DndsDriverLikelihoodSupplier.oncoLikelihood();
        double value = OncoDrivers.missenseProbabilityDriverVariant(27742, dnds.get("HIST2H3D"));
        assertEquals(0.7042, value, EPSILON);
    }

    @Test
    public void testABL1() throws IOException {
        final Map<String, DndsDriverLikelihood> dnds = DndsDriverLikelihoodSupplier.oncoLikelihood();
        double value = OncoDrivers.missenseProbabilityDriverVariant(996698, dnds.get("ABL1"));
        assertEquals(0.0057, value, EPSILON);
    }

    @Test
    public void favourHotspot() {
        final DriverCatalog victim = OncoDrivers.geneDriver(10000, likelihood, Lists.newArrayList(frameshiftHotspot, inframe, missense));
        assertEquals(DriverType.HOTSPOT, victim.driver());
        assertEquals(1, victim.driverLikelihood(), 0.01);
        assertEquals(UNADJUSTED_LIKELIHOOD, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void favourValidInframe() {
        final DriverCatalog victim = OncoDrivers.geneDriver(10000, likelihood, Lists.newArrayList(inframe, missense));
        assertEquals(DriverType.INFRAME, victim.driver());
        assertEquals(1, victim.driverLikelihood(), 0.01);
        assertEquals(UNADJUSTED_LIKELIHOOD, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void ignoreFrameshiftAndInvalidInframe() {
        final DriverCatalog victim = OncoDrivers.geneDriver(10000, likelihood, Lists.newArrayList(invalidInframe, frameshift));
        assertEquals(DriverType.SINGLE_HIT, victim.driver());
        assertEquals(0, victim.driverLikelihood(), 0.01);
        assertEquals(0, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void singleMissense() {
        final DriverCatalog victim = OncoDrivers.geneDriver(10000, likelihood, Lists.newArrayList(missense));
        assertEquals(DriverType.SINGLE_HIT, victim.driver());
        assertEquals(0.68, victim.driverLikelihood(), 0.01);
        assertEquals(0.5, victim.dndsLikelihood(), 0.01);
    }

    @Test
    public void multiMissense() {
        final DriverCatalog victim = OncoDrivers.geneDriver(10000, likelihood, Lists.newArrayList(missense, missense));
        assertEquals(DriverType.MULTI_HIT, victim.driver());
        assertEquals(0.68, victim.driverLikelihood(), 0.01);
        assertEquals(0.5, victim.dndsLikelihood(), 0.01);
    }

    @NotNull
    private DndsDriverLikelihood createLikelihood() {
        return ImmutableDndsDriverLikelihood.builder()
                .gene("GENE")
                .missenseUnadjustedDriverLikelihood(UNADJUSTED_LIKELIHOOD)
                .missenseProbabilityDriver(0.002)
                .missenseProbabilityVariantNonDriverFactor(9.26e-08)
                .build();
    }

    @NotNull
    private EnrichedSomaticVariant create(@NotNull final VariantType type, @NotNull final CodingEffect codingEffect, int repeatCount,
            boolean hotspot) {
        return SomaticVariantTestBuilderFactory.createEnriched()
                .type(type)
                .canonicalCodingEffect(codingEffect)
                .repeatCount(repeatCount)
                .hotspot(hotspot)
                .build();
    }

}
