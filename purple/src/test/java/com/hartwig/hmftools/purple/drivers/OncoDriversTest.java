package com.hartwig.hmftools.purple.drivers;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFactory;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.ImmutableDndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.ModifiableDndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.test.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class OncoDriversTest {

    private static final int SNV_SAMPLE_COUNT = 10000;
    private static final int INDEL_SAMPLE_COUNT = 1000;
    private static final double PASSENGERS_PER_MUTATION = 9.26e-08;

    private SomaticVariant frameshiftHotspot;
    private SomaticVariant frameshiftNearHotspot;
    private SomaticVariant inframe;
    private SomaticVariant unKnownInframe;
    private SomaticVariant frameshift;
    private SomaticVariant missense;
    private DndsDriverGeneLikelihood likelihood;

    @Before
    public void setup() {
        likelihood = createGeneLikelihood(0.002, 0.003, 0.002, 0.001);
        frameshiftHotspot = create(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, 3, Hotspot.HOTSPOT);
        frameshiftNearHotspot = create(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, 3, Hotspot.NEAR_HOTSPOT);

        inframe = create(VariantType.INDEL, CodingEffect.MISSENSE, 3, Hotspot.NON_HOTSPOT);
        unKnownInframe = create(VariantType.INDEL, CodingEffect.MISSENSE, OncoDrivers.MAX_REPEAT_COUNT + 1, Hotspot.NON_HOTSPOT);
        frameshift = create(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, 3, Hotspot.NON_HOTSPOT);
        missense = create(VariantType.SNP, CodingEffect.MISSENSE, 0, Hotspot.NON_HOTSPOT);
    }

    @Test
    public void favourHotspot() {
        DriverCatalog victim = OncoDrivers.geneDriver(10000,
                0,
                "Gene",
                likelihood,
                Lists.newArrayList(frameshiftHotspot, frameshiftNearHotspot, inframe, missense),
                null);
        Assert.assertEquals(LikelihoodMethod.HOTSPOT, victim.likelihoodMethod());
        assertEquals(1, victim.driverLikelihood(), 0.01);
    }

    @Test
    public void doNotFavourNearHotspot() {
        DriverCatalog victim =
                OncoDrivers.geneDriver(10000, 0, "Gene", likelihood, Lists.newArrayList(frameshiftNearHotspot, inframe, missense), null);
        assertEquals(LikelihoodMethod.INFRAME, victim.likelihoodMethod());
        assertEquals(1, victim.driverLikelihood(), 0.01);
    }

    @Test
    public void favourValidInframe() {
        DriverCatalog victim = OncoDrivers.geneDriver(10000, 0, "Gene", likelihood, Lists.newArrayList(inframe, missense), null);
        assertEquals(LikelihoodMethod.INFRAME, victim.likelihoodMethod());
        assertEquals(1, victim.driverLikelihood(), 0.01);
    }

    @Test
    public void favourHighest() {
        DriverCatalog victim = OncoDrivers.geneDriver(SNV_SAMPLE_COUNT,
                INDEL_SAMPLE_COUNT,
                "Gene",
                likelihood,
                Lists.newArrayList(frameshift, missense),
                null);
        assertEquals(LikelihoodMethod.DNDS, victim.likelihoodMethod());
        assertLikelihood(0.001, INDEL_SAMPLE_COUNT, victim.driverLikelihood());
    }

    @Test
    public void ignoreFrameshiftAndInvalidInframe() {
        DriverCatalog victim =
                OncoDrivers.geneDriver(SNV_SAMPLE_COUNT, INDEL_SAMPLE_COUNT, "Gene", likelihood, Lists.newArrayList(unKnownInframe), null);
        assertEquals(LikelihoodMethod.DNDS, victim.likelihoodMethod());
        assertLikelihood(0.001, INDEL_SAMPLE_COUNT, victim.driverLikelihood());
    }

    @Test
    public void singleMissense() {
        DriverCatalog victim =
                OncoDrivers.geneDriver(SNV_SAMPLE_COUNT, INDEL_SAMPLE_COUNT, "Gene", likelihood, Lists.newArrayList(missense), null);
        assertEquals(LikelihoodMethod.DNDS, victim.likelihoodMethod());
        assertLikelihood(0.002, SNV_SAMPLE_COUNT, victim.driverLikelihood());
    }

    @Test
    public void multiMissense() {
        DriverCatalog victim = OncoDrivers.geneDriver(SNV_SAMPLE_COUNT,
                INDEL_SAMPLE_COUNT,
                "Gene",
                likelihood,
                Lists.newArrayList(missense, missense),
                null);
        assertEquals(LikelihoodMethod.DNDS, victim.likelihoodMethod());
        assertLikelihood(0.002, SNV_SAMPLE_COUNT, victim.driverLikelihood());
    }

    @NotNull
    private static DndsDriverImpactLikelihood createLikelihood(double pDriver) {
        return ImmutableDndsDriverImpactLikelihood.builder()
                .driversPerSample(pDriver)
                .passengersPerMutation(PASSENGERS_PER_MUTATION)
                .build();
    }

    @NotNull
    private static SomaticVariant create(@NotNull VariantType type, @NotNull CodingEffect codingEffect, int repeatCount,
            @NotNull Hotspot hotspot) {
        return SomaticVariantTestBuilderFactory.create()
                .type(type)
                .canonicalCodingEffect(codingEffect)
                .repeatCount(repeatCount)
                .hotspot(hotspot)
                .biallelic(false)
                .build();
    }

    @NotNull
    private static ModifiableDndsDriverGeneLikelihood createGeneLikelihood(double missense, double nonsense, double splice, double indel) {
        final DndsDriverImpactLikelihood missenseLikelihood = createLikelihood(missense);
        return ModifiableDndsDriverGeneLikelihood.create()
                .setGene("Gene")
                .setMissense(missenseLikelihood)
                .setNonsense(createLikelihood(nonsense))
                .setSplice(createLikelihood(splice))
                .setIndel(createLikelihood(indel));
    }

    private static void assertLikelihood(double expectedProbability, int sampleCount, double value) {
        double expectedLikelihood = DriverCatalogFactory.probabilityDriverVariant(sampleCount, createLikelihood(expectedProbability));
        assertEquals(expectedLikelihood, value, 1e-10);
    }
}
