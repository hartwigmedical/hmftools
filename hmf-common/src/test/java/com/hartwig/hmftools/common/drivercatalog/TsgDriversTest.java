package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.ImmutableDndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.ImmutableDndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class TsgDriversTest {

    private DndsDriverGeneLikelihood geneLikelihood;
    private SomaticVariant missense;
    private SomaticVariant nonsense;
    private SomaticVariant indel;
    private SomaticVariant hotspot;
    private SomaticVariant biallelic;

    @Before
    public void setup() {
        DndsDriverImpactLikelihood missenseLikelihood = createLikelihood(0.0166790589295988, 9.89844130535209e-08);
        DndsDriverImpactLikelihood nonsenseLikelihood = createLikelihood(0.00730132326080717, 7.40370091523547e-09);
        DndsDriverImpactLikelihood spliceLikelihood = createLikelihood(0.00540540540540541, 1e-9);
        DndsDriverImpactLikelihood indelLikelihood = createLikelihood(0.0137214137214137, 1e-9);
        geneLikelihood = ImmutableDndsDriverGeneLikelihood.builder()
                .gene("TP53")
                .missense(missenseLikelihood)
                .nonsense(nonsenseLikelihood)
                .splice(spliceLikelihood)
                .indel(indelLikelihood)
                .build();

        missense = create(VariantType.SNP, CodingEffect.MISSENSE, false, 0.5);
        nonsense = create(VariantType.SNP, CodingEffect.NONSENSE_OR_FRAMESHIFT, false, 0.5);
        indel = create(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, false, 0.5);
        hotspot = create(VariantType.INDEL, CodingEffect.MISSENSE, true, 0.5);
        biallelic = create(VariantType.MNP, CodingEffect.NONSENSE_OR_FRAMESHIFT, false, 0.9);
    }

    @Test
    public void testHotspotFirst() {
        Map<VariantType, Long> counts = countMap(5161, 10000);
        DriverCatalog victim =
                TsgDrivers.geneDriver(geneLikelihood, Lists.newArrayList(hotspot, biallelic, missense), counts, counts, null);
        assertEquals(LikelihoodMethod.HOTSPOT, victim.likelihoodMethod());
        assertEquals(1, victim.driverLikelihood(), 0.01);
    }

    @Test
    public void testBiallelicSecond() {
        Map<VariantType, Long> counts = countMap(5161, 10000);
        DriverCatalog victim = TsgDrivers.geneDriver(geneLikelihood, Lists.newArrayList(biallelic, missense), counts, counts, null);
        assertEquals(LikelihoodMethod.BIALLELIC, victim.likelihoodMethod());
        assertEquals(1, victim.driverLikelihood(), 0.01);
    }

    @Test
    public void testSingleMissense() {
        Map<VariantType, Long> counts = countMap(351610, 10000);
        DriverCatalog victim = TsgDrivers.geneDriver(geneLikelihood, Lists.newArrayList(missense), counts, counts, null);
        assertEquals(LikelihoodMethod.DNDS, victim.likelihoodMethod());
        assertEquals(0.33, victim.driverLikelihood(), 0.01);
    }

    @Test
    public void testMultiMissense() {
        Map<VariantType, Long> counts = countMap(351610, 10000);
        DriverCatalog victim = TsgDrivers.geneDriver(geneLikelihood, Lists.newArrayList(missense, missense), counts, counts, null);
        assertEquals(LikelihoodMethod.DNDS, victim.likelihoodMethod());
        assertEquals(0.97, victim.driverLikelihood(), 0.01);
    }

    @Test
    public void testSingleNonsense() {
        Map<VariantType, Long> counts = countMap(351610, 10000);
        DriverCatalog victim = TsgDrivers.geneDriver(geneLikelihood, Lists.newArrayList(nonsense), counts, counts, null);
        assertEquals(LikelihoodMethod.DNDS, victim.likelihoodMethod());
        assertEquals(0.74, victim.driverLikelihood(), 0.01);
    }

    @Test
    public void testMixed() {
        Map<VariantType, Long> counts = countMap(351610, 10000);
        DriverCatalog victim = TsgDrivers.geneDriver(geneLikelihood, Lists.newArrayList(missense, nonsense), counts, counts, null);
        assertEquals(LikelihoodMethod.DNDS, victim.likelihoodMethod());
        assertEquals(1, victim.driverLikelihood(), 0.01);
    }

    @Test
    public void testIndel() {
        Map<VariantType, Long> counts = countMap(351610, 10000);
        DriverCatalog victim = TsgDrivers.geneDriver(geneLikelihood, Lists.newArrayList(indel), counts, counts, null);
        assertEquals(LikelihoodMethod.DNDS, victim.likelihoodMethod());
        assertEquals(1, victim.driverLikelihood(), 0.01);
    }

    @NotNull
    static DndsDriverImpactLikelihood createLikelihood(double pDriver, double pVariantNonDriver) {
        return ImmutableDndsDriverImpactLikelihood.builder()
                .driversPerSample(pDriver)
                .passengersPerMutation(pVariantNonDriver)
                .build();
    }

    @NotNull
    static SomaticVariant create(@NotNull VariantType type, @NotNull CodingEffect codingEffect, boolean hotspot, double vaf) {
        boolean biallelic = Doubles.greaterOrEqual(2 * vaf, 1.5);
        return SomaticVariantTestBuilderFactory.create()
                .type(type)
                .canonicalCodingEffect(codingEffect)
                .hotspot(hotspot ? Hotspot.HOTSPOT : Hotspot.NON_HOTSPOT)
                .adjustedCopyNumber(2)
                .adjustedVAF(vaf)
                .variantCopyNumber(2 * vaf)
                .biallelic(biallelic)
                .build();
    }

    @NotNull
    static Map<VariantType, Long> countMap(int snp, int indel) {
        Map<VariantType, Long> countMap = Maps.newHashMap();
        countMap.put(VariantType.SNP, (long) snp);
        countMap.put(VariantType.INDEL, (long) indel);
        return countMap;
    }
}
