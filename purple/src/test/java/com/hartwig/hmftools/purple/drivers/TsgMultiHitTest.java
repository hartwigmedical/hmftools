package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.variant.HotspotType.NON_HOTSPOT;
import static com.hartwig.hmftools.purple.MiscTestUtils.createVariant;
import static com.hartwig.hmftools.purple.drivers.OncoDriversTest.createGeneCopyNumber;
import static com.hartwig.hmftools.purple.drivers.TsgDriversTest.countMap;
import static com.hartwig.hmftools.purple.drivers.TsgDriversTest.createLikelihood;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.driver.dnds.ImmutableDndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import org.junit.Before;
import org.junit.Test;

public class TsgMultiHitTest {

    private DndsDriverGeneLikelihood geneLikelihood;

    @Before
    public void setup()
    {
        DndsDriverImpactLikelihood missenseLikelihood = createLikelihood(0.0045, 3.6E-7);
        DndsDriverImpactLikelihood nonsenseLikelihood = createLikelihood(2.5E-4, 2.5E-8);
        DndsDriverImpactLikelihood spliceLikelihood = createLikelihood(2.1E-4, 9.2E-9);
        DndsDriverImpactLikelihood indelLikelihood = createLikelihood(0.002, 2.5E-7);

        geneLikelihood = ImmutableDndsDriverGeneLikelihood.builder()
                .gene("JAK1")
                .missense(missenseLikelihood)
                .nonsense(nonsenseLikelihood)
                .splice(spliceLikelihood)
                .indel(indelLikelihood)
                .build();
    }

    @Test
    public void testMultiHitIsNeverLessThanEquivalentMissense()
    {
        Map<VariantType,Integer> counts = countMap(83135, 241917);

        SomaticVariant missense = createVariant(VariantType.SNP, CodingEffect.MISSENSE, 0, NON_HOTSPOT, 0.5);
        SomaticVariant frameshift = createVariant(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, 0, NON_HOTSPOT, 0.5);

        GeneCopyNumber geneCopyNumber = createGeneCopyNumber(GENE_NAME_1);

        DriverCatalog multiFrameshift = TsgDrivers.buildDriverCatalog(
                Lists.newArrayList(frameshift, frameshift), counts, counts, geneCopyNumber, geneLikelihood,
                LikelihoodMethod.DNDS, ReportedStatus.REPORTED);

        DriverCatalog missenseAndFrameshift = TsgDrivers.buildDriverCatalog(
                        Lists.newArrayList(missense, frameshift), counts, counts, geneCopyNumber, geneLikelihood,
                        LikelihoodMethod.DNDS, ReportedStatus.REPORTED);

        DriverCatalog singleFrameshift = TsgDrivers.buildDriverCatalog(
                Lists.newArrayList(frameshift), counts, counts, geneCopyNumber, geneLikelihood,
                LikelihoodMethod.DNDS, ReportedStatus.REPORTED);

        DriverCatalog multiMissense = TsgDrivers.buildDriverCatalog(
                Lists.newArrayList(missense, missense), counts, counts, geneCopyNumber, geneLikelihood,
                LikelihoodMethod.DNDS, ReportedStatus.REPORTED);

        assertTrue(Doubles.lessThan(singleFrameshift.driverLikelihood(), missenseAndFrameshift.driverLikelihood()));
        assertTrue(Doubles.equal(multiMissense.driverLikelihood(), missenseAndFrameshift.driverLikelihood()));
        assertTrue(Doubles.equal(multiMissense.driverLikelihood(), multiFrameshift.driverLikelihood()));
    }

    @Test
    public void testMultiHitWithNonImpactingVariant()
    {
        // a variant in the splice region can get a driver catalog record but not be reported nor have a coding impact
        // and should then not affect DNDS multi-hit logic
        Map<VariantType,Integer> counts = countMap(83135, 241917);

        SomaticVariant missense = createVariant(VariantType.SNP, CodingEffect.MISSENSE, 0, NON_HOTSPOT, 0.5);
        SomaticVariant splaceIntronVariant = createVariant(VariantType.SNP, CodingEffect.UNDEFINED, 0, NON_HOTSPOT, 0.5);

        GeneCopyNumber geneCopyNumber = createGeneCopyNumber(GENE_NAME_1);

        DriverCatalog missensePlusNoImpact = TsgDrivers.buildDriverCatalog(
                Lists.newArrayList(missense, splaceIntronVariant), counts, counts, geneCopyNumber, geneLikelihood,
                LikelihoodMethod.DNDS, ReportedStatus.REPORTED);

        assertEquals(0.13, missensePlusNoImpact.driverLikelihood(), 0.1);
    }
}
