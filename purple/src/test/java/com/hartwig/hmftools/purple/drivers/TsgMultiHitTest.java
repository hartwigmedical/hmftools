package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.common.variant.Hotspot.NON_HOTSPOT;
import static com.hartwig.hmftools.purple.MiscTestUtils.createVariant;
import static com.hartwig.hmftools.purple.drivers.OncoDriversTest.createGeneCopyNumber;
import static com.hartwig.hmftools.purple.drivers.TsgDriversTest.countMap;
import static com.hartwig.hmftools.purple.drivers.TsgDriversTest.createLikelihood;

import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.ImmutableDndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

import org.junit.Before;
import org.junit.Test;

public class TsgMultiHitTest {

    private DndsDriverGeneLikelihood geneLikelihood;

    @Before
    public void setup() {
        DndsDriverImpactLikelihood missenseLikelihood = createLikelihood(0.00458696050642149, 3.5842380991213754E-7);
        DndsDriverImpactLikelihood nonsenseLikelihood = createLikelihood(2.549386282565305E-4, 2.467222584510642E-8);
        DndsDriverImpactLikelihood spliceLikelihood = createLikelihood(2.0913513727648957E-4, 9.252084691914908E-9);
        DndsDriverImpactLikelihood indelLikelihood = createLikelihood(0.0020706936081327882, 2.4595609375053004E-7);

        geneLikelihood = ImmutableDndsDriverGeneLikelihood.builder()
                .gene("JAK1")
                .missense(missenseLikelihood)
                .nonsense(nonsenseLikelihood)
                .splice(spliceLikelihood)
                .indel(indelLikelihood)
                .build();
    }

    @Test
    public void testMultiHitIsNeverLessThanEquivalentMissense() {
        Map<VariantType,Integer> counts = countMap(83135, 241917);

        SomaticVariant missense = createVariant(VariantType.SNP, CodingEffect.MISSENSE, 0, NON_HOTSPOT, 0.5);
        SomaticVariant frameshift = createVariant(VariantType.INDEL, CodingEffect.NONSENSE_OR_FRAMESHIFT, 0, NON_HOTSPOT, 0.5);

        GeneCopyNumber geneCopyNumber = createGeneCopyNumber(GENE_NAME_1);

        DriverCatalog multiFrameshift = TsgDrivers.createTsgDriver(
                Lists.newArrayList(frameshift, frameshift), counts, counts, geneCopyNumber, geneLikelihood);

        DriverCatalog missenseAndFrameshift =
                TsgDrivers.createTsgDriver(Lists.newArrayList(missense, frameshift), counts, counts, geneCopyNumber, geneLikelihood);

        DriverCatalog singleFrameshift = TsgDrivers.createTsgDriver(
                Lists.newArrayList(frameshift), counts, counts, geneCopyNumber, geneLikelihood);

        DriverCatalog multiMissense = TsgDrivers.createTsgDriver(
                Lists.newArrayList(missense, missense), counts, counts, geneCopyNumber, geneLikelihood);

        assertTrue(Doubles.lessThan(singleFrameshift.driverLikelihood(), missenseAndFrameshift.driverLikelihood()));
        assertTrue(Doubles.equal(multiMissense.driverLikelihood(), missenseAndFrameshift.driverLikelihood()));
        assertTrue(Doubles.equal(multiMissense.driverLikelihood(), multiFrameshift.driverLikelihood()));
    }

}
