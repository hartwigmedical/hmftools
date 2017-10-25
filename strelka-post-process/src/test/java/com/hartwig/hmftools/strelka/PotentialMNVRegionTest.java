package com.hartwig.hmftools.strelka;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Streams;
import com.google.common.io.Resources;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class PotentialMNVRegionTest {
    private static final Logger LOGGER = LogManager.getLogger(PotentialMNVRegionTest.class);

    private static final File VCF_FILE = new File(Resources.getResource("mnvs.vcf").getPath());
    private static final VCFFileReader VCF_FILE_READER = new VCFFileReader(VCF_FILE, false);
    private static final List<VariantContext> VARIANTS = Streams.stream(VCF_FILE_READER).collect(Collectors.toList());
    private static final VariantContext SNP = VARIANTS.get(7);
    private static final VariantContext DEL = VARIANTS.get(8);
    private static final VariantContext MULTI_ALT = VARIANTS.get(9);
    private static final VariantContext SNP_AFTER_DEL = VARIANTS.get(10);
    private static final VariantContext CONSECUTIVE_SNP1 = VARIANTS.get(11);
    private static final VariantContext CONSECUTIVE_SNP2 = VARIANTS.get(12);
    private static final VariantContext CONSECUTIVE_SNP3 = VARIANTS.get(13);
    private static final VariantContext CONSECUTIVE_SNP4 = VARIANTS.get(14);

    // MIVO: variants at positions: 1  2  3  =>  possible mnvs: (1,2) (1,3) (2,3) (1,2,3)
    @Test
    public void correctlyMakesCombinationsOfConsecutiveVariants() {
        final PotentialMNVRegion oneVariantRegion = PotentialMNVRegion.fromVariant(CONSECUTIVE_SNP1);
        final PotentialMNVRegion twoVariantRegion = PotentialMNVRegion.addVariant(oneVariantRegion, CONSECUTIVE_SNP2);
        final PotentialMNVRegion threeVariantRegion = PotentialMNVRegion.addVariant(twoVariantRegion, CONSECUTIVE_SNP3);
        assertEquals(3, twoVariantRegion.mnvs().size());
        assertEquals(1, twoVariantRegion.potentialMnvs().size());
        assertEquals(7, threeVariantRegion.mnvs().size());
        assertEquals(4, threeVariantRegion.potentialMnvs().size());
    }

    // MIVO: variants at positions: 1  2  4  =>  possible mnvs: (1,2) (2,4) (1,2,4)
    @Test
    public void correctlyMakesCombinationsOfVariants() {
        final PotentialMNVRegion oneVariantRegion = PotentialMNVRegion.fromVariant(CONSECUTIVE_SNP1);
        final PotentialMNVRegion twoVariantRegion = PotentialMNVRegion.addVariant(oneVariantRegion, CONSECUTIVE_SNP2);
        final PotentialMNVRegion threeVariantRegion = PotentialMNVRegion.addVariant(twoVariantRegion, CONSECUTIVE_SNP4);
        assertEquals(3, twoVariantRegion.mnvs().size());
        assertEquals(1, twoVariantRegion.potentialMnvs().size());
        assertEquals(6, threeVariantRegion.mnvs().size());
        assertEquals(3, threeVariantRegion.potentialMnvs().size());
    }

    // MIVO: variants at positions: 1  3(alts: A,T)  =>  possible mnvs: (1,3A) (1,3T)
    @Test
    public void correctlyMakesCombinationsWithMultipleAlts() {
        final PotentialMNVRegion oneVariantRegion = PotentialMNVRegion.fromVariant(SNP);
        final PotentialMNVRegion twoVariantRegion = PotentialMNVRegion.addVariant(oneVariantRegion, MULTI_ALT);
        assertEquals(5, twoVariantRegion.mnvs().size());
        assertEquals(2, twoVariantRegion.potentialMnvs().size());
    }

    // MIVO: variants at positions: 1  1(del of 1)  4  =>  possible mnvs: (1del,4)
    @Test
    public void correctlyMakesCombinationsWithSNPafterDEL() {
        final PotentialMNVRegion oneVariantRegion = PotentialMNVRegion.fromVariant(SNP);
        final PotentialMNVRegion twoVariantRegion = PotentialMNVRegion.addVariant(oneVariantRegion, DEL);
        final PotentialMNVRegion threeVariantRegion = PotentialMNVRegion.addVariant(twoVariantRegion, SNP_AFTER_DEL);
        assertEquals(2, twoVariantRegion.mnvs().size());
        assertEquals(0, twoVariantRegion.potentialMnvs().size());
        assertEquals(4, threeVariantRegion.mnvs().size());
        assertEquals(1, threeVariantRegion.potentialMnvs().size());
    }
}
