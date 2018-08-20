package com.hartwig.hmftools.strelka.mnv;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Streams;
import com.google.common.io.Resources;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class PotentialMNVRegionTest {
    private static final int GAP_SIZE = 1;
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
        final PotentialMNVRegion twoVariantRegion = PotentialMNVRegion.addVariant(oneVariantRegion, CONSECUTIVE_SNP2, GAP_SIZE);
        final PotentialMNVRegion threeVariantRegion = PotentialMNVRegion.addVariant(twoVariantRegion, CONSECUTIVE_SNP3, GAP_SIZE);
        assertEquals(3, twoVariantRegion.mnvs().size());
        assertEquals(1, twoVariantRegion.potentialMnvs().size());
        assertEquals(7, threeVariantRegion.mnvs().size());
        assertEquals(4, threeVariantRegion.potentialMnvs().size());
        assertEquals(Lists.newArrayList(CONSECUTIVE_SNP1, CONSECUTIVE_SNP2, CONSECUTIVE_SNP3), threeVariantRegion.variants());
    }

    // MIVO: variants at positions: 1  2  4  =>  possible mnvs: (1,2) (2,4) (1,2,4)
    @Test
    public void correctlyMakesCombinationsOfVariants() {
        final PotentialMNVRegion oneVariantRegion = PotentialMNVRegion.fromVariant(CONSECUTIVE_SNP1);
        final PotentialMNVRegion twoVariantRegion = PotentialMNVRegion.addVariant(oneVariantRegion, CONSECUTIVE_SNP2, GAP_SIZE);
        final PotentialMNVRegion threeVariantRegion = PotentialMNVRegion.addVariant(twoVariantRegion, CONSECUTIVE_SNP4, GAP_SIZE);
        assertEquals(3, twoVariantRegion.mnvs().size());
        assertEquals(1, twoVariantRegion.potentialMnvs().size());
        assertEquals(6, threeVariantRegion.mnvs().size());
        assertEquals(3, threeVariantRegion.potentialMnvs().size());
        assertEquals(Lists.newArrayList(CONSECUTIVE_SNP1, CONSECUTIVE_SNP2, CONSECUTIVE_SNP4), threeVariantRegion.variants());
    }

    // MIVO: variants at positions: 1  3(alts: A,T)  =>  possible mnvs: (1,3A) (1,3T)
    @Test
    public void correctlyMakesCombinationsWithMultipleAlts() {
        final PotentialMNVRegion oneVariantRegion = PotentialMNVRegion.fromVariant(SNP);
        final PotentialMNVRegion twoVariantRegion = PotentialMNVRegion.addVariant(oneVariantRegion, MULTI_ALT, GAP_SIZE);
        assertEquals(5, twoVariantRegion.mnvs().size());
        assertEquals(2, twoVariantRegion.potentialMnvs().size());
        assertEquals(Lists.newArrayList(SNP, MULTI_ALT), twoVariantRegion.variants());
    }

    // MIVO: variants at positions: 1  1(del of 1)  4  =>  possible mnvs: (1del,4)
    @Test
    public void correctlyMakesCombinationsWithSNPAfterDEL() {
        final PotentialMNVRegion oneVariantRegion = PotentialMNVRegion.fromVariant(SNP);
        final PotentialMNVRegion twoVariantRegion = PotentialMNVRegion.addVariant(oneVariantRegion, DEL, GAP_SIZE);
        final PotentialMNVRegion threeVariantRegion = PotentialMNVRegion.addVariant(twoVariantRegion, SNP_AFTER_DEL, GAP_SIZE);
        assertEquals(2, twoVariantRegion.mnvs().size());
        assertEquals(0, twoVariantRegion.potentialMnvs().size());
        assertEquals(4, threeVariantRegion.mnvs().size());
        assertEquals(1, threeVariantRegion.potentialMnvs().size());
        assertEquals(Lists.newArrayList(SNP, DEL, SNP_AFTER_DEL), threeVariantRegion.variants());
    }

    // MIVO: variants at positions: 1  1(del of 1)  3(alts: A,T)  4  =>
    //      possible mnvs: (1,3A) (1,3T) (1,3A,4) (1,3T,4) (1del,3A) (1del,3T) (1del,4) (1del,3A,4) (1del,3T,4) (3A,4) (3T,4)
    @Test
    public void correctlyMakesCombinationsWithSNPAfterDELAndMultiAlt() {
        final PotentialMNVRegion region = PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(),
                Lists.newArrayList(SNP, DEL, MULTI_ALT, SNP_AFTER_DEL),
                GAP_SIZE);
        assertEquals(16, region.mnvs().size());
        assertEquals(11, region.potentialMnvs().size());
        assertEquals(Lists.newArrayList(SNP, DEL, MULTI_ALT, SNP_AFTER_DEL), region.variants());
    }
}
