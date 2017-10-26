package com.hartwig.hmftools.strelka;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Streams;
import com.google.common.io.Resources;
import com.hartwig.hmftools.strelka.scores.ImmutableVariantScore;
import com.hartwig.hmftools.strelka.scores.ReadType;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class MNVDetectorTest {
    private static final File VCF_FILE = new File(Resources.getResource("mnvs.vcf").getPath());
    private static final VCFFileReader VCF_FILE_READER = new VCFFileReader(VCF_FILE, false);
    private static final List<VariantContext> VARIANTS = Streams.stream(VCF_FILE_READER).collect(Collectors.toList());

    // MIVO: variants at positions: 1  3(alts: A,T)  =>  possible mnvs: (1,3A) (1,3T) => actual mnv: (1,3A) => non-mnv variant: 3T
    @Test
    public void correctlyDetectsNonMnvVariant() {
        final PotentialMNVRegion region =
                PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(), Lists.newArrayList(VARIANTS.get(7), VARIANTS.get(9)));
        assertEquals(2, region.potentialMnvs().size());
        Map<PotentialMNV, MNVScore> scores = region.potentialMnvs()
                .stream()
                .collect(Collectors.toMap(Function.identity(), potentialMNV -> MNVScore.of(potentialMNV.variants())));
        final MNVScore score = TestUtils.build2VariantScores(region.potentialMnvs().get(0).variants(), Lists.newArrayList(
                ImmutablePair.of(ImmutableVariantScore.of(ReadType.ALT, 15), ImmutableVariantScore.of(ReadType.ALT, 15))));
        scores.put(region.potentialMnvs().get(0), score);
        assertTrue(scores.get(region.potentialMnvs().get(0)).isMNV());
        assertEquals(1, MNVDetector.nonMnvVariants(region, scores).size());
        final VariantContext nonMnv = MNVDetector.nonMnvVariants(region, scores).get(0);
        assertEquals("G", nonMnv.getReference().getBaseString());
        assertEquals(1, nonMnv.getAlternateAlleles().size());
        assertEquals("T", nonMnv.getAlternateAllele(0).getBaseString());
    }

    // MIVO: variants at positions: 1  3(alts: A,T)  =>  possible mnvs: (1,3A) (1,3T) => actual mnv: none => non-mnv variants: 1  3(alts: A,T)
    @Test
    public void correctlyDetectsMultiAltNonMnvVariant() {
        final PotentialMNVRegion region =
                PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(), Lists.newArrayList(VARIANTS.get(7), VARIANTS.get(9)));
        assertEquals(2, region.potentialMnvs().size());
        Map<PotentialMNV, MNVScore> scores = region.potentialMnvs()
                .stream()
                .collect(Collectors.toMap(Function.identity(), potentialMNV -> MNVScore.of(potentialMNV.variants())));
        assertEquals(2, MNVDetector.nonMnvVariants(region, scores).size());
        final List<VariantContext> nonMnvs = MNVDetector.nonMnvVariants(region, scores);
        assertEquals("C", nonMnvs.get(0).getReference().getBaseString());
        assertEquals(1, nonMnvs.get(0).getAlternateAlleles().size());
        assertEquals("T", nonMnvs.get(0).getAlternateAllele(0).getBaseString());
        assertEquals("G", nonMnvs.get(1).getReference().getBaseString());
        assertEquals(2, nonMnvs.get(1).getAlternateAlleles().size());
        assertEquals("A", nonMnvs.get(1).getAlternateAllele(0).getBaseString());
        assertEquals("T", nonMnvs.get(1).getAlternateAllele(1).getBaseString());
    }
}
