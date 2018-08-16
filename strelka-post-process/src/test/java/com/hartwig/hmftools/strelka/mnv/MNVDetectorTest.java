package com.hartwig.hmftools.strelka.mnv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Streams;
import com.google.common.io.Resources;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class MNVDetectorTest {
    private static final File VCF_FILE = new File(Resources.getResource("mnvs.vcf").getPath());
    private static final VCFFileReader VCF_FILE_READER = new VCFFileReader(VCF_FILE, false);
    private static final List<VariantContext> VARIANTS = Streams.stream(VCF_FILE_READER).collect(Collectors.toList());

    @Test
    public void variantsWithin4BasesFormRegionFor4GapSize() {
        final PotentialMNVRegion region = PotentialMNVRegion.fromVariant(VARIANTS.get(28));
        final Pair<PotentialMNVRegion, Optional<PotentialMNVRegion>> output = MNVDetector.addMnvToRegion(region, VARIANTS.get(29), 4);
        assertTrue(MNVDetector.variantFitsRegion(region, VARIANTS.get(29), 4));
        assertEquals(Optional.empty(), output.getRight());
        assertEquals(2, output.getLeft().variants().size());
        assertEquals(1, output.getLeft().potentialMnvs().size());
        assertEquals(Lists.newArrayList(1221315, 1221316, 1221317, 1221318), output.getLeft().potentialMnvs().get(0).gapPositions());
    }

    @Test
    public void variantsWithin4BasesDoNotFormRegionFor3GapSize() {
        final PotentialMNVRegion region = PotentialMNVRegion.fromVariant(VARIANTS.get(28));
        final Pair<PotentialMNVRegion, Optional<PotentialMNVRegion>> output = MNVDetector.addMnvToRegion(region, VARIANTS.get(29), 3);
        assertFalse(MNVDetector.variantFitsRegion(region, VARIANTS.get(29), 3));
        assertEquals(1, output.getRight().get().variants().size());
        assertEquals(0, output.getRight().get().potentialMnvs().size());
        assertEquals(1, output.getLeft().variants().size());
        assertEquals(0, output.getLeft().potentialMnvs().size());
    }

    @Test
    public void consecutiveVariantsFormRegion() {
        final PotentialMNVRegion region = PotentialMNVRegion.fromVariant(VARIANTS.get(11));
        final Pair<PotentialMNVRegion, Optional<PotentialMNVRegion>> output = MNVDetector.addMnvToRegion(region, VARIANTS.get(12), 1);
        assertTrue(MNVDetector.variantFitsRegion(region, VARIANTS.get(12), 1));
        assertEquals(Optional.empty(), output.getRight());
        assertEquals(2, output.getLeft().variants().size());
        assertEquals(1, output.getLeft().potentialMnvs().size());
        assertEquals(Lists.newArrayList(), output.getLeft().potentialMnvs().get(0).gapPositions());
    }
}
