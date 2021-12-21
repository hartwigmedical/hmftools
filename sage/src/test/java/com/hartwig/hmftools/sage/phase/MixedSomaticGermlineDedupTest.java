package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.sage.phase.BufferedPostProcessorTest.create;
import static com.hartwig.hmftools.sage.phase.MixedSomaticGermlineDedup.codonDifferences;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.junit.Test;

public class MixedSomaticGermlineDedupTest
{

    @Test
    public void testCodonDifferences()
    {
        final BaseRegion codon = new BaseRegion(100, 102);

        final VariantHotspot snvBefore = create(99, "A", "T");
        final VariantHotspot snvOne = create(100, "A", "T");
        final VariantHotspot snvTwo = create(101, "A", "T");
        final VariantHotspot snvThree = create(102, "A", "T");
        final VariantHotspot snvAfter = create(103, "A", "T");

        assertEquals(0, codonDifferences(codon, snvBefore));
        assertEquals(0, codonDifferences(codon, snvAfter));

        assertEquals(1, codonDifferences(codon, snvOne));
        assertEquals(1, codonDifferences(codon, snvTwo));
        assertEquals(1, codonDifferences(codon, snvThree));

        assertEquals(1, codonDifferences(codon, create(102, "AA", "TT")));
        assertEquals(1, codonDifferences(codon, create(101, "AAA", "TAT")));
        assertEquals(2, codonDifferences(codon, create(100, "AAA", "TAT")));
        assertEquals(3, codonDifferences(codon, create(100, "AAA", "TTT")));
    }
}
