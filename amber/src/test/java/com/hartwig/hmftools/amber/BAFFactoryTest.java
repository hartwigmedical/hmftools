package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberApplication.DEFAULT_MAX_DEPTH_PERCENTAGE;
import static com.hartwig.hmftools.amber.AmberApplication.DEFAULT_MAX_HET_AF_PERCENTAGE;
import static com.hartwig.hmftools.amber.AmberApplication.DEFAULT_MIN_DEPTH_PERCENTAGE;
import static com.hartwig.hmftools.amber.AmberApplication.DEFAULT_MIN_HET_AF_PERCENTAGE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.pileup.PileupFile;

import org.junit.Test;

public class BAFFactoryTest {

    private static final double EPSILON = 1e-10;

    private final static BAFFactory VICTIM = new BAFFactory(DEFAULT_MIN_HET_AF_PERCENTAGE,
            DEFAULT_MAX_HET_AF_PERCENTAGE,
            DEFAULT_MIN_DEPTH_PERCENTAGE,
            DEFAULT_MAX_DEPTH_PERCENTAGE);

    private final static Pileup GOOD_NORMAL = PileupFile.fromString("seq2\t156\tA\t10\t...G..GGGG\t<975;:<<<<<");

    @Test
    public void testDepths() {
        final Pileup normal = PileupFile.fromString("seq2\t156\tA\t11\t.$......+2AG.+2AG.+2AGGG\t<975;:<<<<<");
        final Pileup tumor = PileupFile.fromString("seq2\t156\tA\t31\t.$......+2AG.+2AG.+2AGGG\t<975;:<<<<<");
        assertNotEquals(normal.readCount(), tumor.readCount());

        final AmberBAF victim = BAFFactory.create('A', normal, tumor);
        assertEquals(normal.readCount(), victim.normalDepth());
        assertEquals(tumor.readCount(), victim.tumorDepth());
    }


    @Test
    public void testWorking() {
        final Pileup tumor = PileupFile.fromString("seq2\t156\tA\t10\t.GGGGGGGGG\t<975;:<<<<<");

        final List<AmberBAF> result = VICTIM.create(Collections.singletonList(GOOD_NORMAL), Collections.singletonList(tumor));
        assertEquals(1, result.size());

        AmberBAF baf = result.get(0);
        assertEquals(10, baf.normalDepth());
        assertEquals(0.5, baf.normalBAF(), EPSILON) ;
        assertEquals(10, baf.tumorDepth());
        assertEquals(0.9, baf.tumorBAF(), EPSILON) ;
    }

    @Test
    public void testNoTumorReadCount() {
        final Pileup tumor = PileupFile.fromString("seq2\t156\tA\t0\t*\t*");

        final List<AmberBAF> result = VICTIM.create(Collections.singletonList(GOOD_NORMAL), Collections.singletonList(tumor));
        assertEquals(0, result.size());
    }


    @Test
    public void testInvalidBAF() {
        final Pileup tumor = PileupFile.fromString("seq2\t156\tA\t10\tCCCCCCCCCC\t*");

        final List<AmberBAF> result = VICTIM.create(Collections.singletonList(GOOD_NORMAL), Collections.singletonList(tumor));
        assertEquals(0, result.size());
    }

    @Test
    public void testUMC() {
        final Pileup normal = PileupFile.fromString("16\t76595025\tA\t1\t*\tf");
        final Pileup tumor = PileupFile.fromString("16\t76595025\tA\t1\t.\tm");
        final List<AmberBAF> result = VICTIM.create(Collections.singletonList(normal), Collections.singletonList(tumor));
        assertEquals(0, result.size());
    }

}
