package com.hartwig.hmftools.common.variant.recovery;

import static com.hartwig.hmftools.common.variant.recovery.RecoverStructuralVariants.UNBALANCED_MIN_DEPTH_WINDOW_COUNT;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurityAdjusterTypicalChromosome;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import mockit.Injectable;

public class RecoverStructuralVariantsTest {

    @Injectable
    private RecoveredVariantFactory recoveredVariantFactory;
    private final PurityAdjuster purityAdjuster = new PurityAdjusterTypicalChromosome(Gender.FEMALE, 1, 0.68);

    @Test
    public void testRecoverUnbalancedSingle() throws IOException {
        List<VariantContext> result = singleTest(UNBALANCED_MIN_DEPTH_WINDOW_COUNT, 3);
        assertFalse(result.isEmpty());
        assertEquals(10000, result.get(0).getStart());
    }

    @Test
    public void testSingleIsAlreadyBalanced() throws IOException {
        List<VariantContext> result = singleTest(UNBALANCED_MIN_DEPTH_WINDOW_COUNT, 0.8);
        assertTrue(result.isEmpty());
    }

    @Test
    public void testSingleHasInsufficientDepth() throws IOException {
        List<VariantContext> result = singleTest(UNBALANCED_MIN_DEPTH_WINDOW_COUNT - 1, 3);
        assertTrue(result.isEmpty());
    }

    @Test
    public void testDelEndIsUnbalanced() throws IOException {
        final List<VariantContext> result = delTest(3, 0.8, 0.8);
        assertFalse(result.isEmpty());
        assertEquals(20000, result.get(0).getStart());
    }

    @Test
    public void testDelStartIsUnbalanced() throws IOException {
        final List<VariantContext> result = delTest(0.8, 0.8, 3);
        assertFalse(result.isEmpty());
        assertEquals(10000, result.get(0).getStart());
    }

    @Test
    public void testDelBothEndsUnbalanced() throws IOException {
        final List<VariantContext> result = delTest(3, 3, 3);
        assertTrue(result.isEmpty());
    }

    @Test
    public void testDelBothEndsUnbalancedButWithInsufficientDepthWindowCountInOneEnd() throws IOException {
        final List<VariantContext> result = delTest(3, 3, 3, UNBALANCED_MIN_DEPTH_WINDOW_COUNT - 1);
        assertTrue(result.isEmpty());
    }

    @Test
    public void testDelAlreadyBalanced() throws IOException {
        final List<VariantContext> result = delTest(3, 0.8, 3);
        assertTrue(result.isEmpty());
    }

    @NotNull
    private List<VariantContext> delTest(double startCopyNumber, double middleCopyNumber, double endCopyNumber) throws IOException {
        return delTest(startCopyNumber, middleCopyNumber, endCopyNumber, UNBALANCED_MIN_DEPTH_WINDOW_COUNT);
    }

    @NotNull
    private List<VariantContext> delTest(double startCopyNumber, double middleCopyNumber, double endCopyNumber, int startDepthWindowCount)
            throws IOException {
        StructuralVariant del = createDel();
        assertEquals(10001, del.start().cnaPosition());
        assertEquals(20000, del.end().cnaPosition());

        PurpleCopyNumber first = create(1, 10000, startCopyNumber, startDepthWindowCount);
        PurpleCopyNumber second = create(10001, 19999, middleCopyNumber, UNBALANCED_MIN_DEPTH_WINDOW_COUNT);
        PurpleCopyNumber third = create(20000, 39999, endCopyNumber, UNBALANCED_MIN_DEPTH_WINDOW_COUNT);

        RecoverStructuralVariants victim =
                new RecoverStructuralVariants(purityAdjuster, recoveredVariantFactory, Lists.newArrayList(first, second, third));

        return victim.recoverFromUnbalancedVariants(Lists.newArrayList(del), Collections.emptyList());
    }

    @NotNull
    private List<VariantContext> singleTest(int depthWindowCount, double endCopyNumber) throws IOException {
        long position = 10000;

        PurpleCopyNumber start = create(1, position, 3, depthWindowCount);
        PurpleCopyNumber end = create(position + 1, 2 * position, endCopyNumber, depthWindowCount);
        RecoverStructuralVariants victim =
                new RecoverStructuralVariants(purityAdjuster, recoveredVariantFactory, Lists.newArrayList(start, end));

        StructuralVariant single = createSingle(position);
        assertEquals(position + 1, single.start().cnaPosition());

        return victim.recoverFromUnbalancedVariants(Lists.newArrayList(single), Collections.emptyList());
    }

    @NotNull
    private static StructuralVariant createSingle(final long startPosition) {
        return PurpleDatamodelTest.createStructuralVariantSingleBreakend("1", startPosition, 0.9).build();
    }

    @NotNull
    private static StructuralVariant createDel() {
        return PurpleDatamodelTest.createStructuralVariant("1", 10000, "1", 20000, StructuralVariantType.DEL, 0.9, 0.9).build();
    }

    @NotNull
    private static PurpleCopyNumber create(final long start, final long end, final double copyNumber, int depthWindowCount) {
        return PurpleTestUtils.createCopyNumber("1", start, end, copyNumber).depthWindowCount(depthWindowCount).build();
    }
}
