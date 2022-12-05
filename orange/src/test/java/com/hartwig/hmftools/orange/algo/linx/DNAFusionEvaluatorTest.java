package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxTestFactory;

import org.junit.Test;

public class DNAFusionEvaluatorTest {

    @Test
    public void canDetermineIfFusionIsPresent() {
        LinxFusion fusion1 = LinxTestFactory.fusionBuilder().geneStart("start 1").geneEnd("end 1").build();
        LinxFusion fusion2 = LinxTestFactory.fusionBuilder().geneStart("start 2").geneEnd("end 1").build();

        List<LinxFusion> fusions = Lists.newArrayList(fusion1, fusion2);
        assertTrue(DNAFusionEvaluator.hasFusion(fusions, "start 2", "end 1"));
        assertFalse(DNAFusionEvaluator.hasFusion(fusions, "start 2", "end 2"));
        assertFalse(DNAFusionEvaluator.hasFusion(fusions, "start 1", "end 2"));
    }
}