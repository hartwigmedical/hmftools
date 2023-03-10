package com.hartwig.hmftools.orange.algo.linx;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import org.junit.Test;

import java.util.List;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class DNAFusionEvaluatorTest {

    @Test
    public void canDetermineIfFusionIsPresent() {
        LinxFusion fusion1 = LinxOrangeTestFactory.fusionBuilder().geneStart("start 1").geneEnd("end 1").build();
        LinxFusion fusion2 = LinxOrangeTestFactory.fusionBuilder().geneStart("start 2").geneEnd("end 1").build();

        List<LinxFusion> fusions = Lists.newArrayList(fusion1, fusion2);
        assertTrue(DNAFusionEvaluator.hasFusion(fusions, "start 2", "end 1"));
        assertFalse(DNAFusionEvaluator.hasFusion(fusions, "start 2", "end 2"));
        assertFalse(DNAFusionEvaluator.hasFusion(fusions, "start 1", "end 2"));
    }
}