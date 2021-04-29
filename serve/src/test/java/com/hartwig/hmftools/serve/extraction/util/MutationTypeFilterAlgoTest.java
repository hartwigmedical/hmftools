package com.hartwig.hmftools.serve.extraction.util;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.serve.DriverGeneTestFactory;

import org.junit.Test;

public class MutationTypeFilterAlgoTest {

    @Test
    public void canDetermineMutationFilter() {
        String tsg = "tsg";
        String onco = "onco";

        MutationTypeFilterAlgo algo = new MutationTypeFilterAlgo(DriverGeneTestFactory.createDriverGenes(tsg, onco));

        assertEquals(MutationTypeFilter.NONSENSE_OR_FRAMESHIFT, algo.determine(onco, "EXON 9 FRAMESHIFT"));
        assertEquals(MutationTypeFilter.SPLICE, algo.determine(onco, "Exon 12 splice site insertion"));
        assertEquals(MutationTypeFilter.SPLICE, algo.determine(onco, "EXON 14 SKIPPING MUTATION"));
        assertEquals(MutationTypeFilter.INFRAME_DELETION, algo.determine(onco, "EGFR exon 19 deletions"));
        assertEquals(MutationTypeFilter.INFRAME_DELETION, algo.determine(onco, "Null (Partial deletion of Exons 2 & 3)"));
        assertEquals(MutationTypeFilter.INFRAME_INSERTION, algo.determine(onco, "Exon 20 insertions"));
        assertEquals(MutationTypeFilter.INFRAME, algo.determine(onco, "Exon 20 insertions/deletions"));
        assertEquals(MutationTypeFilter.INFRAME, algo.determine(onco, "Exon 19 deletion/insertion"));
        assertEquals(MutationTypeFilter.MISSENSE, algo.determine(onco, "mut"));
        assertEquals(MutationTypeFilter.INFRAME_INSERTION, algo.determine(onco, "insertion"));
        assertEquals(MutationTypeFilter.ANY, algo.determine(tsg, "mut"));
        assertEquals(MutationTypeFilter.NONSENSE_OR_FRAMESHIFT, algo.determine(tsg, "frameshift"));

        assertEquals(MutationTypeFilter.ANY, algo.determine("NOT-A-GENE", "abcd"));
    }
}