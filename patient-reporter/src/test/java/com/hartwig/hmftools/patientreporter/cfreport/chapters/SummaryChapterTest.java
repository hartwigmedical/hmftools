package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Sets;

import org.junit.Test;

public class SummaryChapterTest {

    @Test
    public void canSortSummaryOfGenesCoorectly() {
        Set<String> genes = Sets.newHashSet("A", "C", "B");
        Set<String> sortedGenes = SummaryChapter.sortSummaryGenes(genes);
        Set<String> correctGenes = Sets.newHashSet("A", "B", "C");
        assertEquals(correctGenes, sortedGenes);
    }

}