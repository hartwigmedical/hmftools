package com.hartwig.hmftools.patientreporter.cfreport.chapters.analysed;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientreporter.cfreport.chapters.analysed.SummaryChapter;

import org.junit.Test;

public class SummaryChapterTest {

    @Test
    public void canSortSummaryOfGenesCorrectly() {
        Set<String> genes = Sets.newHashSet("A", "C", "B");
        Set<String> sortedGenes = SummaryChapter.sortGenes(genes);
        Set<String> correctGenes = Sets.newHashSet("A", "B", "C");
        assertEquals(correctGenes, sortedGenes);
    }
}