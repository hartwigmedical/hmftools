package com.hartwig.hmftools.serve.sources.ckb.filter;

import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.classification.CkbConstants;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.serve.sources.ckb.CkbTestFactory;

import org.junit.Test;

public class CkbFilterTest {

    @Test
    public void canFilterOnKeywords() {
        CkbFilter filter = new CkbFilter();

        String firstFilterKeyword = FilterFactory.VARIANT_KEYWORDS_TO_FILTER.iterator().next();
        CkbEntry entry = CkbTestFactory.createEntryWithVariant(firstFilterKeyword + " filter me!");

        List<CkbEntry> entries = filter.run(Lists.newArrayList(entry));
        assertTrue(entries.isEmpty());

        filter.reportUnusedFilterEntries();
    }

    @Test
    public void canFilterOnGenes() {
        CkbFilter filter = new CkbFilter();

        String firstFilterGene = FilterFactory.GENES_FOR_WHICH_TO_FILTER_ALL.iterator().next();
        CkbEntry filterEntry = CkbTestFactory.createEntryWithGene(firstFilterGene);
        assertTrue(filter.run(Lists.newArrayList(filterEntry)).isEmpty());

        String firstExonFilterGene = FilterFactory.GENES_FOR_WHICH_TO_FILTER_EXON_EVENTS.iterator().next();
        CkbEntry filterExonEntry = CkbTestFactory.createEntryWithGeneAndVariant(firstExonFilterGene, "exon 1");
        assertTrue(filter.run(Lists.newArrayList(filterExonEntry)).isEmpty());

        filter.reportUnusedFilterEntries();
    }

    @Test
    public void canRemoveUnmappableGenes() {
        CkbFilter filter = new CkbFilter();

        String firstUnmappableGene = CkbConstants.UNMAPPABLE_GENES.iterator().next();
        CkbEntry unmappableEntry = CkbTestFactory.createEntryWithGene(firstUnmappableGene);
        assertTrue(filter.run(Lists.newArrayList(unmappableEntry)).isEmpty());
    }

    @Test
    public void canRemoveUnresolvableFusionLegs() {
        CkbFilter filter = new CkbFilter();

        String firstUnresolvableLeg = CkbConstants.UNRESOLVABLE_FUSION_LEGS.iterator().next();
        CkbEntry unresolvableEntry = CkbTestFactory.createEntryWithGeneAndVariant("BRAF", firstUnresolvableLeg + "-BRAF");
        assertTrue(filter.run(Lists.newArrayList(unresolvableEntry)).isEmpty());

        String firstUnmappableLeg = CkbConstants.UNMAPPABLE_GENES.iterator().next();
        CkbEntry unmappableEntry = CkbTestFactory.createEntryWithGeneAndVariant("BRAF", firstUnmappableLeg + "-BRAF");
        assertTrue(filter.run(Lists.newArrayList(unmappableEntry)).isEmpty());
    }

    @Test
    public void canRemoveExclusiveFusionGenes() {
        CkbFilter filter = new CkbFilter();

        String firstExclusiveFusionGene = CkbConstants.EXCLUSIVE_FUSION_GENES.iterator().next();
        CkbEntry exclusiveFusionEntry =
                CkbTestFactory.createEntryWithGeneAndVariant(firstExclusiveFusionGene, firstExclusiveFusionGene + " mutant");
        assertTrue(filter.run(Lists.newArrayList(exclusiveFusionEntry)).isEmpty());
    }
}