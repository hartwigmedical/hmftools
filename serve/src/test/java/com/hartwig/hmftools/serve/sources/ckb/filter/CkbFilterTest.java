package com.hartwig.hmftools.serve.sources.ckb.filter;

import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.serve.sources.ckb.CkbTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CkbFilterTest {

    @Test
    public void canFilterOnKeywords() {
        CkbFilter filter = new CkbFilter(createFilterEntrySet(CkbFilterType.FILTER_ANY_VARIANT_WITH_KEYWORD, "benign"));
        CkbEntry entry = CkbTestFactory.createEntryWithVariant("filter benign me!");

        List<CkbEntry> entries = filter.run(Lists.newArrayList(entry));
        assertTrue(entries.isEmpty());

        filter.reportUnusedFilterEntries();
    }

    @Test
    public void canFilterOnGenes() {
        CkbFilter fullFilter = new CkbFilter(createFilterEntrySet(CkbFilterType.FILTER_ALL_EVIDENCE_ON_GENE, "gene"));
        CkbEntry filterEntry = CkbTestFactory.createEntryWithGene("gene");
        assertTrue(fullFilter.run(Lists.newArrayList(filterEntry)).isEmpty());
        fullFilter.reportUnusedFilterEntries();

        CkbFilter exonFilter = new CkbFilter(createFilterEntrySet(CkbFilterType.FILTER_EVIDENCE_FOR_EXONS_ON_GENE, "gene"));
        CkbEntry filterExonEntry = CkbTestFactory.createEntryWithGeneAndVariant("gene", "exon 1");
        assertTrue(exonFilter.run(Lists.newArrayList(filterExonEntry)).isEmpty());
        exonFilter.reportUnusedFilterEntries();
    }

    @Test
    public void canRemoveUnresolvableFusionLegs() {
        CkbFilter filter = new CkbFilter(createFilterEntrySet(CkbFilterType.FILTER_SECONDARY_GENE_WHEN_FUSION_LEG, "FILT"));
        CkbEntry unresolvableEntry = CkbTestFactory.createEntryWithGeneAndVariant("BRAF", "FILT-BRAF");
        assertTrue(filter.run(Lists.newArrayList(unresolvableEntry)).isEmpty());
    }

    @Test
    public void canRemoveExclusiveFusionGenes() {
        CkbFilter filter = new CkbFilter(createFilterEntrySet(CkbFilterType.ALLOW_GENE_IN_FUSIONS_EXCLUSIVELY, "gene"));
        CkbEntry exclusiveFusionEntry = CkbTestFactory.createEntryWithGeneAndVariant("gene", "gene mutant");
        assertTrue(filter.run(Lists.newArrayList(exclusiveFusionEntry)).isEmpty());
    }

    @NotNull
    private static Set<CkbFilterEntry> createFilterEntrySet(@NotNull CkbFilterType type, @NotNull String value) {
        return Sets.newHashSet(ImmutableCkbFilterEntry.builder().type(type).value(value).build());
    }
}