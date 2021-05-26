package com.hartwig.hmftools.serve.sources.ckb.filter;

import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.serve.sources.ckb.CkbTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CkbFilterTest {

    @Test
    public void canFilterOnKeywords() {
        CkbFilter filter = new CkbFilter(createFilterEntryList(CkbFilterType.FILTER_EVENT_WITH_KEYWORD, "benign"));
        CkbEntry entry = CkbTestFactory.createEntryWithVariant("filter benign me!");
        assertTrue(filter.run(Lists.newArrayList(entry)).isEmpty());

        filter.reportUnusedFilterEntries();
    }

    @Test
    public void canFilterOnGenes() {
        CkbFilter fullFilter = new CkbFilter(createFilterEntryList(CkbFilterType.FILTER_ALL_EVIDENCE_ON_GENE, "gene"));
        CkbEntry filterEntry = CkbTestFactory.createEntryWithGene("gene");
        assertTrue(fullFilter.run(Lists.newArrayList(filterEntry)).isEmpty());
        fullFilter.reportUnusedFilterEntries();

        CkbFilter exonFilter = new CkbFilter(createFilterEntryList(CkbFilterType.FILTER_EVIDENCE_FOR_EXONS_ON_GENE, "gene"));
        CkbEntry filterExonEntry = CkbTestFactory.createEntryWithGeneAndVariant("gene", "exon 1");
        assertTrue(exonFilter.run(Lists.newArrayList(filterExonEntry)).isEmpty());
        exonFilter.reportUnusedFilterEntries();
    }

    @Test
    public void canRemoveUnresolvableFusionLegs() {
        CkbFilter filter = new CkbFilter(createFilterEntryList(CkbFilterType.FILTER_SECONDARY_GENE_WHEN_FUSION_LEG, "FILT"));
        CkbEntry unresolvableEntry = CkbTestFactory.createEntryWithGeneAndVariant("BRAF", "FILT-BRAF");
        assertTrue(filter.run(Lists.newArrayList(unresolvableEntry)).isEmpty());
    }

    @Test
    public void canRemoveExclusiveFusionGenes() {
        CkbFilter filter = new CkbFilter(createFilterEntryList(CkbFilterType.ALLOW_GENE_IN_FUSIONS_EXCLUSIVELY, "gene"));
        CkbEntry exclusiveFusionEntry = CkbTestFactory.createEntryWithGeneAndVariant("gene", "gene mutant");
        assertTrue(filter.run(Lists.newArrayList(exclusiveFusionEntry)).isEmpty());
    }

    @Test
    public void canFilterOnFullNames() {
        CkbFilter filter = new CkbFilter(createFilterEntryList(CkbFilterType.FILTER_EXACT_VARIANT_FULLNAME, "BRAF V600E"));
        CkbEntry entry = CkbTestFactory.createEntryWithFullName("BRAF V600E");
        assertTrue(filter.run(Lists.newArrayList(entry)).isEmpty());
    }

    @NotNull
    private static List<CkbFilterEntry> createFilterEntryList(@NotNull CkbFilterType type, @NotNull String value) {
        return Lists.newArrayList(ImmutableCkbFilterEntry.builder().type(type).value(value).build());
    }
}