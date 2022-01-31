package com.hartwig.hmftools.serve.sources.actin.filter;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.sources.actin.ActinTestFactory;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ActinFilterTest {

    @Test
    public void canFilter() {
        List<ActinFilterEntry> filterEntries = Lists.newArrayList();
        filterEntries.add(create(ActinFilterType.FILTER_VARIANT_ON_GENE, "BRAF V600E"));
        filterEntries.add(create(ActinFilterType.FILTER_EVERYTHING_FOR_GENE, "KRAS"));
        filterEntries.add(create(ActinFilterType.FILTER_EVERYTHING_FOR_RULE, ActinRule.WILDTYPE_OF_GENE_X.toString()));

        ActinFilter filter = new ActinFilter(filterEntries);

        ActinEntry brafV600E =
                ImmutableActinEntry.builder().from(ActinTestFactory.createTestEntry()).gene("BRAF").mutation("V600E").build();
        assertTrue(filter.run(Lists.newArrayList(brafV600E)).isEmpty());

        ActinEntry kras = ImmutableActinEntry.builder().from(ActinTestFactory.createTestEntry()).gene("KRAS").build();
        assertTrue(filter.run(Lists.newArrayList(kras)).isEmpty());

        ActinEntry wildtype =
                ImmutableActinEntry.builder().from(ActinTestFactory.createTestEntry()).rule(ActinRule.WILDTYPE_OF_GENE_X).build();
        assertTrue(filter.run(Lists.newArrayList(wildtype)).isEmpty());

        ActinEntry tmb = ImmutableActinEntry.builder().from(ActinTestFactory.createTestEntry()).rule(ActinRule.TMB_OF_AT_LEAST_X).build();
        assertFalse(filter.run(Lists.newArrayList(tmb)).isEmpty());

        filter.reportUnusedFilterEntries();
    }

    @NotNull
    private static ActinFilterEntry create(@NotNull ActinFilterType type, @NotNull String value) {
        return ImmutableActinFilterEntry.builder().type(type).value(value).build();
    }
}
