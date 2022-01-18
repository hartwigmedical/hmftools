package com.hartwig.hmftools.serve.sources.actin.curation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.sources.actin.ActinTestFactory;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.junit.Test;

public class ActinCuratorTest {

    @Test
    public void canCurateActinTrialsOnGenes() {
        String firstGene = CurationFactory.GENE_MAPPINGS.keySet().iterator().next();

        ActinEntry actinEntry = ImmutableActinEntry.builder().from(ActinTestFactory.createTestEntry()).gene(firstGene).build();

        ActinCurator curator = new ActinCurator();
        ActinEntry actinEntryCurated = curator.run(Lists.newArrayList(actinEntry)).get(0);
        assertNotNull(actinEntryCurated);

        String expectedGene = CurationFactory.GENE_MAPPINGS.get(firstGene);
        assertEquals(expectedGene, actinEntryCurated.gene());

        curator.reportUnusedCurationEntries();
    }
}