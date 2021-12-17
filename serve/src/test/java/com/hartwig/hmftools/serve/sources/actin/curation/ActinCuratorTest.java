package com.hartwig.hmftools.serve.sources.actin.curation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;
import com.hartwig.hmftools.serve.sources.actin.reader.ImmutableActinEntry;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class ActinCuratorTest {

    @Test
    public void canCurateActinTrialsOnGenes() {
        String firstGene = CurationFactory.GENE_MAPPINGS.keySet().iterator().next();

        ActinEntry actinEntry = ImmutableActinEntry.builder()
                .trial("trial1")
                .rule(ActinRule.ACTIVATION_OF_GENE_X)
                .parameters(Lists.newArrayList("RABAC1"))
                .build();


        ActinCurator curator = new ActinCurator();
        ActinEntry actinEntryCurated = curator.run(Lists.newArrayList(actinEntry)).get(0);
        assertNotNull(actinEntryCurated);

        String expectedGene = CurationFactory.GENE_MAPPINGS.get(firstGene);

        List<String> parameters = actinEntry.parameters();
        String variant = Strings.EMPTY;
        String gene = Strings.EMPTY;

        if (parameters.size() == 2) {
            gene = parameters.get(0);
            variant = parameters.get(1);
        }

        if (parameters.size() == 1) {
            gene = parameters.get(0);
        }

        assertEquals(expectedGene, gene);

        curator.reportUnusedCurationEntries();
    }
}