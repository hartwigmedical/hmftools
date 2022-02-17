package com.hartwig.hmftools.serve.sources.ckb.curation;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.serve.sources.ckb.CkbTestFactory;

import org.junit.Test;

public class CkbCuratorTest {

    @Test
    public void canCurateVariants() {
        CkbCurator curator = new CkbCurator();

        CurationEntry firstCurationKey = CurationFactory.VARIANT_MAPPINGS.keySet().iterator().next();
        CkbEntry entry = CkbTestFactory.createEntryWithGeneAndVariant(firstCurationKey.geneSymbol(), firstCurationKey.variant());

        List<CkbEntry> entries = curator.run(Lists.newArrayList(entry));

        Variant firstVariant = entries.get(0).variants().get(0);

        CurationEntry firstCuratedValue = CurationFactory.VARIANT_MAPPINGS.get(firstCurationKey);
        assertEquals(firstCuratedValue.geneSymbol(), firstVariant.gene().geneSymbol());
        assertEquals(firstCuratedValue.variant(), firstVariant.variant());

        curator.reportUnusedCurationEntries();
    }
}