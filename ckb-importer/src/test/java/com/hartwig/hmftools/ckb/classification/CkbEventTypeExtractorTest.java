package com.hartwig.hmftools.ckb.classification;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.CkbTestFactory;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.common.serve.classification.EventType;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class CkbEventTypeExtractorTest {

    @Test
    public void entriesWithMultipleVariantsAreCombined() {
        Variant variant1 = CkbTestFactory.createVariant();
        Variant variant2 = CkbTestFactory.createVariant();
        CkbEntry entry = CkbTestFactory.createEntry(Lists.newArrayList(variant1, variant2));

        assertEquals(EventType.COMBINED, CkbEventTypeExtractor.classify(entry));
    }

    @Test
    public void canClassifyCharacteristics() {
        Variant characteristic = CkbTestFactory.createVariant("-", "MSI neg", "MSI neg", Strings.EMPTY);
        assertEquals(EventType.CHARACTERISTIC, CkbEventTypeExtractor.classify(CkbTestFactory.createEntry(characteristic)));
    }
}