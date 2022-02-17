package com.hartwig.hmftools.serve.transvar;

import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.transvar.datamodel.ImmutableTransvarRecord;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarAnnotation;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarRecord;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TransvarTest {

    @Test
    public void noHotspotsWhenTransvarProcessReturnsEmpty() {
        Transvar transvar = returnsNoTransvarRecord();
        assertTrue(transvar.resolve("BRAF", "ENST00000288602", "V600E").isEmpty());
    }

    @Test
    public void noHotspotsWhenGeneIsUnknown() {
        Transvar transvar = returnsSingleTransvarRecord(createTestRecord());
        assertTrue(transvar.resolve("DoesNotExist", null, Strings.EMPTY).isEmpty());
    }

    @Test
    public void noHotspotsWhenRecordIsNotOnSpecificTranscript() {
        Transvar transvar = returnsSingleTransvarRecord(createTestRecord());
        assertTrue(transvar.resolve("BRAF", "DoesNotExist", Strings.EMPTY).isEmpty());
    }

    @NotNull
    private static TransvarRecord createTestRecord() {
        return ImmutableTransvarRecord.builder()
                .chromosome("7")
                .gdnaPosition(10)
                .transcript("X")
                .variantSpanMultipleExons(false)
                .annotation(new TransvarAnnotation() {
                })
                .build();
    }

    @NotNull
    private static Transvar returnsSingleTransvarRecord(@NotNull TransvarRecord record) {
        return TransvarTestFactory.testTransvar37((gene, proteinAnnotation) -> Lists.newArrayList(record));
    }

    @NotNull
    private static Transvar returnsNoTransvarRecord() {
        return TransvarTestFactory.testTransvar37((gene, proteinAnnotation) -> Lists.newArrayList());
    }
}