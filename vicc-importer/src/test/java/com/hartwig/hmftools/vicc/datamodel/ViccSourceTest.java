package com.hartwig.hmftools.vicc.datamodel;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.junit.Test;

public class ViccSourceTest {

    @Test
    public void canResolveSourcesFromDisplayString() {
        for (ViccSource source : ViccSource.values()) {
            assertEquals(source, ViccSource.fromViccKnowledgebaseString(source.display()));
        }
    }

    @Test
    public void canConvertToKnowledgebase() {
        assertEquals(Knowledgebase.VICC_CGI, ViccSource.toKnowledgebase(ViccSource.CGI));
    }

    @Test
    public void unresolvedSourceLeadsToUnknown() {
        assertEquals(ViccSource.UNKNOWN, ViccSource.fromViccKnowledgebaseString("i don't know this source"));
    }
}