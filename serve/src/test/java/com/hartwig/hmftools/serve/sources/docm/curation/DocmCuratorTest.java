package com.hartwig.hmftools.serve.sources.docm.curation;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.sources.docm.DocmEntry;
import com.hartwig.hmftools.serve.sources.docm.ImmutableDocmEntry;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class DocmCuratorTest {

    @Test
    public void canCurateDocmEntries() {
        DocmEntry validEntry = ImmutableDocmEntry.builder().gene("gene").transcript("transcript").proteinAnnotation("annotation").build();

        CurationKey firstBlacklistKey = firstBlacklistKey();
        DocmEntry blacklistEntry = ImmutableDocmEntry.builder()
                .gene(firstBlacklistKey.gene())
                .transcript(firstBlacklistKey.transcript())
                .proteinAnnotation(firstBlacklistKey.proteinAnnotation())
                .build();

        List<DocmEntry> entries = Lists.newArrayList(validEntry, blacklistEntry);

        DocmCurator curator = new DocmCurator();
        List<DocmEntry> curated = curator.curate(entries);

        assertEquals(1, curated.size());
        assertEquals(validEntry, curated.get(0));

        curator.reportUnusedBlacklistEntries();
    }

    @NotNull
    private static CurationKey firstBlacklistKey() {
        return CurationFactory.ENTRY_BLACKLIST.iterator().next();
    }

}