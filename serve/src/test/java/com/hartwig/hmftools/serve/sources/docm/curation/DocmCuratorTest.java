package com.hartwig.hmftools.serve.sources.docm.curation;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.sources.docm.DocmEntry;
import com.hartwig.hmftools.serve.sources.docm.ImmutableDocmEntry;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class DocmCuratorTest {

    @Test
    public void canRemoveBlacklistedDocmEntries() {
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

    @Test
    public void canMapGenes() {
        String firstMappableGene = firstGeneForMapping();
        DocmEntry entry =
                ImmutableDocmEntry.builder().gene(firstMappableGene).transcript(Strings.EMPTY).proteinAnnotation(Strings.EMPTY).build();

        DocmCurator curator = new DocmCurator();
        List<DocmEntry> curated = curator.curate(Lists.newArrayList(entry));

        assertEquals(1, curated.size());
        assertEquals(CurationFactory.GENE_MAPPINGS.get(firstMappableGene), curated.get(0).gene());
    }

    @NotNull
    private static String firstGeneForMapping() {
        return CurationFactory.GENE_MAPPINGS.keySet().iterator().next();
    }

    @NotNull
    private static CurationKey firstBlacklistKey() {
        return CurationFactory.ENTRY_BLACKLIST.iterator().next();
    }
}