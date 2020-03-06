package com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.knowledgebasegenerator.AllGenomicEvents;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.ImmutableKnownAmplificationDeletion;
import com.hartwig.hmftools.knowledgebasegenerator.cnv.KnownAmplificationDeletion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CollapseSources {

    private static final Logger LOGGER = LogManager.getLogger(CollapseSources.class);
    private static final String SOURCE_LINK_SEPARATOR = ";";

    public static List<KnownAmplificationDeletion> collapeeSourseInformation(@NotNull AllGenomicEvents allGenomicEvents) {

        Map<SourceKey, List<KnownAmplificationDeletion>> collapseSource = Maps.newHashMap();
        for (KnownAmplificationDeletion amplification : allGenomicEvents.knownAmplifications()) {
            SourceKey key = new SourceKey(amplification.gene());
            List<KnownAmplificationDeletion> items = collapseSource.get(key);
            if (items == null) {
                items = Lists.newArrayList();
            }
            items.add(amplification);
            collapseSource.put(key, items);
        }

        List<KnownAmplificationDeletion> sourceMerged = Lists.newArrayList();
        for (Map.Entry<SourceKey, List<KnownAmplificationDeletion>> entry : collapseSource.entrySet()) {
            List<KnownAmplificationDeletion> itemsForKey = entry.getValue();
            for (KnownAmplificationDeletion amps : itemsForKey) {
                String gene = amps.gene();
                String eventType = amps.eventType();
                //TODO add delimiter
                String source = amps.source() + SOURCE_LINK_SEPARATOR + amps.source();
                String sourceLink = amps.sourceLink() + SOURCE_LINK_SEPARATOR + amps.sourceLink();
                sourceMerged.add(ImmutableKnownAmplificationDeletion.builder()
                        .gene(gene)
                        .eventType(eventType)
                        .source(source)
                        .sourceLink(sourceLink)
                        .build());
            }
        }
        return sourceMerged;
    }

    private static class SourceKey {

        @NotNull
        private final String gene;

        private SourceKey(@NotNull final String gene) {
            this.gene = gene;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final SourceKey key = (SourceKey) o;
            return Objects.equals(gene, key.gene);
        }

        @Override
        public int hashCode() {
            return Objects.hash(gene);
        }
    }
}
