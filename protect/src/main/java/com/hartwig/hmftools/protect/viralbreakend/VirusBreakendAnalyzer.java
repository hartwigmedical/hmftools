package com.hartwig.hmftools.protect.viralbreakend;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class VirusBreakendAnalyzer {

    private static final Set<String> EXCLUDED_VIRAL_INSERTIONS =
            Sets.newHashSet("Human immunodeficiency virus", "Human immunodeficiency virus 1", "Human immunodeficiency virus 2");

    private VirusBreakendAnalyzer() {
    }

    @NotNull
    public static List<VirusBreakend> analyzeVirusBreakends(@NotNull List<VirusBreakend> virusBreakends) {
        Map<VirusBreakendAnalyzer.VirusKey, List<VirusBreakend>> itemsPerKey = Maps.newHashMap();
        for (VirusBreakend virusBreakend : virusBreakends) {
            VirusBreakendAnalyzer.VirusKey key = new VirusBreakendAnalyzer.VirusKey(virusBreakend.Reference());
            List<VirusBreakend> items = itemsPerKey.get(key);

            if (items == null) {
                items = Lists.newArrayList();
            }
            items.add(virusBreakend);
            itemsPerKey.put(key, items);
        }

        List<VirusBreakend> virusBreakendFiltered = Lists.newArrayList();
        for (Map.Entry<VirusBreakendAnalyzer.VirusKey, List<VirusBreakend>> entry : itemsPerKey.entrySet()) {
            List<VirusBreakend> itemsForKey = entry.getValue();

            int count = itemsForKey.size();
            assert count > 0;
            String virusName = itemsForKey.get(0).nameAssigned();
            if (!EXCLUDED_VIRAL_INSERTIONS.contains(virusName)) {
                virusBreakendFiltered.add(ImmutableVirusBreakend.builder()
                        .from(itemsForKey.get(0))
                        .integrations(Integer.toString(count)) // TODO: can we use integrations of file?
                        .build());
            }
        }

        return virusBreakendFiltered;
    }

    private static class VirusKey {

        @NotNull
        private final String virusId;

        private VirusKey(@NotNull final String virusId) {
            this.virusId = virusId;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final VirusBreakendAnalyzer.VirusKey key = (VirusBreakendAnalyzer.VirusKey) o;
            return Objects.equals(virusId, key.virusId);
        }

        @Override
        public int hashCode() {
            return Objects.hash(virusId);
        }
    }
}
