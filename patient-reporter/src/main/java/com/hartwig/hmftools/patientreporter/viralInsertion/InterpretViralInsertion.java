package com.hartwig.hmftools.patientreporter.viralInsertion;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertFile;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class InterpretViralInsertion {
    private static final Logger LOGGER = LogManager.getLogger(InterpretViralInsertion.class);

    private InterpretViralInsertion() {

    }

    @NotNull
    public static List<ViralInsertion> interpretVirals(@NotNull String viralInsertTsv) throws IOException {
        List<LinxViralInsertFile> viralInsertFileList = LinxViralInsertFile.read(viralInsertTsv);
        LOGGER.info("Loaded {} viral insertions from {}", viralInsertFileList.size(), viralInsertTsv);

        Map<InterpretViralInsertion.VirusKey, List<LinxViralInsertFile>> itemsPerKey = Maps.newHashMap();
        for (LinxViralInsertFile viralInsertion : viralInsertFileList) {
            InterpretViralInsertion.VirusKey key = new InterpretViralInsertion.VirusKey(viralInsertion.VirusId);
            List<LinxViralInsertFile> items = itemsPerKey.get(key);

            if (items == null) {
                items = Lists.newArrayList();
            }
            items.add(viralInsertion);
            itemsPerKey.put(key, items);
        }

        List<ViralInsertion> viralInsertions = Lists.newArrayList();
        String virusName = Strings.EMPTY;
        int count = 0;
        for (Map.Entry<InterpretViralInsertion.VirusKey, List<LinxViralInsertFile>> entry : itemsPerKey.entrySet()) {
            List<LinxViralInsertFile> itemsForKey = entry.getValue();
            for (LinxViralInsertFile virus : itemsForKey) {
                virusName = virus.VirusName;
            }

            count = itemsForKey.size();

            viralInsertions.add(ImmutableViralInsertion.builder().virus(virusName).countVirus(count).build());
        }
        return viralInsertions;
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
            final InterpretViralInsertion.VirusKey key = (InterpretViralInsertion.VirusKey) o;
            return Objects.equals(virusId, key.virusId);
        }

        @Override
        public int hashCode() {
            return Objects.hash(virusId);
        }

    }
}
