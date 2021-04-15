package com.hartwig.hmftools.serve.sources.ckb;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CkbUtil {

    private static final Logger LOGGER = LogManager.getLogger(CkbUtil.class);
    private static final String FIELD_DELIMITER = "\t";

    private CkbUtil() {
    }

    public static void printExtractionResults(@NotNull ExtractionResult result) {
        LOGGER.info("Analysis performed on CKB extraction result");
        LOGGER.info(" {} actionable hotspot records generated", result.actionableHotspots().size());
        LOGGER.info(" {} actionable range records generated", result.actionableRanges().size());
        LOGGER.info(" {} actionable gene records generated", result.actionableGenes().size());
        LOGGER.info(" {} actionable fusion records generated", result.actionableFusions().size());
        LOGGER.info(" {} actionable tumor characteristics records generated", result.actionableCharacteristics().size());
    }

    public static void writeEventsToTsv(@NotNull String eventTsv, @NotNull List<CkbEntry> entries) throws IOException {
        List<String> lines = Lists.newArrayList();
        String header = new StringJoiner(FIELD_DELIMITER).add("profileName").add("firstGene").add("firstVariant").add("type").toString();
        lines.add(header);

        Set<CkbEvent> events = Sets.newHashSet();
        for (CkbEntry entry : entries) {
            Variant firstVariant = entry.variants().get(0);
            events.add(new CkbEvent(entry.profileName(),
                    firstVariant.gene().geneSymbol(),
                    firstVariant.variant(),
                    entry.type().toString()));
        }

        for (CkbEvent event : events) {
            lines.add(new StringJoiner(FIELD_DELIMITER).add(event.profileName())
                    .add(event.gene())
                    .add(event.variant())
                    .add(event.type())
                    .toString());
        }
        LOGGER.info("Writing {} unique CKB events to {}", lines.size() - 1, eventTsv);
        Files.write(new File(eventTsv).toPath(), lines);
    }

    private static class CkbEvent {

        @NotNull
        private final String profileName;
        @Nullable
        private final String gene;
        @Nullable
        private final String variant;
        @NotNull
        private final String type;

        public CkbEvent(@NotNull final String profileName, @Nullable final String gene, @Nullable final String variant,
                @NotNull final String type) {
            this.profileName = profileName;
            this.gene = gene;
            this.variant = variant;
            this.type = type;
        }

        @NotNull
        public String profileName() {
            return profileName;
        }

        @Nullable
        public String gene() {
            return gene;
        }

        @Nullable
        public String variant() {
            return variant;
        }

        @NotNull
        public String type() {
            return type;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final CkbEvent ckbEvent = (CkbEvent) o;
            return profileName.equals(ckbEvent.profileName) && Objects.equals(gene, ckbEvent.gene) && Objects.equals(variant,
                    ckbEvent.variant) && type.equals(ckbEvent.type);
        }

        @Override
        public int hashCode() {
            return Objects.hash(profileName, gene, variant, type);
        }
    }
}

