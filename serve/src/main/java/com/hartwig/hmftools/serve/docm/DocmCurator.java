package com.hartwig.hmftools.serve.docm;

import java.util.List;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.serve.hotspot.ProteinKeyFormatter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class DocmCurator {

    private static final Logger LOGGER = LogManager.getLogger(DocmCurator.class);

    private static final Set<CurationKey> ENTRY_BLACKLIST = Sets.newHashSet();

    static {
        // Not clear what the "minus" means, so ignoring. Could be DEL?
        ENTRY_BLACKLIST.add(new CurationKey("BRAF", "ENST00000288602", "K601-"));
        ENTRY_BLACKLIST.add(new CurationKey("PIK3R1", "ENST00000396611", "T576-"));
        ENTRY_BLACKLIST.add(new CurationKey("FLT3", "ENST00000380982", "D835-"));
        ENTRY_BLACKLIST.add(new CurationKey("KIT", "ENST00000288135", "V559-"));
        ENTRY_BLACKLIST.add(new CurationKey("FLT3", "ENST00000380982", "I836-"));

        // Unclear what these variants means. Maybe these are INS?
        ENTRY_BLACKLIST.add(new CurationKey("NPM1", "ENST00000296930", "W288SX"));
        ENTRY_BLACKLIST.add(new CurationKey("NPM1", "ENST00000296930", "W288PX"));
        ENTRY_BLACKLIST.add(new CurationKey("NPM1", "ENST00000296930", "W288HX"));
        ENTRY_BLACKLIST.add(new CurationKey("WT1", "ENST00000332351", "S381LYG"));
        ENTRY_BLACKLIST.add(new CurationKey("ERBB2", "ENST00000269571", "G778GSPX"));
        ENTRY_BLACKLIST.add(new CurationKey("ERBB2", "ENST00000269571", "G776CX"));

        // Variants that don't exist on the configured transcript TODO verify
        ENTRY_BLACKLIST.add(new CurationKey("BRAF", "ENST00000496384", "F203L"));
        ENTRY_BLACKLIST.add(new CurationKey("BRAF", "ENST00000496384", "T207I"));
        ENTRY_BLACKLIST.add(new CurationKey("BRAF", "ENST00000496384", "G77L"));
        ENTRY_BLACKLIST.add(new CurationKey("BRAF", "ENST00000496384", "G77V"));
    }

    @NotNull
    public static List<DocmEntry> curate(@NotNull List<DocmEntry> entries) {
        List<DocmEntry> curatedEntries = Lists.newArrayList();
        for (DocmEntry entry : entries) {
            CurationKey key = new CurationKey(entry.gene(), entry.transcript(), entry.proteinAnnotation());
            if (ENTRY_BLACKLIST.contains(key)) {
                LOGGER.debug("Removing DocmEntry '{}' because of blacklist curation.",
                        ProteinKeyFormatter.toProteinKey(entry.gene(), entry.transcript(), entry.proteinAnnotation()));
            } else {
                curatedEntries.add(entry);
            }
        }
        return curatedEntries;
    }

    private static final class CurationKey {

        @NotNull
        private final String gene;
        @NotNull
        private final String transcript;
        @NotNull
        private final String proteinAnnotation;

        public CurationKey(@NotNull final String gene, @NotNull final String transcript, @NotNull final String proteinAnnotation) {
            this.gene = gene;
            this.transcript = transcript;
            this.proteinAnnotation = proteinAnnotation;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final CurationKey that = (CurationKey) o;
            return gene.equals(that.gene) && transcript.equals(that.transcript) && proteinAnnotation.equals(that.proteinAnnotation);
        }

        @Override
        public int hashCode() {
            return Objects.hash(gene, transcript, proteinAnnotation);
        }
    }
}
