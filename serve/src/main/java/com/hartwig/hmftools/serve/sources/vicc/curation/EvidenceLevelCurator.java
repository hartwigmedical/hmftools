package com.hartwig.hmftools.serve.sources.vicc.curation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class EvidenceLevelCurator {

    private static final Logger LOGGER = LogManager.getLogger(EvidenceLevelCurator.class);

    @NotNull
    private final Set<EvidenceLevelCurationKey> evaluatedEvidenceLevelCurationKeys = Sets.newHashSet();

    public EvidenceLevelCurator() {
    }

    @NotNull
    public EvidenceLevel curate(@NotNull ViccSource source, @NotNull List<String> genes, @NotNull String treatment,
            @NotNull EvidenceLevel level, @NotNull EvidenceDirection direction) {
        EvidenceLevelCurationKey key = ImmutableEvidenceLevelCurationKey.builder()
                .source(source)
                .genes(genes)
                .treatment(treatment)
                .level(level)
                .direction(direction)
                .build();

        evaluatedEvidenceLevelCurationKeys.add(key);

        EvidenceLevel curated = level;
        if (EvidenceLevelCurationFactory.EVIDENCE_LEVEL_MAPPINGS.containsKey(key)) {
            curated = EvidenceLevelCurationFactory.EVIDENCE_LEVEL_MAPPINGS.get(key);
            LOGGER.debug("Mapping VICC evidence level from '{}' to '{}' for {}", level, curated, key);
        }

        return curated;
    }

    public void reportUnusedCurationKeys() {
        int unusedKeyCount = 0;

        for (EvidenceLevelCurationKey key : EvidenceLevelCurationFactory.EVIDENCE_LEVEL_MAPPINGS.keySet()) {
            if (!evaluatedEvidenceLevelCurationKeys.contains(key)) {
                unusedKeyCount++;
                LOGGER.debug("Key '{}' hasn't been used during VICC evidence level mapping", key);
            }
        }

        int totalKeyCount = EvidenceLevelCurationFactory.EVIDENCE_LEVEL_MAPPINGS.keySet().size();
        LOGGER.debug("Found {} unused VICC evidence level curation entries. {} keys have been requested against {} entries",
                unusedKeyCount,
                evaluatedEvidenceLevelCurationKeys.size(),
                totalKeyCount);
    }
}
