package com.hartwig.hmftools.serve.sources.ckb.curation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class DrugCurator {

    private static final Logger LOGGER = LogManager.getLogger(DrugCurator.class);

    @NotNull
    private final Set<DrugCurationKey> evaluatedDrugCurationKeys = Sets.newHashSet();

    public DrugCurator() {
    }

    @NotNull
    public List<List<String>> curate(@NotNull Knowledgebase source, @NotNull EvidenceLevel level, @NotNull String treatment) {
        DrugCurationKey
                key = new DrugCurationKey(source, level, treatment);
        evaluatedDrugCurationKeys.add(key);

        List<List<String>> curated = Lists.newArrayList();
        if (DrugCurationFactory.DRUG_MAPPINGS.containsKey(key)) {
            curated = DrugCurationFactory.DRUG_MAPPINGS.get(key).drugs();
            LOGGER.debug("Mapping CKB drug '{}' to '{}'", key, curated);
        } else if (DrugCurationFactory.DRUG_BLACKLIST.contains(key)) {
            LOGGER.debug("Blacklisting CKB drug '{}'", key);
        } else {
            curated.add(Lists.newArrayList(treatment.split(",")));
        }

        return curated;
    }

    public void reportUnusedCurationKeys() {
        int unusedKeyCount = 0;
        for (DrugCurationKey key : DrugCurationFactory.DRUG_BLACKLIST) {
            if (!evaluatedDrugCurationKeys.contains(key)) {
                unusedKeyCount++;
                LOGGER.debug("Key '{}' hasn't been used during CKB drug blacklisting", key);
            }
        }

        for (DrugCurationKey key : DrugCurationFactory.DRUG_MAPPINGS.keySet()) {
            if (!evaluatedDrugCurationKeys.contains(key)) {
                unusedKeyCount++;
                LOGGER.debug("Key '{}' hasn't been used during CKB drug mapping", key);
            }
        }

        int totalKeyCount = DrugCurationFactory.DRUG_BLACKLIST.size() + DrugCurationFactory.DRUG_MAPPINGS.keySet().size();
        LOGGER.debug("Found {} unused CKB drug curation entries. {} keys have been requested against {} entries",
                unusedKeyCount,
                evaluatedDrugCurationKeys.size(),
                totalKeyCount);
    }
}
