package com.hartwig.hmftools.serve.sources.vicc.curation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ViccDrugCurator {

    private static final Logger LOGGER = LogManager.getLogger(ViccDrugCurator.class);

    @NotNull
    private final Set<DrugCurationKey> evaluatedDrugCurationKeys = Sets.newHashSet();

    public ViccDrugCurator() {
    }

    @NotNull
    public List<List<String>> curate(@NotNull ViccSource source, @NotNull EvidenceLevel level, @NotNull String treatment) {
        DrugCurationKey key = new DrugCurationKey(source, level, treatment);
        evaluatedDrugCurationKeys.add(key);

        List<List<String>> curated = Lists.newArrayList();
        if (DrugCurationFactory.DRUG_MAPPINGS.containsKey(key)) {
            curated = DrugCurationFactory.DRUG_MAPPINGS.get(key).drugs();
            LOGGER.debug("Mapping VICC drug '{}' to '{}'", key, curated);
        } else if (DrugCurationFactory.DRUG_BLACKLIST.contains(key)) {
            LOGGER.debug("Blacklisting VICC drug '{}'", key);
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
                LOGGER.debug("Key '{}' hasn't been used during VICC drug blacklisting", key);
            }
        }

        for (DrugCurationKey key : DrugCurationFactory.DRUG_MAPPINGS.keySet()) {
            if (!evaluatedDrugCurationKeys.contains(key)) {
                unusedKeyCount++;
                LOGGER.debug("Key '{}' hasn't been used during VICC drug mapping", key);
            }
        }

        int totalKeyCount = DrugCurationFactory.DRUG_BLACKLIST.size() + DrugCurationFactory.DRUG_MAPPINGS.keySet().size();
        LOGGER.debug("Found {} unused VICC drug curation entries. {} keys have been requested against {} entries",
                unusedKeyCount,
                evaluatedDrugCurationKeys.size(),
                totalKeyCount);
    }
}
