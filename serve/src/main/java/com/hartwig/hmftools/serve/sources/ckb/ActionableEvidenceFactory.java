package com.hartwig.hmftools.serve.sources.ckb;

import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.evidence.Evidence;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.sources.vicc.curation.DrugCurator;
import com.hartwig.hmftools.serve.sources.vicc.curation.EvidenceLevelCurator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.immutables.value.internal.$guava$.annotations.$VisibleForTesting;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class ActionableEvidenceFactory {
    private static final Logger LOGGER = LogManager.getLogger(ActionableEvidenceFactory.class);

    @NotNull
    private final DoidLookup missingDoidLookup;
    @NotNull
    private final DrugCurator drugCurator;
    @NotNull
    private final EvidenceLevelCurator evidenceLevelCurator;

    ActionableEvidenceFactory(@NotNull final DoidLookup missingDoidLookup, @NotNull final DrugCurator drugCurator,
            @NotNull final EvidenceLevelCurator evidenceLevelCurator) {
        this.missingDoidLookup = missingDoidLookup;
        this.drugCurator = drugCurator;
        this.evidenceLevelCurator = evidenceLevelCurator;
    }

    @NotNull
    public Set<ActionableEvent> toActionableEvents(@NotNull CkbEntry entry) {
        Set<ActionableEvent> actionableEvents = Sets.newHashSet();

        EvidenceLevel level = EvidenceLevel.A;
        EvidenceDirection direction = EvidenceDirection.RESPONSIVE;
        Set<String> urls = Sets.newHashSet();
        String cancerType = Strings.EMPTY;
        String therapyName = Strings.EMPTY;
        for (Evidence evidence : entry.evidences()) {
            therapyName = evidence.therapy().therapyName();
            level = resolveLevel(evidence.ampCapAscoEvidenceLevel());
            direction = resolveDirection(evidence.responseType());
            cancerType = evidence.indication().name();

            for (Reference reference : evidence.references()) {
                urls.add(reference.url());
            }

        }

        ImmutableActionableEvidence.Builder builder =
                ImmutableActionableEvidence.builder().source(Knowledgebase.CKB).level(level).direction(direction).urls(urls);
        actionableEvents.add(builder.cancerType(cancerType).doid(Strings.EMPTY).treatment(therapyName).build());

        return actionableEvents;
    }

    @Nullable
    @$VisibleForTesting
    static EvidenceLevel resolveLevel(@Nullable String evidenceLabel) {
        if (evidenceLabel == null) {
            return null;
        }

        EvidenceLevel level = EvidenceLevel.fromString(evidenceLabel);
        if (level == null) {
            LOGGER.warn("Could not resolve evidence label '{}'", evidenceLabel);
        }
        return level;
    }

    @Nullable
    @VisibleForTesting
    static EvidenceDirection resolveDirection(@Nullable String direction) {

        return null;
    }
}
