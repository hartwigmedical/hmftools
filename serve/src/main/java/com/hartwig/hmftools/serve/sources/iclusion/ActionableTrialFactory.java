package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTumorLocation;
import com.hartwig.hmftools.serve.cancertype.CancerType;
import com.hartwig.hmftools.serve.cancertype.ImmutableCancerType;
import com.hartwig.hmftools.serve.curation.DoidLookup;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ActionableTrialFactory {

    private static final String DELIMITER = ",";

    private static final Logger LOGGER = LogManager.getLogger(ActionableTrialFactory.class);

    @NotNull
    private final DoidLookup missingDoidLookup;

    public ActionableTrialFactory(@NotNull final DoidLookup missingDoidLookup) {
        this.missingDoidLookup = missingDoidLookup;
    }

    @NotNull
    public List<ActionableTrial> toActionableTrials(@NotNull IclusionTrial trial, @NotNull String sourceEvent) {
        Set<CancerType> blacklistedCancerTypes = Sets.newHashSet();

        // Source event is appended later
        ImmutableActionableTrial.Builder actionableBuilder = ImmutableActionableTrial.builder()
                .source(Knowledgebase.ICLUSION)
                .sourceEvent(sourceEvent)
                .sourceUrls(Sets.newHashSet("https://trial-eye.com/hmf/" + trial.id()))
                .treatment(trial.acronym())
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Sets.newHashSet());

        List<ActionableTrial> actionableTrials = Lists.newArrayList();
        Set<String> blacklistTumorLocations = Sets.newHashSet();
        Set<String> blacklistedDoids = Sets.newHashSet();
        for (IclusionTumorLocation blacklistTumorLocation : trial.blacklistedTumorLocations()) {
            List<String> blacklistDoids = blacklistTumorLocation.doids();
            blacklistTumorLocations.add(blacklistTumorLocation.primaryTumorLocation());
            if (blacklistDoids.isEmpty()) {
                Set<String> manualDoidsBlacklist =
                        missingDoidLookup.lookupDoidsForCancerType(blacklistTumorLocation.primaryTumorLocation());
                if (manualDoidsBlacklist == null) {
                    LOGGER.warn("No doids could be derived for iClusion primary tumor location '{}'",
                            blacklistTumorLocation.primaryTumorLocation());
                } else {
                    LOGGER.debug("Resolved doids to '{}' for iClusion primary tumor location '{}'",
                            manualDoidsBlacklist,
                            blacklistTumorLocation.primaryTumorLocation());
                    blacklistDoids = Lists.newArrayList(manualDoidsBlacklist.iterator());
                }
            }
            for (String doid : blacklistDoids) {
                String doidCorrected = extractDoid(doid);
                blacklistedDoids.add(doidCorrected);
            }
        }

        if (!blacklistTumorLocations.isEmpty() && !blacklistedDoids.isEmpty()) {
            blacklistedCancerTypes.add(ImmutableCancerType.builder()
                    .name(urlsToString(blacklistTumorLocations))
                    .doid(urlsToString(blacklistedDoids))
                    .build());
        }

        for (IclusionTumorLocation tumorLocation : trial.tumorLocations()) {
            List<String> doids = tumorLocation.doids();
            if (doids.isEmpty()) {
                Set<String> manualDoids = missingDoidLookup.lookupDoidsForCancerType(tumorLocation.primaryTumorLocation());
                if (manualDoids == null) {
                    LOGGER.warn("No doids could be derived for iClusion primary tumor location '{}'", tumorLocation.primaryTumorLocation());
                } else {
                    LOGGER.debug("Resolved doids to '{}' for iClusion primary tumor location '{}'",
                            manualDoids,
                            tumorLocation.primaryTumorLocation());
                    doids = Lists.newArrayList(manualDoids.iterator());
                }
            }
            for (String doid : doids) {
                String doidCorrected = extractDoid(doid);

                if (doidCorrected.equals("162")) {
                    blacklistedCancerTypes.add(ImmutableCancerType.builder().name("Hematologic cancer").doid("2531").build());
                }

                actionableTrials.add(actionableBuilder.applicableCancerType(ImmutableCancerType.builder()
                        .name(tumorLocation.primaryTumorLocation())
                        .doid(doidCorrected)
                        .build()).blacklistCancerTypes(blacklistedCancerTypes).build());
            }
        }
        return actionableTrials;
    }

    @NotNull
    @VisibleForTesting
    static String extractDoid(@NotNull String doid) {
        String doidCorrected = doid;
        if (doidCorrected.equals("0050686")) {
            doidCorrected = "162";
        }
        return doidCorrected;
    }

    @NotNull
    private static String urlsToString(@NotNull Set<String> urls) {
        StringJoiner joiner = new StringJoiner(DELIMITER);
        for (String url : urls) {
            joiner.add(url);
        }
        return joiner.toString();
    }
}