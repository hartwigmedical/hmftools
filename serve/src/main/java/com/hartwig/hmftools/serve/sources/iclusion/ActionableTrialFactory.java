package com.hartwig.hmftools.serve.sources.iclusion;

import java.util.List;
import java.util.Set;

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
import org.jetbrains.annotations.NotNull;

public class ActionableTrialFactory {

    private static final Logger LOGGER = LogManager.getLogger(ActionableTrialFactory.class);

    static final String CANCER_DOID = "162";
    static final String ORGAN_SYSTEM_CANCER_DOID = "0050686";

    static final CancerType HEMATOLOGIC_CANCER_TYPE = ImmutableCancerType.builder().name("Hematologic cancer").doid("2531").build();

    @NotNull
    private final DoidLookup missingDoidLookup;

    public ActionableTrialFactory(@NotNull final DoidLookup missingDoidLookup) {
        this.missingDoidLookup = missingDoidLookup;
    }

    @NotNull
    public List<ActionableTrial> toActionableTrials(@NotNull IclusionTrial trial, @NotNull String sourceEvent) {
        ImmutableActionableTrial.Builder actionableBuilder = toActionableBuilder(trial, sourceEvent);

        Set<CancerType> blacklistCancerTypes = determineBlacklistedCancerTypes(trial);

        List<ActionableTrial> actionableTrials = Lists.newArrayList();
        for (IclusionTumorLocation tumorLocation : trial.tumorLocations()) {
            for (String doid : determineDoids(tumorLocation)) {
                String doidCurated = curateDoid(doid);
                CancerType applicable = ImmutableCancerType.builder().name(tumorLocation.primaryTumorLocation()).doid(doidCurated).build();

                Set<CancerType> finalBlacklistedCancerTypes = Sets.newHashSet(blacklistCancerTypes.iterator());
                if (doidCurated.equals(CANCER_DOID)) {
                    finalBlacklistedCancerTypes.add(HEMATOLOGIC_CANCER_TYPE);
                }

                actionableTrials.add(actionableBuilder.applicableCancerType(applicable)
                        .blacklistCancerTypes(finalBlacklistedCancerTypes)
                        .build());
            }
        }
        return actionableTrials;
    }

    @NotNull
    private static ImmutableActionableTrial.Builder toActionableBuilder(@NotNull IclusionTrial trial, @NotNull String sourceEvent) {
        return ImmutableActionableTrial.builder()
                .source(Knowledgebase.ICLUSION)
                .sourceEvent(sourceEvent)
                .sourceUrls(Sets.newHashSet("https://trial-eye.com/hmf/" + trial.id()))
                .treatment(trial.acronym())
                .level(EvidenceLevel.B)
                .direction(EvidenceDirection.RESPONSIVE)
                .evidenceUrls(Sets.newHashSet());
    }

    @NotNull
    private Set<CancerType> determineBlacklistedCancerTypes(@NotNull IclusionTrial trial) {
        Set<CancerType> blacklistedCancerTypes = Sets.newHashSet();
        for (IclusionTumorLocation blacklistTumorLocation : trial.blacklistedTumorLocations()) {
            for (String doid : determineDoids(blacklistTumorLocation)) {
                blacklistedCancerTypes.add(ImmutableCancerType.builder()
                        .name(blacklistTumorLocation.primaryTumorLocation())
                        .doid(curateDoid(doid))
                        .build());
            }
        }

        return blacklistedCancerTypes;
    }

    @NotNull
    private List<String> determineDoids(@NotNull IclusionTumorLocation tumorLocation) {
        if (!tumorLocation.doids().isEmpty()) {
            return tumorLocation.doids();
        }

        Set<String> manualDoids = missingDoidLookup.lookupDoidsForCancerType(tumorLocation.primaryTumorLocation());
        if (manualDoids == null) {
            LOGGER.warn("No doids could be derived for iClusion primary tumor location '{}'", tumorLocation.primaryTumorLocation());
            manualDoids = Sets.newHashSet();
        } else {
            LOGGER.debug("Resolved doids to '{}' for iClusion primary tumor location '{}'",
                    manualDoids,
                    tumorLocation.primaryTumorLocation());
        }

        return Lists.newArrayList(manualDoids.iterator());
    }

    @NotNull
    @VisibleForTesting
    static String curateDoid(@NotNull String doid) {
        if (doid.equals("0050686")) {
            return CANCER_DOID;
        } else if (doid.equals("MESH: D009382") || doid.equals("UNKNOWN")) {
            return ORGAN_SYSTEM_CANCER_DOID;
        } else {
            return doid;
        }
    }
}