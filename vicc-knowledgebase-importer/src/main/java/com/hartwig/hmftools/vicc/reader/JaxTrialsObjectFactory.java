package com.hartwig.hmftools.vicc.reader;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.ImmutableJaxTrials;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.ImmutableJaxTrialsIndications;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.ImmutableJaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.ImmutableJaxTrialsTherapies;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.ImmutableJaxTrialsVariantRequirementDetails;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsIndications;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsTherapies;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsVariantRequirementDetails;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class JaxTrialsObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(JaxTrialsObjectFactory.class);

    private static final List<Integer> EXPECTED_JAX_TRIALS_ELEMENT_SIZES = Lists.newArrayList(11);
    private static final List<Integer> EXPECTED_JAX_TRIALS_INDICATIONS_ELEMENT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_JAX_TRIALS_VARIANTREQUIREMENTDETAILS_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_JAX_TRIALS_MOLECULAIRPROFILE_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_JAX_TRIALS_THERAPIES_ELEMENT_SIZES = Lists.newArrayList(2);

    private JaxTrialsObjectFactory() {
    }

    @NotNull
    static JaxTrials create(@NotNull JsonObject objectJaxTrials) {
        Set<String> keysJaxTrials = objectJaxTrials.keySet();
        if (!EXPECTED_JAX_TRIALS_ELEMENT_SIZES.contains(keysJaxTrials.size())) {
            LOGGER.warn("Found {} in jax trials rather than the expected {}", keysJaxTrials.size(), EXPECTED_JAX_TRIALS_ELEMENT_SIZES);
            LOGGER.warn(keysJaxTrials);
        }

        return ImmutableJaxTrials.builder()
                .indications(createJaxTrialsIndications(objectJaxTrials.getAsJsonArray("indications")))
                .title(objectJaxTrials.getAsJsonPrimitive("title").getAsString())
                .gender(objectJaxTrials.get("gender").isJsonNull() ? null : objectJaxTrials.getAsJsonPrimitive("gender").getAsString())
                .nctId(objectJaxTrials.getAsJsonPrimitive("nctId").getAsString())
                .sponsors(objectJaxTrials.getAsJsonPrimitive("sponsors").getAsString())
                .recruitment(objectJaxTrials.getAsJsonPrimitive("recruitment").getAsString())
                .variantRequirements(objectJaxTrials.getAsJsonPrimitive("variantRequirements").getAsString())
                .updateDate(objectJaxTrials.getAsJsonPrimitive("updateDate").getAsString())
                .phase(objectJaxTrials.getAsJsonPrimitive("phase").getAsString())
                .variantRequirementDetails(createJaxTrialsVariantRequirementsDetails(objectJaxTrials.getAsJsonArray(
                        "variantRequirementDetails")))
                .therapies(createJaxTrialsTherapies(objectJaxTrials.getAsJsonArray("therapies")))
                .build();
    }

    @NotNull
    private static List<JaxTrialsIndications> createJaxTrialsIndications(@NotNull JsonArray arrayIndications) {
        List<JaxTrialsIndications> indicationsList = Lists.newArrayList();

        for (JsonElement indications : arrayIndications) {
            Set<String> keysIndications = indications.getAsJsonObject().keySet();
            if (!EXPECTED_JAX_TRIALS_INDICATIONS_ELEMENT_SIZES.contains(keysIndications.size())) {
                LOGGER.warn("Found {} in jax trials indications rather than the expected {}",
                        keysIndications.size(),
                        EXPECTED_JAX_TRIALS_INDICATIONS_ELEMENT_SIZES);
                LOGGER.warn(keysIndications);
            }

            indicationsList.add(ImmutableJaxTrialsIndications.builder()
                    .source(indications.getAsJsonObject().getAsJsonPrimitive("source").getAsString())
                    .id(indications.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(indications.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return indicationsList;
    }

    @NotNull
    private static List<JaxTrialsVariantRequirementDetails> createJaxTrialsVariantRequirementsDetails(
            @NotNull JsonArray arrayVariantRequirementDetails) {
        List<JaxTrialsVariantRequirementDetails> variantRequirementDetailsList = Lists.newArrayList();
        for (JsonElement variantRequirementDetails : arrayVariantRequirementDetails) {
            Set<String> keysRequirementDetails = variantRequirementDetails.getAsJsonObject().keySet();
            if (!EXPECTED_JAX_TRIALS_VARIANTREQUIREMENTDETAILS_ELEMENT_SIZES.contains(keysRequirementDetails.size())) {
                LOGGER.warn("Found {} in jax trials variant requirement details rather than the expected {}",
                        keysRequirementDetails.size(),
                        EXPECTED_JAX_TRIALS_VARIANTREQUIREMENTDETAILS_ELEMENT_SIZES);
                LOGGER.warn(keysRequirementDetails);
            }

            variantRequirementDetailsList.add(ImmutableJaxTrialsVariantRequirementDetails.builder()
                    .molecularProfiles(createJaxTrialsMolecularProfile(variantRequirementDetails.getAsJsonObject()
                            .getAsJsonObject("molecularProfile")))
                    .requirementType(variantRequirementDetails.getAsJsonObject().getAsJsonPrimitive("requirementType").getAsString())
                    .build());
        }
        return variantRequirementDetailsList;
    }

    @NotNull
    private static List<JaxTrialsMolecularProfile> createJaxTrialsMolecularProfile(@NotNull JsonObject objectMolecularProfile) {
        List<JaxTrialsMolecularProfile> molecularProfileList = Lists.newArrayList();

        Set<String> keysMolecularProfile = objectMolecularProfile.keySet();
        if (!EXPECTED_JAX_TRIALS_MOLECULAIRPROFILE_ELEMENT_SIZES.contains(keysMolecularProfile.size())) {
            LOGGER.warn("Found {}  in jax trials molecular profile rather than the expected {}",
                    keysMolecularProfile.size(),
                    EXPECTED_JAX_TRIALS_MOLECULAIRPROFILE_ELEMENT_SIZES);
            LOGGER.warn(keysMolecularProfile);
        }

        molecularProfileList.add(ImmutableJaxTrialsMolecularProfile.builder()
                .profileName(objectMolecularProfile.getAsJsonPrimitive("profileName").getAsString())
                .id(objectMolecularProfile.getAsJsonPrimitive("id").getAsString())
                .build());
        return molecularProfileList;
    }

    @NotNull
    private static List<JaxTrialsTherapies> createJaxTrialsTherapies(@NotNull JsonArray arrayTherapies) {
        List<JaxTrialsTherapies> jaxTrialsTherapiesList = Lists.newArrayList();
        for (JsonElement therapies : arrayTherapies) {
            Set<String> keysTherapies = therapies.getAsJsonObject().keySet();
            if (!EXPECTED_JAX_TRIALS_THERAPIES_ELEMENT_SIZES.contains(keysTherapies.size())) {
                LOGGER.warn("Found {} in jax trials therapies rather than the expected {}",
                        keysTherapies.size(),
                        EXPECTED_JAX_TRIALS_THERAPIES_ELEMENT_SIZES);
                LOGGER.warn(keysTherapies);
            }

            jaxTrialsTherapiesList.add(ImmutableJaxTrialsTherapies.builder()
                    .id(therapies.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .therapyName(therapies.getAsJsonObject().getAsJsonPrimitive("therapyName").getAsString())
                    .build());
        }
        return jaxTrialsTherapiesList;
    }
}
