package com.hartwig.hmftools.vicc.reader;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.jax.ImmutableJax;
import com.hartwig.hmftools.vicc.datamodel.jax.ImmutableJaxIndications;
import com.hartwig.hmftools.vicc.datamodel.jax.ImmutableJaxMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.jax.ImmutableJaxReferences;
import com.hartwig.hmftools.vicc.datamodel.jax.ImmutableJaxTherapy;
import com.hartwig.hmftools.vicc.datamodel.jax.Jax;
import com.hartwig.hmftools.vicc.datamodel.jax.JaxIndications;
import com.hartwig.hmftools.vicc.datamodel.jax.JaxMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.jax.JaxReferences;
import com.hartwig.hmftools.vicc.datamodel.jax.JaxTherapy;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class JaxObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(JaxObjectFactory.class);

    private static final List<Integer> EXPECTED_JAX_ELEMENT_SIZES = Lists.newArrayList(9);
    private static final List<Integer> EXPECTED_JAX_THERAPY_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_JAX_INDICATIONS_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_JAX_REFERENCES_ELEMENT_SIZES = Lists.newArrayList(4);
    private static final List<Integer> EXPECTED_JAX_MOLECULAR_PROFILE_ELEMENT_SIZES = Lists.newArrayList(2);

    private JaxObjectFactory() {
    }

    @NotNull
    static Jax create(@NotNull JsonObject objectJax) {
        Set<String> keysJax = objectJax.keySet();
        if (!EXPECTED_JAX_ELEMENT_SIZES.contains(keysJax.size())) {
            LOGGER.warn("Found {}  in jax rather than the expected {}", keysJax.size(), EXPECTED_JAX_ELEMENT_SIZES);
            LOGGER.warn(keysJax);
        }

        return ImmutableJax.builder()
                .responseType(objectJax.getAsJsonPrimitive("responseType").getAsString())
                .approvalStatus(objectJax.getAsJsonPrimitive("approvalStatus").getAsString())
                .molecularProfile(createMolecularProfile(objectJax.getAsJsonObject("molecularProfile")))
                .therapy(createJaxTherapy(objectJax.getAsJsonObject("therapy")))
                .evidenceType(objectJax.getAsJsonPrimitive("evidenceType").getAsString())
                .indications(createJaxIndications(objectJax.getAsJsonObject("indication")))
                .efficacyEvidence(objectJax.getAsJsonPrimitive("efficacyEvidence").getAsString())
                .references(objectJax.get("references").isJsonNull() ? null : createJaxReferences(objectJax.getAsJsonArray("references")))
                .id(objectJax.getAsJsonPrimitive("id").getAsString())
                .build();
    }

    @NotNull
    private static JaxMolecularProfile createMolecularProfile(@NotNull JsonObject objectMolecularProfile) {
        Set<String> keysMolecularProfile = objectMolecularProfile.keySet();
        if (!EXPECTED_JAX_MOLECULAR_PROFILE_ELEMENT_SIZES.contains(keysMolecularProfile.size())) {
            LOGGER.warn("Found {} in jax molecular profile rather than the expected {}",
                    keysMolecularProfile.size(),
                    EXPECTED_JAX_MOLECULAR_PROFILE_ELEMENT_SIZES);
            LOGGER.warn(keysMolecularProfile);
        }

        return ImmutableJaxMolecularProfile.builder()
                .profileName(objectMolecularProfile.getAsJsonPrimitive("profileName").getAsString())
                .id(objectMolecularProfile.getAsJsonPrimitive("id").getAsString())
                .build();
    }

    @NotNull
    private static JaxTherapy createJaxTherapy(@NotNull JsonObject objectTherapy) {
        Set<String> keysTherapy = objectTherapy.keySet();
        if (!EXPECTED_JAX_THERAPY_ELEMENT_SIZES.contains(keysTherapy.size())) {
            LOGGER.warn("Found {} in jax therapy rather than the expected {}", keysTherapy.size(), EXPECTED_JAX_THERAPY_ELEMENT_SIZES);
            LOGGER.warn(keysTherapy);
        }

        return ImmutableJaxTherapy.builder()
                .id(objectTherapy.getAsJsonPrimitive("id").getAsString())
                .therapyName(objectTherapy.getAsJsonPrimitive("therapyName").getAsString())
                .build();
    }

    @NotNull
    private static JaxIndications createJaxIndications(@NotNull JsonObject objectIndications) {
        Set<String> keysIndications = objectIndications.keySet();
        if (!EXPECTED_JAX_INDICATIONS_SIZES.contains(keysIndications.size())) {
            LOGGER.warn("Found {} in jax indications rather than the expected {}", keysIndications.size(), EXPECTED_JAX_INDICATIONS_SIZES);
            LOGGER.warn(keysIndications);
        }

        return ImmutableJaxIndications.builder()
                .source(objectIndications.getAsJsonPrimitive("source").getAsString())
                .id(objectIndications.getAsJsonPrimitive("id").getAsString())
                .name(objectIndications.getAsJsonPrimitive("name").getAsString())
                .build();
    }

    @NotNull
    private static List<JaxReferences> createJaxReferences(@NotNull JsonArray objectReferences) {
        List<JaxReferences> listReferences = Lists.newArrayList();
        for (JsonElement references : objectReferences) {
            Set<String> keysReferences = references.getAsJsonObject().keySet();
            if (!EXPECTED_JAX_REFERENCES_ELEMENT_SIZES.contains(keysReferences.size())) {
                LOGGER.warn("Found {} in jax references rather than the expected {}",
                        keysReferences.size(),
                        EXPECTED_JAX_REFERENCES_ELEMENT_SIZES);
                LOGGER.warn(keysReferences);
            }

            listReferences.add(ImmutableJaxReferences.builder()
                    .url(references.getAsJsonObject().getAsJsonPrimitive("url").getAsString())
                    .id(references.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .pubMedId(references.getAsJsonObject().get("pubMedId").isJsonNull()
                            ? null
                            : references.getAsJsonObject().getAsJsonPrimitive("pubMedId").getAsString())
                    .title(references.getAsJsonObject().getAsJsonPrimitive("title").getAsString())
                    .build());
        }
        return listReferences;
    }
}
