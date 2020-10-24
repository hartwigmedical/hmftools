package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.common.utils.json.JsonFunctions.nullableString;
import static com.hartwig.hmftools.common.utils.json.JsonFunctions.string;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.ImmutableJaxTrials;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.ImmutableJaxTrialsIndication;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.ImmutableJaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.ImmutableJaxTrialsTherapy;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsIndication;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsTherapy;

import org.jetbrains.annotations.NotNull;

final class JaxTrialsObjectFactory {

    private JaxTrialsObjectFactory() {
    }

    @NotNull
    static JaxTrials create(@NotNull JsonObject jaxTrialsObject) {
        ViccDatamodelCheckerFactory.jaxTrialsEntryChecker().check(jaxTrialsObject);

        return ImmutableJaxTrials.builder()
                .nctId(string(jaxTrialsObject, "nctId"))
                .title(string(jaxTrialsObject, "title"))
                .variantRequirements(string(jaxTrialsObject, "variantRequirements"))
                .molecularProfiles(createMolecularProfiles(jaxTrialsObject.getAsJsonArray("variantRequirementDetails")))
                .indications(createIndications(jaxTrialsObject.getAsJsonArray("indications")))
                .therapies(createTherapies(jaxTrialsObject.getAsJsonArray("therapies")))
                .gender(nullableString(jaxTrialsObject, "gender"))
                .recruitment(string(jaxTrialsObject, "recruitment"))
                .phase(string(jaxTrialsObject, "phase"))
                .sponsors(string(jaxTrialsObject, "sponsors"))
                .updateDate(string(jaxTrialsObject, "updateDate"))
                .build();
    }

    @NotNull
    private static List<JaxTrialsMolecularProfile> createMolecularProfiles(@NotNull JsonArray variantRequirementDetailsArray) {
        List<JaxTrialsMolecularProfile> molecularProfileList = Lists.newArrayList();
        DatamodelChecker variantRequirementDetailsChecker = ViccDatamodelCheckerFactory.jaxTrialsVariantRequirementDetailsChecker();
        DatamodelChecker molecularProfileChecker = ViccDatamodelCheckerFactory.jaxTrialsMolecularProfileChecker();

        for (JsonElement variantRequirementDetailsElement : variantRequirementDetailsArray) {
            JsonObject variantRequirementDetailsObject = variantRequirementDetailsElement.getAsJsonObject();
            variantRequirementDetailsChecker.check(variantRequirementDetailsObject);

            JsonObject molecularProfileObject = variantRequirementDetailsObject.getAsJsonObject("molecularProfile");
            molecularProfileChecker.check(molecularProfileObject);

            molecularProfileList.add(ImmutableJaxTrialsMolecularProfile.builder()
                    .requirementType(string(variantRequirementDetailsObject, "requirementType"))
                    .profileName(string(molecularProfileObject, "profileName"))
                    .id(string(molecularProfileObject, "id"))
                    .build());
        }
        return molecularProfileList;
    }

    @NotNull
    private static List<JaxTrialsIndication> createIndications(@NotNull JsonArray indicationArray) {
        List<JaxTrialsIndication> indicationList = Lists.newArrayList();
        DatamodelChecker indicationChecker = ViccDatamodelCheckerFactory.jaxTrialsIndicationChecker();

        for (JsonElement indicationElement : indicationArray) {
            JsonObject indicationObject = indicationElement.getAsJsonObject();
            indicationChecker.check(indicationObject);

            indicationList.add(ImmutableJaxTrialsIndication.builder()
                    .name(string(indicationObject, "name"))
                    .source(string(indicationObject, "source"))
                    .id(string(indicationObject, "id"))
                    .build());
        }
        return indicationList;
    }

    @NotNull
    private static List<JaxTrialsTherapy> createTherapies(@NotNull JsonArray therapyArray) {
        List<JaxTrialsTherapy> therapyList = Lists.newArrayList();
        DatamodelChecker therapyChecker = ViccDatamodelCheckerFactory.jaxTrialsTherapyChecker();

        for (JsonElement therapyElement : therapyArray) {
            JsonObject therapyObject = therapyElement.getAsJsonObject();
            therapyChecker.check(therapyObject);

            therapyList.add(ImmutableJaxTrialsTherapy.builder()
                    .therapyName(string(therapyObject, "therapyName"))
                    .id(string(therapyObject, "id"))
                    .build());
        }
        return therapyList;
    }
}
