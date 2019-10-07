package com.hartwig.hmftools.vicc.reader;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.sage.ImmutableSage;
import com.hartwig.hmftools.vicc.datamodel.sage.Sage;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class SageObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(SageObjectFactory.class);

    private static final List<Integer> EXPECTED_SAGE_ELEMENT_SIZES = Lists.newArrayList(8);

    private SageObjectFactory() {
    }

    @NotNull
    static Sage create(@NotNull JsonObject objectSage) {
        Set<String> keysSage = objectSage.keySet();

        if (!EXPECTED_SAGE_ELEMENT_SIZES.contains(keysSage.size())) {
            LOGGER.warn("Found " + keysSage.size() + " in sage rather than the expected " + EXPECTED_SAGE_ELEMENT_SIZES);
            LOGGER.warn(keysSage);
        }

        return ImmutableSage.builder()
                .entrezId(objectSage.getAsJsonPrimitive("entrez_id").getAsString())
                .clinicalManifestation(objectSage.getAsJsonPrimitive("clinical_manifestation").getAsString())
                .publicationUrl(objectSage.getAsJsonPrimitive("publication_url").getAsString())
                .germlineOrSomatic(objectSage.getAsJsonPrimitive("germline_or_somatic").getAsString())
                .evidenceLabel(objectSage.getAsJsonPrimitive("evidence_label").getAsString())
                .drugLabels(objectSage.getAsJsonPrimitive("drug_labels").getAsString())
                .responseType(objectSage.getAsJsonPrimitive("response_type").getAsString())
                .gene(objectSage.getAsJsonPrimitive("gene").getAsString())
                .build();
    }
}
