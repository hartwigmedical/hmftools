package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.jsonArrayToStringList;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.cgi.Cgi;
import com.hartwig.hmftools.vicc.datamodel.cgi.ImmutableCgi;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class CgiObjectFactory {

    private static final Logger LOGGER = LogManager.getLogger(CgiObjectFactory.class);

    private static final List<Integer> EXPECTED_CGI_ELEMENT_SIZES = Lists.newArrayList(23);

    private CgiObjectFactory() {
    }

    @NotNull
    static Cgi create(@NotNull JsonObject objectCgi) {
        Set<String> keysCgi = objectCgi.keySet();

        if (!EXPECTED_CGI_ELEMENT_SIZES.contains(keysCgi.size())) {
            LOGGER.warn("Found " + keysCgi.size() + " in cgi rather than the expected " + EXPECTED_CGI_ELEMENT_SIZES);
            LOGGER.warn(keysCgi);
        }

        return ImmutableCgi.builder()
                .targeting(objectCgi.getAsJsonPrimitive("Targeting").getAsString())
                .source(objectCgi.getAsJsonPrimitive("Source").getAsString())
                .cDNA(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("cDNA"))))
                .primary_tumor_type(objectCgi.getAsJsonPrimitive("Primary Tumor type").getAsString())
                .individual_mutation(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("individual_mutation"))))
                .drugsFullName(objectCgi.getAsJsonPrimitive("Drug full name").getAsString())
                .curator(objectCgi.getAsJsonPrimitive("Curator").getAsString())
                .drug_family(objectCgi.getAsJsonPrimitive("Drug family").getAsString())
                .alteration(objectCgi.getAsJsonPrimitive("Alteration").getAsString())
                .drug(objectCgi.getAsJsonPrimitive("Drug").getAsString())
                .biomarker(objectCgi.getAsJsonPrimitive("Biomarker").getAsString())
                .gDNA(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("gDNA"))))
                .drug_status(objectCgi.getAsJsonPrimitive("Drug status").getAsString())
                .gene(objectCgi.getAsJsonPrimitive("Gene").getAsString())
                .transcript(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("transcript"))))
                .strand(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("strand"))))
                .info(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("info"))))
                .assay_type(objectCgi.getAsJsonPrimitive("Assay type").getAsString())
                .alteration_type(objectCgi.getAsJsonPrimitive("Alteration type").getAsString())
                .region(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("region"))))
                .evidence_level(objectCgi.getAsJsonPrimitive("Evidence level").getAsString())
                .association(objectCgi.getAsJsonPrimitive("Association").getAsString())
                .metastatic_Tumor_Type(objectCgi.getAsJsonPrimitive("Metastatic Tumor Type").getAsString())
                .build();
    }
}
