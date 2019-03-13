package com.hartwig.hmftools.common.vicc;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public abstract class ViccFactory {
    private static final Logger LOGGER = LogManager.getLogger(ViccFactory.class);

    private ViccFactory() {
    }

    public static void extractAllFile(@NotNull String allJsonPath) throws IOException {
    }

    public static void extractBRCAFile(@NotNull String brcaJsonPath) throws IOException {
        final String csvFileName = "/data/experiments/knowledgebase_vicckb/brca.csv";
        PrintWriter writer = new PrintWriter(new File(csvFileName));
        writer.append("Variant_frequency_LOVD;ClinVarAccession_ENIGMA");
        final Object object = new JsonParser().parse(new FileReader(brcaJsonPath));
        JsonArray jsonObject = (JsonArray) object;

        for (Object dataJson : jsonObject) {
            JsonObject brca = (JsonObject) dataJson;

            LOGGER.info(brca.getAsJsonObject("brca").get("Variant_frequency_LOVD"));
            LOGGER.info(brca.getAsJsonObject("brca").get("ClinVarAccession_ENIGMA"));
            LOGGER.info(brca.getAsJsonObject("gene_identifiers").get("symbol"));
            LOGGER.info(brca.getAsJsonObject("gene_identifiers").get("entrez_id"));
            LOGGER.info(brca.getAsJsonObject("gene_identifiers").get("ensembl_gene_id"));
            BRCA brca1 =  ImmutableBRCA.builder()
                    .variantFrequencyLOVD(brca.getAsJsonObject("brca").get("Variant_frequency_LOVD").toString())
                    .clinVarAccessionENIGMA(brca.getAsJsonObject("brca").get("ClinVarAccession_ENIGMA").toString()).build();
            LOGGER.info(brca1.toString());
            writeToCsvFile(brca1.toString(), writer);
        }
    }

    private static void writeToCsvFile(@NotNull String csvFormatString, @NotNull PrintWriter writer) throws IOException {
        writer.append(csvFormatString);



    }

    public static void extractCgiFile(@NotNull String cgiJsonPath) throws IOException {
    }

    public static void extractCivicFile(@NotNull String civicJsonPath) throws IOException {
    }

    public static void extractJaxFile(@NotNull String jaxJsonPath) throws IOException {
    }

    public static void extractJaxTrialsFile(@NotNull String jaxTrialsJsonPath) throws IOException {
    }

    public static void extractMolecularMatchFile(@NotNull String molecularMatchJsonPath) throws IOException {
    }

    public static void extractMolecularMatchTrailsFile(@NotNull String molecularMatchTrialsJsonPath) throws IOException {
    }

    public static void extractOncokbFile(@NotNull String oncokbJsonPath) throws IOException {
    }

    public static void extractPmkbFile(@NotNull String pmkbJsonPath) throws IOException {
    }

    public static void extractSageFile(@NotNull String sageJsonPath) throws IOException {
    }

}
