package com.hartwig.hmftools.common.vicc;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;

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
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(brcaJsonPath));
        reader.setLenient(true);

        int index = 0;

        try {
            while (reader.hasNext()) {
                LOGGER.info(index);
                JsonObject object = parser.parse(reader).getAsJsonObject();
                LOGGER.info(object.keySet());
                if (index == 0) {
                    //                    LOGGER.info(object.getAsJsonObject("features").keySet()); // Set of keys
                    writer.append(object.getAsJsonObject("brca").keySet().toString());
                    writer.append(",genes");
                    writer.append(",tags");
                    writer.append(",source");
                }
                StringBuilder StringToCSVBRCA = new StringBuilder();

                //BRCA
                for (int i = 0; i < object.getAsJsonObject("brca").keySet().size(); i++) {
                    List<String> keysOfBRCAObject = new ArrayList<>(object.getAsJsonObject("brca").keySet());
                    StringToCSVBRCA.append(object.getAsJsonObject("brca").get(keysOfBRCAObject.get(i)))
                            .append(";"); // merge 1 object to string
                }

                //Genes
                JsonArray arrayGenes = object.getAsJsonArray("genes");
                String genes = arrayGenes.toString();
                LOGGER.info(genes); // TODO: remove [] from string
                StringToCSVBRCA.append(genes).append(";");// merge 1 object to string

                //Tags
                JsonArray arrayTags = object.getAsJsonArray("tags");
                String tags = arrayTags.toString();
                StringToCSVBRCA.append(tags).append(";");// merge 1 object to string

                //Source
                StringToCSVBRCA.append(object.getAsJsonPrimitive("source")).append(";");// merge 1 object to string

                //gene_identifiers
                JsonArray arrayGeneIdentifiers = object.getAsJsonArray("gene_identifiers");
                LOGGER.info(arrayGeneIdentifiers);
                JsonObject objectGeneIdentiefiers = (JsonObject) arrayGeneIdentifiers.iterator().next();
                for (int i = 0; i < objectGeneIdentiefiers.keySet().size(); i++) {
                    List<String> keysOfGeneIdentifiersObject = new ArrayList<>(objectGeneIdentiefiers.keySet());
                    StringToCSVBRCA.append(objectGeneIdentiefiers.get(keysOfGeneIdentifiersObject.get(i))).append(";");
                }

                //association
                for (int i = 0; i < object.getAsJsonObject("association").keySet().size(); i++) {
                    List<String> keysOfAssocationObject = new ArrayList<>(object.getAsJsonObject("association").keySet());
                    if (keysOfAssocationObject.get(i).equals("description")) {
                        JsonElement description = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                        LOGGER.info(description);
                    } else if (keysOfAssocationObject.get(i).equals("evidence")) {
                        JsonElement elementEvidence = object.getAsJsonObject("association").get("evidence");
                        JsonArray arrayEvidence = elementEvidence.getAsJsonArray();
                        JsonObject objectEvidence = (JsonObject) arrayEvidence.iterator().next();
                        for (int a = 0; a < objectEvidence.keySet().size(); a++) {
                            List<String> keysOfEvidenceObject = new ArrayList<>(objectEvidence.keySet());
                            if (keysOfEvidenceObject.get(a).equals("evidenceType")) {
                                for (int b = 0; b < objectEvidence.get("evidenceType").getAsJsonObject().keySet().size(); b++) {
                                    List<String> keysOfEvidenceTypeObject =
                                            new ArrayList<>(objectEvidence.get("evidenceType").getAsJsonObject().keySet());
                                    LOGGER.info(objectEvidence.get("evidenceType").getAsJsonObject().get(keysOfEvidenceTypeObject.get(b)));
                                }
                            } else {
                                JsonElement info = objectEvidence.get(keysOfEvidenceObject.get(a));
                                LOGGER.info(info);
                            }
                        }
                    } else if (keysOfAssocationObject.get(i).equals("environmentalContexts")) {
                        JsonElement elementEnvironmentalContexts = object.getAsJsonObject("association").get("environmentalContexts");
                        LOGGER.info(elementEnvironmentalContexts);
                    } else if (keysOfAssocationObject.get(i).equals("evidence_label")) {
                        JsonElement evidence_label = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                        LOGGER.info(evidence_label);
                    } else if (keysOfAssocationObject.get(i).equals("phenotype")) {
                        JsonElement elementPhenotype = object.getAsJsonObject("association").get("phenotype");
                        LOGGER.info(elementPhenotype.getAsJsonObject().keySet());
                        for (int a = 0; a < elementPhenotype.getAsJsonObject().keySet().size(); a++) {
                            List<String> keysOfPhenotypeObject = new ArrayList<>(elementPhenotype.getAsJsonObject().keySet());
                            LOGGER.info(keysOfPhenotypeObject.get(a));
                            if (keysOfPhenotypeObject.get(a).equals("type")) {
                                List<String> keysOfPhenotypeTypeObject =
                                        new ArrayList<>(elementPhenotype.getAsJsonObject().get("type").getAsJsonObject().keySet());
                                for (int c = 0; c < keysOfPhenotypeObject.size(); c++) {
                                    LOGGER.info(elementPhenotype.getAsJsonObject()
                                            .get("type")
                                            .getAsJsonObject()
                                            .get(keysOfPhenotypeTypeObject.get(c)));
                                }
                            } else {
                                JsonElement elemtPhenotype = elementPhenotype.getAsJsonObject().get(keysOfPhenotypeObject.get(a));
                            }
                        }
                    } else if (keysOfAssocationObject.get(i).equals("oncogenic")) {
                        JsonElement oncogenic = object.getAsJsonObject("association").get(keysOfAssocationObject.get(i));
                        LOGGER.info(oncogenic);
                    }
                }

                //feature_names
                LOGGER.info(object.getAsJsonPrimitive("feature_names"));

                //dev_tags
                JsonArray arrayDevTags = object.getAsJsonArray("dev_tags");
                StringToCSVBRCA.append(arrayDevTags).append(";");// merge 1 object to string

                //features
                JsonArray arrayFeatures = object.getAsJsonArray("features");
                JsonObject objectFeatures = (JsonObject) arrayFeatures.iterator().next();
                for (int i = 0; i < objectFeatures.keySet().size(); i++) {
                    List<String> keysOfFeaturesObject = new ArrayList<>(objectFeatures.keySet());

                    if (keysOfFeaturesObject.get(i).equals("sequence_ontology")) {
                        List<String> keysOfSequenceOntologyObject =
                                new ArrayList<>(objectFeatures.getAsJsonObject("sequence_ontology").keySet());

                        for (int e = 0; e < objectFeatures.getAsJsonObject("sequence_ontology").keySet().size(); e++) {
                            LOGGER.info(objectFeatures.get("sequence_ontology").getAsJsonObject().get(keysOfSequenceOntologyObject.get(e)));
                        }
                    } else {
                        LOGGER.info(objectFeatures.get(keysOfFeaturesObject.get(i)));
                    }
                }
                index++;
                writer.append(StringToCSVBRCA);
                writer.append("\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        writer.close();
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
