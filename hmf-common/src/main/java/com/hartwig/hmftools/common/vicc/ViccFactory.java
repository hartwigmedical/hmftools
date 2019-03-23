package com.hartwig.hmftools.common.vicc;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public abstract class ViccFactory {
    private static final Logger LOGGER = LogManager.getLogger(ViccFactory.class);

    private ViccFactory() {
    }

    private static StringBuilder readObjectBRCA(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //BRCA object
        StringBuilder stringToCSVBRCA = new StringBuilder();
        for (int i = 0; i < object.getAsJsonObject("brca").keySet().size(); i++) {
            List<String> keysOfBRCAObject = new ArrayList<>(object.getAsJsonObject("brca").keySet());
            stringToCSVBRCA.append(object.getAsJsonObject("brca").get(keysOfBRCAObject.get(i))).append(";"); // brca data
        }
        headerCSV.append(object.getAsJsonObject("brca").keySet()).append(";"); // header brca
        return stringToCSVBRCA;
    }

    private static StringBuilder readObjectCGI(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //CGI object
        StringBuilder stringToCSVCGI = new StringBuilder();
        return stringToCSVCGI;
    }

    private static StringBuilder readObjectCIVIC(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //CIVIC object
        StringBuilder stringToCSVCIVIC = new StringBuilder();
        return stringToCSVCIVIC;
    }

    private static StringBuilder readObjectJax(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //Jax object
        StringBuilder stringToCSVJax = new StringBuilder();
        return stringToCSVJax;
    }

    private static StringBuilder readObjectJaxTrials(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //Jax_trails object
        StringBuilder stringToCSVJaxTrials = new StringBuilder();
        return stringToCSVJaxTrials;
    }

    private static StringBuilder readObjectMolecularMatch(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //MolecularMatch object
        StringBuilder stringToCSVMolecularMatch = new StringBuilder();
        return stringToCSVMolecularMatch;
    }

    private static StringBuilder readObjectMolecularMatchTrials(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //MolecularMatchTrials object
        StringBuilder stringToCSVMolecularMatchTrials = new StringBuilder();
        return stringToCSVMolecularMatchTrials;
    }

    private static StringBuilder readObjectOncokb(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //ONCOKB object
        StringBuilder stringToCSVOncoKb = new StringBuilder();
        return stringToCSVOncoKb;
    }

    private static StringBuilder readObjectPmkb(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //PMKB object
        StringBuilder stringToCSVPmkb = new StringBuilder();
        return stringToCSVPmkb;
    }

    private static StringBuilder readObjectSage(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //SAGE object
        StringBuilder stringToCSVSage = new StringBuilder();
        return stringToCSVSage;
    }

    private static StringBuilder readObjectGenes(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //Genes object
        StringBuilder stringToCSVGenes = new StringBuilder();
        JsonArray arrayGenes = object.getAsJsonArray("genes");
        String genes = arrayGenes.toString();
        stringToCSVGenes.append(genes).append(";"); // genes data
        headerCSV.append("genes").append(";"); // header genes
        return stringToCSVGenes;
    }

    private static StringBuilder readObjectTags(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //Tags object
        StringBuilder stringToCSVTags = new StringBuilder();
        JsonArray arrayTags = object.getAsJsonArray("tags");
        String tags = arrayTags.toString();
        stringToCSVTags.append(tags).append(";"); // tags data
        headerCSV.append("tags").append(";"); // header tags
        return stringToCSVTags;
    }

    private static StringBuilder readObjectDevTags(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //dev_tags
        StringBuilder stringToCSVDevTags = new StringBuilder();
        JsonArray arrayDevTags = object.getAsJsonArray("dev_tags");
        stringToCSVDevTags.append(arrayDevTags).append(";"); // dev tags data
        headerCSV.append("dev_tags").append(";"); // header tags data
        return stringToCSVDevTags;
    }

    private static StringBuilder readObjectSource(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //Source object
        StringBuilder stringToCSVSource = new StringBuilder();
        stringToCSVSource.append(object.getAsJsonPrimitive("source")).append(";"); // source data
        headerCSV.append("source").append(";"); // header source
        return stringToCSVSource;
    }

    private static StringBuilder readObjectGeneIdentifiers(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //gene_identifiers object
        StringBuilder stringToCSVGeneIdentifiers = new StringBuilder();
        JsonArray arrayGeneIdentifiers = object.getAsJsonArray("gene_identifiers");
        JsonObject objectGeneIdentiefiers = (JsonObject) arrayGeneIdentifiers.iterator().next();
        for (int i = 0; i < objectGeneIdentiefiers.keySet().size(); i++) {
            List<String> keysOfGeneIdentifiersObject = new ArrayList<>(objectGeneIdentiefiers.keySet());
            stringToCSVGeneIdentifiers.append(objectGeneIdentiefiers.get(keysOfGeneIdentifiersObject.get(i)))
                    .append(";"); // gene identifiers data
        }
        headerCSV.append(objectGeneIdentiefiers.keySet()).append(";"); // header gene identifiers
        return stringToCSVGeneIdentifiers;
    }

    private static StringBuilder readObjectAssociation(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //association object
        int indexValue;
        List<String> headerEvidence = Lists.newArrayList();
        List<String> headerPhenotype = Lists.newArrayList();
        List<String> keysOfAssocationObject = Lists.newArrayList();

        StringBuilder stringToCSVAssociation = new StringBuilder();
        for (int i = 0; i < object.getAsJsonObject("association").keySet().size(); i++) {
            keysOfAssocationObject = new ArrayList<>(object.getAsJsonObject("association").keySet());

            if (keysOfAssocationObject.get(i).equals("description")) {
                stringToCSVAssociation.append(object.getAsJsonObject("association").get(keysOfAssocationObject.get(i)))
                        .append(";"); // association data
            } else if (keysOfAssocationObject.get(i).equals("evidence")) {
                JsonElement elementEvidence = object.getAsJsonObject("association").get("evidence");
                JsonArray arrayEvidence = elementEvidence.getAsJsonArray();
                JsonObject objectEvidence = (JsonObject) arrayEvidence.iterator().next();

                for (int a = 0; a < objectEvidence.keySet().size(); a++) {
                    List<String> keysOfEvidenceObject = new ArrayList<>(objectEvidence.keySet());
                    if (keysOfEvidenceObject.get(a).equals("evidenceType")) {
                        indexValue = keysOfEvidenceObject.indexOf("evidenceType");
                        keysOfEvidenceObject.remove(indexValue);
                        keysOfEvidenceObject.add(indexValue,
                                String.join(",", objectEvidence.get("evidenceType").getAsJsonObject().keySet()));
                        headerEvidence = keysOfEvidenceObject;
                        for (int b = 0; b < objectEvidence.get("evidenceType").getAsJsonObject().keySet().size(); b++) {
                            List<String> keysOfEvidenceTypeObject =
                                    new ArrayList<>(objectEvidence.get("evidenceType").getAsJsonObject().keySet());
                            stringToCSVAssociation.append(objectEvidence.get("evidenceType")
                                    .getAsJsonObject()
                                    .get(keysOfEvidenceTypeObject.get(b))).append(";"); // association data
                        }
                    } else {
                        stringToCSVAssociation.append(objectEvidence.get(keysOfEvidenceObject.get(a))).append(";"); // association data
                    }
                }
            } else if (keysOfAssocationObject.get(i).equals("environmentalContexts")) {
                stringToCSVAssociation.append(object.getAsJsonObject("association").get("environmentalContexts"))
                        .append(";"); // association data
            } else if (keysOfAssocationObject.get(i).equals("evidence_label")) {
                stringToCSVAssociation.append(object.getAsJsonObject("association").get(keysOfAssocationObject.get(i)))
                        .append(";"); // association data
            } else if (keysOfAssocationObject.get(i).equals("phenotype")) {
                JsonElement elementPhenotype = object.getAsJsonObject("association").get("phenotype");
                for (int a = 0; a < elementPhenotype.getAsJsonObject().keySet().size(); a++) {
                    List<String> keysOfPhenotypeObject = new ArrayList<>(elementPhenotype.getAsJsonObject().keySet());
                    if (keysOfPhenotypeObject.get(a).equals("type")) {

                        List<String> keysOfPhenotypeTypeObject =
                                new ArrayList<>(elementPhenotype.getAsJsonObject().get("type").getAsJsonObject().keySet());

                        indexValue = keysOfPhenotypeObject.indexOf("type");
                        keysOfPhenotypeObject.remove(indexValue);
                        keysOfPhenotypeObject.add(indexValue, String.join(",", keysOfPhenotypeTypeObject));
                        headerPhenotype = keysOfPhenotypeObject;

                        for (int c = 0; c < keysOfPhenotypeObject.size(); c++) {
                            stringToCSVAssociation.append(elementPhenotype.getAsJsonObject() // association data
                                    .get("type").getAsJsonObject().get(keysOfPhenotypeTypeObject.get(c))).append(";");
                        }
                    } else {
                        stringToCSVAssociation.append(elementPhenotype.getAsJsonObject().get(keysOfPhenotypeObject.get(a)))
                                .append(";"); // association data
                    }
                }
            } else if (keysOfAssocationObject.get(i).equals("oncogenic")) {
                stringToCSVAssociation.append(object.getAsJsonObject("association").get(keysOfAssocationObject.get(i)))
                        .append(";"); // association data
            }
            keysOfAssocationObject.set(keysOfAssocationObject.indexOf("evidence"), String.join(",", headerEvidence));
            keysOfAssocationObject.set(keysOfAssocationObject.indexOf("phenotype"), String.join(",", headerPhenotype));
        }
        headerCSV.append(keysOfAssocationObject).append(";"); // header features
        return stringToCSVAssociation;
    }

    private static StringBuilder readObjectFeaturesNames(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //feature_names object
        StringBuilder stringToCSVFeaturesNames = new StringBuilder();
        stringToCSVFeaturesNames.append(object.getAsJsonPrimitive("feature_names")).append(";"); // features names data
        headerCSV.append("feature_names").append(";"); // header features names
        return stringToCSVFeaturesNames;
    }

    private static StringBuilder readObjectFeatures(@NotNull JsonObject object, @NotNull StringBuilder headerCSV) {
        //features
        StringBuilder stringToCSVFeatures = new StringBuilder();
        JsonArray arrayFeatures = object.getAsJsonArray("features");
        JsonObject objectFeatures = (JsonObject) arrayFeatures.iterator().next();
        int indexValue;
        List<String> keysOfFeaturesObject = new ArrayList<>(objectFeatures.keySet());
        for (int i = 0; i < objectFeatures.keySet().size(); i++) {

            if (keysOfFeaturesObject.get(i).equals("sequence_ontology")) {
                List<String> keysOfSequenceOntologyObject = new ArrayList<>(objectFeatures.getAsJsonObject("sequence_ontology").keySet());
                indexValue = keysOfFeaturesObject.indexOf("sequence_ontology");
                keysOfFeaturesObject.set(indexValue, String.join(",", keysOfSequenceOntologyObject));

                for (int e = 0; e < objectFeatures.getAsJsonObject("sequence_ontology").keySet().size(); e++) {
                    stringToCSVFeatures.append(objectFeatures.get("sequence_ontology")
                            .getAsJsonObject()
                            .get(keysOfSequenceOntologyObject.get(e))).append(";"); // features data
                }
            } else {
                stringToCSVFeatures.append(objectFeatures.get(keysOfFeaturesObject.get(i))).append(";"); // features data
            }
        }
        headerCSV.append(keysOfFeaturesObject).append(";"); // header features
        return stringToCSVFeatures;
    }

    public static void extractAllFile(@NotNull String allJsonPath) throws IOException {
        final String csvFileName = "/Users/liekeschoenmaker/hmf/tmp/all.csv";
        PrintWriter writer = new PrintWriter(new File(csvFileName));
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(csvFileName));
        reader.setLenient(true);

        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject object = parser.parse(reader).getAsJsonObject();

            StringBuilder stringToCSVAll = new StringBuilder();
            StringBuilder headerCSV = new StringBuilder();
            StringBuilder StringToCSVSource = readObjectSource(object, headerCSV);
            StringBuilder StringToCSVGenes = readObjectGenes(object, headerCSV);

            stringToCSVAll.append(StringToCSVSource);
            stringToCSVAll.append(StringToCSVGenes);


            writer.append(headerCSV);
            writer.append("/n");
            writer.append(stringToCSVAll);
            writer.append("/n");



        }
        reader.close();
        writer.close();
    }

    public static void extractBRCAFile(@NotNull String brcaJsonPath) throws IOException {
        final String csvFileName = "/Users/liekeschoenmaker/hmf/tmp/brca.csv";
        PrintWriter writer = new PrintWriter(new File(csvFileName));
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(brcaJsonPath));
        reader.setLenient(true);

        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject object = parser.parse(reader).getAsJsonObject();

            StringBuilder stringToCSVAll = new StringBuilder();
            StringBuilder headerCSV = new StringBuilder();

            StringBuilder StringToCSVBRCA = readObjectBRCA(object, headerCSV);
            StringBuilder StringToCSVGenes = readObjectGenes(object, headerCSV);
            StringBuilder StringToCSVTags = readObjectTags(object, headerCSV);
            StringBuilder StringToCSVDevTags = readObjectDevTags(object, headerCSV);
            StringBuilder StringToCSVSource = readObjectSource(object, headerCSV);
            StringBuilder StringToCSVGeneIdentifiers = readObjectGeneIdentifiers(object, headerCSV);
            StringBuilder StringToCSVAssociation = readObjectAssociation(object, headerCSV);
            StringBuilder StringToCSVFeaturesNames = readObjectFeaturesNames(object, headerCSV);
            StringBuilder StringToCSVFeatures = readObjectFeatures(object, headerCSV);

            stringToCSVAll.append(StringToCSVBRCA)
                    .append(StringToCSVGenes)
                    .append(StringToCSVTags)
                    .append(StringToCSVDevTags)
                    .append(StringToCSVSource)
                    .append(StringToCSVGeneIdentifiers)
                    .append(StringToCSVAssociation)
                    .append(StringToCSVFeaturesNames)
                    .append(StringToCSVFeatures);

            writer.append(headerCSV);
            writer.append("\n");
            writer.append(stringToCSVAll);
            writer.append("\n");
        }
        reader.close();
        writer.close();
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
