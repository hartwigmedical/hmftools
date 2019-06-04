package com.hartwig.hmftools.vicc;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ViccFactory {
    private static final Logger LOGGER = LogManager.getLogger(ViccFactory.class);

    private ViccFactory() {
    }

    public static void extractAllFile(@NotNull String allJsonPath) throws IOException {
        final String csvFileName = "/Users/liekeschoenmaker/hmf/tmp/all.csv";
        PrintWriter writer = new PrintWriter(new File(csvFileName));
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(allJsonPath));
        reader.setLenient(true);
        int index = 1;

        while (reader.peek() != JsonToken.END_DOCUMENT && index < 5) {
            JsonObject object = parser.parse(reader).getAsJsonObject();
            source.readObjectSource(object);
            genes.readObjectGenes(object);
            tags.readObjectTags(object);
            devTags.readObjectDevTags(object);
            geneIdentifiers.readObjectGeneIdentifiers(object);
            featuresNames.readObjectFeaturesNames(object);

            sage.readObjectSage(object);
            StringBuilder stringToCSVPmkb = pmkb.readObjectPmkb(object);
            StringBuilder stringToCSVBrca = brca.readObjectBRCA(object);
            StringBuilder stringToCSVCGI = cgi.readObjectCGI(object);
            StringBuilder stringToCSVOncokb = oncokb.readObjectOncokb(object);
            StringBuilder stringToCSVJax = jax.readObjectJax(object);
            StringBuilder stringToCSVJaxTrials = jaxTrials.readObjectJaxTrials(object);

            StringBuilder stringToCSVAssociation = association.readObjectAssociation(object); //TODO check fields for every db
            StringBuilder stringToCSVFeatures = features.readObjectFeatures(object); //TODO check fields for every db
            StringBuilder stringToCSVCIVIC = civic.readObjectCIVIC(object); //TODO check fields for every db
            StringBuilder stringToCSVMolecularMatch = molecularMatch.readObjectMolecularMatch(object); //TODO check fields for every db
            StringBuilder stringToCSVMolecularMatchTrials =
                    molecularMatchTrial.readObjectMolecularMatchTrials(object); //TODO check fields for every db
            index++;

        }
        reader.close();
        writer.close();
    }
}
