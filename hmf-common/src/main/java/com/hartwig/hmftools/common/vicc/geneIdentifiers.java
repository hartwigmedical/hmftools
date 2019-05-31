package com.hartwig.hmftools.common.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class geneIdentifiers {

    public static void readObjectGeneIdentifiers(@NotNull JsonObject object) {
        JsonArray arrayGeneIdentifiers = object.getAsJsonArray("gene_identifiers");

        for (int j = 0; j < arrayGeneIdentifiers.size(); j++) {
            JsonObject objectGeneIdentiefiers = (JsonObject) arrayGeneIdentifiers.get(j);
            for (int i = 0; i < objectGeneIdentiefiers.keySet().size(); i++) {
                List<String> keysOfGeneIdentifiersObject = new ArrayList<>(objectGeneIdentiefiers.keySet());
                if (keysOfGeneIdentifiersObject.get(i).equals("symbol")) {
                    JsonElement symbol = objectGeneIdentiefiers.get(keysOfGeneIdentifiersObject.get(i));
                    // TODO: add to sql
                } else if (keysOfGeneIdentifiersObject.get(i).equals("entrez_id")) {
                    JsonElement entrezId = objectGeneIdentiefiers.get(keysOfGeneIdentifiersObject.get(i));
                    // TODO: add to sql
                } else if (keysOfGeneIdentifiersObject.get(i).equals("ensembl_gene_id")) {
                    JsonElement ensemblGeneID = objectGeneIdentiefiers.get(keysOfGeneIdentifiersObject.get(i));
                    // TODO: add to sql
                }
            }
        }
    }
}
