package com.hartwig.hmftools.common.vicc;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class genes {

    public static void readObjectGenes(@NotNull JsonObject object) {
        //Genes object
        StringBuilder stringToCSVGenes = new StringBuilder();
        JsonArray arrayGenes = object.getAsJsonArray("genes");
        String genes = arrayGenes.toString();
        genes = genes.substring(1, genes.length() - 1);

        // TODO: add to sql
    }
}
