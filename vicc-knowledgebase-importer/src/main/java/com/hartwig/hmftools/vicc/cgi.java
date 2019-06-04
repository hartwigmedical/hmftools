package com.hartwig.hmftools.vicc;

import java.util.ArrayList;
import java.util.List;

import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;

public class cgi {

    public static StringBuilder readObjectCGI(@NotNull JsonObject object) {
        //CGI object
        StringBuilder stringToCSVCGI = new StringBuilder();
        StringBuilder stringInfo = new StringBuilder();
        if (object.getAsJsonObject("cgi") != null) {
            List<String> keysOfCGI = new ArrayList<>(object.getAsJsonObject("cgi").keySet());
            for (int i = 0; i < keysOfCGI.size(); i++) {
                if (keysOfCGI.get(i).equals("info")) {
                    JsonElement elementInfo = object.getAsJsonObject("cgi").get(keysOfCGI.get(i));
                    for (int z = 0; z < elementInfo.getAsJsonArray().size(); z++) {
                        String info = elementInfo.getAsJsonArray().get(z).toString().replaceAll(";", ":");
                        stringInfo.append(info).append(",");
                    }
                    stringToCSVCGI.append(stringInfo).append(";");
                } else {
                    stringToCSVCGI.append(object.getAsJsonObject("cgi").get(keysOfCGI.get(i))).append(";");
                }
            }
        } else {
            stringToCSVCGI.append(";;;;;;;;;;;;;;;;;;;;;;; ");

        }
        return stringToCSVCGI;
    }

    public static StringBuilder readObjectCGISpecificFields(@NotNull JsonObject object) {

        //removed data fields:
        // Targeting;cDNA;individual_mutation;Curator;Alteration;Drug;gDNA;Drug status;cgi.transcript;strand;
        // info;Assay type;Alteration type;region;Metastatic Tumor Type;

        //CGI object
        StringBuilder stringToCSVCGI = new StringBuilder();
        StringBuilder stringInfo = new StringBuilder();
        if (object.getAsJsonObject("cgi") != null) {
            List<String> keysOfCGI = new ArrayList<>(object.getAsJsonObject("cgi").keySet());
            for (int i = 0; i < keysOfCGI.size(); i++) {
                 if (keysOfCGI.get(i).equals("Source") || keysOfCGI.get(i).equals("Primary Tumor type") || keysOfCGI.get(i)
                        .equals("Drug full name") || keysOfCGI.get(i).equals("Drug family") || keysOfCGI.get(i).equals("Alteration")
                        || keysOfCGI.get(i).equals("Biomarker") || keysOfCGI.get(i).equals("Gene") || keysOfCGI.get(i)
                        .equals("Evidence level") || keysOfCGI.get(i).equals("Association")) {
                     stringToCSVCGI.append(object.getAsJsonObject("cgi").get(keysOfCGI.get(i))).append(";");
                 }
            }
        }
        return stringToCSVCGI;
    }
}
