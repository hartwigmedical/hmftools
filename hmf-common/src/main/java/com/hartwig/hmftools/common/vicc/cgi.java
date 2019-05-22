package com.hartwig.hmftools.common.vicc;

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
}
