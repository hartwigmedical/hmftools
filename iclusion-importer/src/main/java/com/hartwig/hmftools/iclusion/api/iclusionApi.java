package com.hartwig.hmftools.iclusion.api;

import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class iclusionApi {

    private static final Logger LOGGER = LogManager.getLogger(iclusionApi.class);

    private iclusionApi() {
    }

    public static String connectWithIclusionApi(@NotNull String iClusionLink, @NotNull String iClusionClientId,
            @NotNull String iClusionClientSecret, @NotNull String iClusionUsername, @NotNull String iClusionPassword) throws IOException {

        LOGGER.info("Connecting with iclusion API on {}", iClusionLink);
        URL url = new URL(iClusionLink); // url iclusion
        HttpURLConnection connection = (HttpURLConnection) url.openConnection();

      //  connection.setRequestMethod("POST");
       // connection.setRequestProperty("Accept", "application/json");
       // connection.setRequestProperty("content-type", "multipart/form-data");
        connection.setRequestProperty("grant_type", "password");
        connection.setRequestProperty("client_id", iClusionClientId);
        connection.setRequestProperty("client_secret", iClusionClientSecret);
        connection.setRequestProperty("username", iClusionUsername);
        connection.setRequestProperty("password", iClusionPassword);

        List<String> setCookie = connection.getHeaderFields().get("Set-Cookie");
        int tokenEnd = setCookie.get(1).indexOf(";");

        //TODO create access with token
        return setCookie.get(1).substring(0, tokenEnd).split("=")[1];

    }
}
