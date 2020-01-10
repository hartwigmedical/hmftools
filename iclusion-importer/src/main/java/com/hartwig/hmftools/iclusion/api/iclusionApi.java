package com.hartwig.hmftools.iclusion.api;

import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.URL;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class iclusionApi {

    private static final Logger LOGGER = LogManager.getLogger(iclusionApi.class);

    private iclusionApi() {
    }

    public static void connectWithIclusionApi(@NotNull String iClusionLink, @NotNull String iClusionClientId,
            @NotNull String iClusionClientSecret, @NotNull String iClusionUsername, @NotNull String iClusionPassword) throws IOException {

        LOGGER.info("Connecting with iclusion API on {}", iClusionLink);
        URL url = new URL(iClusionLink); // url iclusion
        HttpURLConnection connection = (HttpURLConnection) url.openConnection();

        connection.setRequestProperty("Accept", "application/json");
        connection.setRequestProperty("content-type", "multipart/form-data");
        connection.setRequestProperty("grant_type", "password");
        connection.setRequestProperty("client_id", iClusionClientId);
        connection.setRequestProperty("client_secret", iClusionClientSecret);
        connection.setRequestProperty("username", iClusionUsername);
        connection.setRequestProperty("password", iClusionPassword);
        //  connection.setRequestProperty("Authorization", "Bearer " + "");
        connection.connect();

        LOGGER.info("headerField set cookie" + connection.getHeaderFields().get("Set-Cookie"));
        LOGGER.info("");
        LOGGER.info("headerField specific" + connection.getHeaderFields().get("Set-Cookie").get(1));



    }
}
