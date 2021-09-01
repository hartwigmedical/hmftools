package com.hartwig.hmftools.orange;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import com.google.gson.GsonBuilder;

import org.jetbrains.annotations.NotNull;

public final class OrangeApi {

    private OrangeApi() {
    }

    public static void write(@NotNull OrangeReport report, @NotNull String jsonFile) throws IOException {
        String json = new GsonBuilder().serializeNulls().serializeSpecialFloatingPointValues().create().toJson(report);
        BufferedWriter writer = new BufferedWriter(new FileWriter(jsonFile));

        writer.write(json);
        writer.close();
    }
}
