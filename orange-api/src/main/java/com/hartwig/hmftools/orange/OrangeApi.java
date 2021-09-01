package com.hartwig.hmftools.orange;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Stream;

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

    @NotNull
    public static OrangeReport read(@NotNull String jsonFile) {
        String json = readLineByLine(jsonFile);
        return new GsonBuilder().serializeNulls().serializeSpecialFloatingPointValues().create().fromJson(json, OrangeReport.class);
    }

    @NotNull
    private static String readLineByLine(String filePath) {
        StringBuilder contentBuilder = new StringBuilder();

        try (Stream<String> stream = Files.lines(Paths.get(filePath), StandardCharsets.UTF_8)) {
            stream.forEach(s -> contentBuilder.append(s).append("\n"));
        } catch (IOException e) {
            e.printStackTrace();
        }

        return contentBuilder.toString();
    }
}
