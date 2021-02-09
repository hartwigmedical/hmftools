package com.hartwig.hmftools.ckb.reader;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.ckb.datamodel.CkbJsonObject;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public abstract class CkbJsonDirectoryReader<T extends CkbJsonObject> {

    private static final Logger LOGGER = LogManager.getLogger(CkbJsonDirectoryReader.class);

    @Nullable
    private final Integer maxFilesToRead;

    public CkbJsonDirectoryReader(@Nullable final Integer maxFilesToRead) {
        this.maxFilesToRead = maxFilesToRead;
    }

    @NotNull
    public List<T> read(@NotNull String dir) throws IOException {
        List<T> entries = Lists.newArrayList();
        File[] files = new File(dir).listFiles();

        LOGGER.debug(" {} files found in directory {}", files.length, dir);

        int currentFileIndex = 0;
        while (currentFileIndex < files.length && (maxFilesToRead == null || currentFileIndex < maxFilesToRead)) {
            JsonParser parser = new JsonParser();
            JsonReader reader = new JsonReader(new FileReader(files[currentFileIndex]));
            reader.setLenient(true);

            while (reader.peek() != JsonToken.END_DOCUMENT) {
                entries.add(read(parser.parse(reader).getAsJsonObject()));
            }
            reader.close();
            currentFileIndex++;
        }

        LOGGER.debug("  Done reading {} files ", currentFileIndex);
        return entries;
    }

    @NotNull
    protected abstract T read(@NotNull JsonObject object);
}
