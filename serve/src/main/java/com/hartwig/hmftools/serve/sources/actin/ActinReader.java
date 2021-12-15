package com.hartwig.hmftools.serve.sources.actin;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinFileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ActinReader {

    private static final Logger LOGGER = LogManager.getLogger(ActinReader.class);

    private ActinReader(){
    }

    @NotNull
    public static List<ActinEntry> read(@NotNull String actinTrialTsv) throws IOException {
        LOGGER.info("Reading ACTIN trial database from {}", actinTrialTsv);
        List<ActinEntry> actinEntries = ActinFileReader.read(actinTrialTsv);
        LOGGER.info(" Read {} entries", actinEntries.size());

        return actinEntries;
    }
}
