package com.hartwig.hmftools.datamodel.finding;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.file.Path;
import java.time.LocalDate;
import java.util.ServiceLoader;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.TypeAdapterFactory;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.datamodel.LocalDateAdapter;
import com.hartwig.hmftools.datamodel.OrangeJson;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class FindingApplication
{
    public static final String APP_NAME = "Finding";

    public static final Logger LOGGER = LogManager.getLogger(FindingApplication.class);

    public static final String VERSION = FindingApplication.class.getPackage().getImplementationVersion();

    public static void main(String[] args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        FindingApplicationConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        FindingApplicationConfig config = new FindingApplicationConfig(configBuilder);
        new FindingApplication(config).run();
    }

    private final FindingApplicationConfig config;

    private FindingApplication(final FindingApplicationConfig config)
    {
        this.config = config;
    }

    private void run() throws Exception
    {
        LOGGER.info("Generating finding json");

        OrangeRecord orangeRecord = OrangeJson.getInstance().read(config.OrangeJsonPath);

        FindingRecord findingRecord = FindingRecordFactory.fromOrangeRecord(orangeRecord,
                config.ClinicalTranscriptsPath != null ? Path.of(config.ClinicalTranscriptsPath) : null,
                Path.of(config.DriverGenePath));

        GsonBuilder gsonBuilder = new GsonBuilder();
        gsonBuilder.setPrettyPrinting();
        for (TypeAdapterFactory factory : ServiceLoader.load(TypeAdapterFactory.class)) {
            gsonBuilder.registerTypeAdapterFactory(factory);
        }

        Gson gson = gsonBuilder.serializeNulls().serializeSpecialFloatingPointValues()
                .registerTypeAdapter(LocalDate.class, new LocalDateAdapter())
                .create();

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(config.FindingJsonPath))) {
            gson.toJson(findingRecord, FindingRecord.class, writer);
        }

        LOGGER.info("Done!");
    }
}
