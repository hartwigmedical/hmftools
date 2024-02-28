package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadVirusInterpreter
{
    private static final String VIRUS_ANNOTATION_TSV = "virus_annotation_tsv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addPath(VIRUS_ANNOTATION_TSV, true, "Path towards the virus annotations TSV file");

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);

        String virusAnnotationTsv = configBuilder.getValue(VIRUS_ANNOTATION_TSV);

        try(DatabaseAccess dbWriter = databaseAccess(configBuilder))
        {
            LOGGER.info("Reading virus annotation TSV {}", virusAnnotationTsv);
            List<AnnotatedVirus> virusAnnotations = AnnotatedVirusFile.read(virusAnnotationTsv);
            LOGGER.info(" Read {} virus annotation", virusAnnotations.size());

            LOGGER.info("Writing virus annotations into database for {}", sample);
            dbWriter.writeVirusInterpreter(sample, virusAnnotations);

            LOGGER.info("Complete");
        }
        catch (Exception e)
        {
            LOGGER.error("Failed to load virus annotations", e);
            System.exit(1);
        }
    }
}
