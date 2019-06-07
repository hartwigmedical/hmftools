package com.hartwig.hmftools.vicc;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.vicc.dao.ViccDAO;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ViccJsonToSQLImporter {
    private static final Logger LOGGER = LogManager.getLogger(ViccJsonToSQLImporter.class);

    public static void main(final String... args) throws IOException, SQLException {
        LOGGER.info("Attempting to load up the VICC json all file into a sql database");

        final String baseDir =
                System.getProperty("user.home") + File.separator + "hmf" + File.separator + "projects" + File.separator + "vicc";
        final String inputFile = baseDir + File.separator + "all.json";

        List<ViccEntry> viccEntries = ViccFactory.readViccKnowledgebaseJsonFile(inputFile);

        ViccDAO viccDAO = ViccDAO.connectToViccDAO("build", "build", "jdbc:mysql://localhost:3306/vicc_db?serverTimezone=CET");

        viccDAO.deleteAll();
        for (ViccEntry viccEntry : viccEntries) {
            viccDAO.writeViccEntry(viccEntry);
        }
    }
}
