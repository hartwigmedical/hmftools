package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class RunsFolderProcessor {
    private static final Logger LOGGER = LogManager.getLogger(PatientDbRunner.class);

    @NotNull
    public static List<CpctRunData> getPatientRunsData(@NotNull File dir) throws IOException, ParseException {
        final List<CpctRunData> runsData = Lists.newArrayList();
        final DateFormat dateFormat = new SimpleDateFormat("yyMMdd");
        final File[] folders = dir.listFiles(File::isDirectory);
        if (folders == null) {
            throw new IOException("List files in " + dir.getName() + " returned null.");
        }
        for (final File folder : folders) {
            final String folderName = folder.getName();
            final String[] names = folderName.split("_");
            final CpctRunData cpctRunData = new CpctRunData(dateFormat.parse(names[0]), names[2], names[3], names[4]);
            runsData.add(cpctRunData);
        }
        return runsData;
    }
}
