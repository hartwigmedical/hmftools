package com.hartwig.hmftools.patientdb.readers;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.data.CpctRunData;

import org.jetbrains.annotations.NotNull;

public class RunsFolderReader {
    @NotNull
    public static List<CpctRunData> getPatientRunsData(@NotNull final File dir) throws IOException, ParseException {
        final List<CpctRunData> runsData = Lists.newArrayList();
        final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyMMdd");
        final File[] folders = dir.listFiles(File::isDirectory);
        if (folders == null) {
            throw new IOException("List files in " + dir.getName() + " returned null.");
        }
        for (final File folder : folders) {
            final String folderName = folder.getName();
            final String[] names = folderName.split("_");
            LocalDate uploadDate = LocalDate.parse(names[0], dateFormatter);
            final CpctRunData cpctRunData = new CpctRunData(folderName, uploadDate, names[2], names[3], names[4]);
            runsData.add(cpctRunData);
        }
        return runsData;
    }
}
