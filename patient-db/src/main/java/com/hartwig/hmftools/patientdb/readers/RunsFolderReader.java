package com.hartwig.hmftools.patientdb.readers;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.io.FolderChecker;
import com.hartwig.hmftools.patientdb.context.ProductionRunContextFactory;
import com.hartwig.hmftools.patientdb.context.RunContext;

import org.jetbrains.annotations.NotNull;

public final class RunsFolderReader {
    
    private RunsFolderReader() {
    }

    @NotNull
    public static List<RunContext> extractRunContexts(@NotNull File dir) throws IOException {
        List<RunContext> runContexts = Lists.newArrayList();
        File[] folders = dir.listFiles(File::isDirectory);
        if (folders == null) {
            throw new IOException("List files in " + dir.getName() + " returned null.");
        }

        for (File folder : folders) {
            String runDirectory = FolderChecker.build().checkFolder(folder.getPath());
            runContexts.add(ProductionRunContextFactory.fromRunDirectory(runDirectory));
        }
        return runContexts;
    }
}
