package com.hartwig.hmftools.amber;

import java.io.IOException;

import com.hartwig.hmftools.common.amber.AmberBAFFile;
import com.hartwig.hmftools.common.pcf.PCFFile;
import com.hartwig.hmftools.common.r.RExecutor;

import org.jetbrains.annotations.NotNull;

public class BAFSegmentation {

    @NotNull
    private final String outputDirectory;

    public BAFSegmentation(@NotNull final String outputDirectory) {
        this.outputDirectory = outputDirectory;
    }

    public void applySegmentation(@NotNull final String tumor) throws InterruptedException, IOException {
        final String ratioFile = AmberBAFFile.generateAmberFilenameForReading(outputDirectory, tumor);
        final String pcfFile = PCFFile.generateBAFFilename(outputDirectory, tumor);
        int result = RExecutor.executeFromClasspath("r/bafSegmentation.R", ratioFile, pcfFile);
        if (result != 0) {
            throw new IOException("R execution failed. Unable to complete segmentation.");
        }
    }
    
}
