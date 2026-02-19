package com.hartwig.hmftools.amber.blacklist;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;

import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.amber.AmberBAFFile;

public class AmberPanelTrainer
{
    private final List<String> mSamples;
    private final File mAmberDirectory;
    private final File mOutputFile;

    public AmberPanelTrainer(final List<String> samples, final File amberDirectory, final File outputFile)
    {
        if(outputFile.exists())
        {
            throw new IllegalArgumentException("Output file already exists: " + outputFile);
        }
        mSamples = samples;
        mAmberDirectory = amberDirectory;
        mOutputFile = outputFile;
    }

    public void run() throws IOException
    {
        AmberBlacklistStatistics statistics = new AmberBlacklistStatistics();
        for(String sample : mSamples)
        {
            String file = AmberBAFFile.generateAmberFilenameForWriting(mAmberDirectory.getAbsolutePath(), sample);
            if(!new File(file).exists())
            {
                file = file.substring(0, file.length() - 3);
            }
            Collection<AmberBAF> data = AmberBAFFile.read(file, true).values();
            data.forEach(baf -> statistics.record(baf, baf.tumorBAF()));
        }
        List<AmberBlacklistPoint> results = statistics.findSuspiciousPoints(mSamples.size());
        File outputFile = new File(mOutputFile.getAbsolutePath());
        AmberBlacklistFile.writeToFile(outputFile, results);
    }
}
