package com.hartwig.hmftools.bachelor;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.stream.Stream;

import org.jetbrains.annotations.Nullable;

class RunDirectory {

    private final Path mSampleDir;
    private final File mGermline;

    private String VCF_FILE_SUFFIX1 = "annotated.vcf";
    private String VCF_FILE_SUFFIX2 = "GoNLv5.vcf";
    private String VCF_FILE_SUFFIX3 = "annotated.vcf.gz";

    RunDirectory(final Path runDirectory)
    {
        mSampleDir = runDirectory;
        mGermline = findGermline();
    }

    public File germline() { return mGermline; }

    public Path sampleDir() {
        return mSampleDir;
    }

    private File findGermline()
    {
        try
        {
            // first try for tbe expected possible germline file names
            final String filePrefix = mSampleDir.toRealPath().toString() + "/" + mSampleDir.getFileName() + ".";

            if(Files.exists(Paths.get(filePrefix + VCF_FILE_SUFFIX1)))
            {
                return new File(filePrefix + VCF_FILE_SUFFIX1);
            }

            if(Files.exists(Paths.get(filePrefix + VCF_FILE_SUFFIX2)))
            {
                return new File(filePrefix + VCF_FILE_SUFFIX2);
            }

            if(Files.exists(Paths.get(filePrefix + VCF_FILE_SUFFIX3)))
            {
                return new File(filePrefix + VCF_FILE_SUFFIX3);
            }

            try (final Stream<Path> stream = Files.walk(mSampleDir.toRealPath(), 1, FileVisitOption.FOLLOW_LINKS))
            {
                return stream.filter(p -> p.toString().endsWith(VCF_FILE_SUFFIX1)
                        || p.toString().endsWith(VCF_FILE_SUFFIX2)
                        || p.toString().endsWith(VCF_FILE_SUFFIX3))
                        .map(Path::toFile)
                        .findFirst()
                        .orElse(null);
            }
        }
        catch (final IOException e)
        {
            return null;
        }
    }

    final String directoryName() { return mSampleDir.getFileName().toString(); }

    final String getPatientID()
    {
        final String[] split = mSampleDir.getFileName().toString().split("_");

        if(split.length < 2)
            return "";

        return split[split.length - 1];
    }
}
