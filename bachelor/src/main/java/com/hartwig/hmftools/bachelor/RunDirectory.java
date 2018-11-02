package com.hartwig.hmftools.bachelor;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.Stream;

import org.jetbrains.annotations.Nullable;

class RunDirectory {

    private final Path mSampleDir;
    private final File mGermline;
    private final File mSomatic;
    private final File mCopyNumber;
    private final File structuralVariants;

    private String VCF_FILE_SUFFIX1 = "annotated.vcf";
    private String VCF_FILE_SUFFIX2 = "GoNLv5.vcf";
    private String VCF_FILE_SUFFIX3 = "annotated.vcf.gz";

    RunDirectory(final Path runDirectory)
    {
        mSampleDir = runDirectory;
        mGermline = findGermline();
        mSomatic = findSomatic();
        mCopyNumber = findCopyNumber();
        structuralVariants = findStructuralVariants();
    }

    public File germline() {
        return mGermline;
    }

    public File somatic() {
        return mSomatic;
    }

    public File copyNumber() {
        return mCopyNumber;
    }

    public File structuralVariants() {
        return structuralVariants;
    }

    public Path prefix() {
        return mSampleDir;
    }

    @Nullable
    private File findGermline()
    {
        try
        {
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

    @Nullable
    private File findSomatic()
    {
        try
        {
            try (final Stream<Path> stream = Files.walk(mSampleDir, FileVisitOption.FOLLOW_LINKS))
            {
                return stream.filter(p -> p.toString().endsWith("post_processed_v2.2.vcf.gz"))
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

    @Nullable
    private File findCopyNumber()
    {
        try
        {
            try (final Stream<Path> stream = Files.walk(mSampleDir, FileVisitOption.FOLLOW_LINKS))
            {
                return stream.filter(p -> p.toString().endsWith("purple.gene.cnv")).map(Path::toFile).findFirst().orElse(null);
            }
        }
        catch (final IOException e)
        {
            return null;
        }
    }

    @Nullable
    private File findStructuralVariants()
    {
        try
        {
            try (final Stream<Path> stream = Files.walk(mSampleDir, FileVisitOption.FOLLOW_LINKS))
            {
                return stream.filter(p -> p.toString().endsWith("bpi.vcf")).map(Path::toFile).findFirst().orElse(null);
            }
        }
        catch (final IOException e)
        {
            return null;
        }
    }

    String getPatientID()
    {
        final String[] split = mSampleDir.getFileName().toString().split("_");
        return split[split.length - 1];
    }
}
