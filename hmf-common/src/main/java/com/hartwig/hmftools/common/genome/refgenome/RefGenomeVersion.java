package com.hartwig.hmftools.common.genome.refgenome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum RefGenomeVersion
{
    V37("37", true),
    V38("38", false);

    @NotNull
    private final String mIdentifier;
    private final boolean mIs37;

    // config option
    public static final String REF_GENOME_VERSION = "ref_genome_version";
    public static final String REF_GENOME_VERSION_CFG_DESC = "Ref genome version, 37 or 38";

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeVersion.class);
    private static final String GZIP_EXTENSION = ".gz";

    @NotNull
    public static RefGenomeVersion from(@NotNull final String version)
    {
        if(version.equals(V37.toString()) || version.equals("37") || version.equals("HG37"))
        {
            return V37;
        }
        else if(version.equals(V38.toString()) || version.equals("38") || version.equals("HG38"))
        {
            return V38;
        }

        throw new IllegalArgumentException("Cannot resolve ref genome version: " + version);
    }

    RefGenomeVersion(@NotNull final String identifier, final boolean is37)
    {
        mIdentifier = identifier;
        mIs37 = is37;
    }

    public boolean is37() { return mIs37; }
    public boolean is38 () { return !mIs37; }

    public String identifier() { return mIdentifier; }

    @NotNull
    public String versionedChromosome(@NotNull String chromosome)
    {
        if(this == V38)
        {
            return RefGenomeFunctions.enforceChrPrefix(chromosome);
        }
        else if(this == V37)
        {
            return RefGenomeFunctions.stripChrPrefix(chromosome);
        }
        else
        {
            LOGGER.warn("Unrecognized ref genome version for making chromosome ref genome specific: {}", this);
            return chromosome;
        }
    }

    public String addVersionToFilePath(final String filePath)
    {
        String modifiedFilePath = filePath;
        if(filePath.endsWith(GZIP_EXTENSION))
        {
            modifiedFilePath = filePath.substring(0, filePath.indexOf(GZIP_EXTENSION));
        }

        if(!modifiedFilePath.contains("."))
        {
            throw new IllegalStateException("Cannot include ref genome version in file path that has no proper extension: " + filePath);
        }

        int extensionStart = modifiedFilePath.lastIndexOf(".");
        String versionedFilePath =
                modifiedFilePath.substring(0, extensionStart) + "." + this.mIdentifier + modifiedFilePath.substring(extensionStart);

        if(filePath.endsWith(GZIP_EXTENSION))
        {
            versionedFilePath = versionedFilePath + GZIP_EXTENSION;
        }

        return versionedFilePath;
    }
}
