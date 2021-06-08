package com.hartwig.hmftools.common.genome.refgenome;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public enum RefGenomeVersion
{
    V37("37", true),
    V38("38", false),
    HG19("37", true); // included to distinguish from GRCh37 since has has the 'chr' prefix

    @NotNull
    private final String mIdentifier;
    private final boolean mIs37;

    // config option
    public static final String REF_GENOME_VERSION = "ref_genome_version";

    private static final Logger LOGGER = LogManager.getLogger(RefGenomeVersion.class);
    private static final String GZIP_EXTENSION = ".gz";

    @NotNull
    public static RefGenomeVersion from(@NotNull final String version)
    {
        // TODO Remove handling of RG per 1st of july 2021
        if (version.startsWith("RG"))
        {
            LOGGER.warn("Avoid using ref genome versions starting with RG: {}", version);
        }

        if (version.equals(V37.toString()) || version.equals("RG_37") || version.equals("37") || version.equals("HG37"))
        {
            return V37;
        }
        else if (version.equals(V38.toString()) || version.equals("RG_38") || version.equals("38") || version.equals("HG38"))
        {
            return V38;
        }
        else if (version.equals(HG19.toString()) || version.equals("RG_19") || version.equals("19") || version.equals("HG19"))
        {
            return HG19;
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

    @NotNull
    public String versionedChromosome(@NotNull String chromosome)
    {
        if (this == V38 || this == HG19)
        {
            return RefGenomeFunctions.enforceChrPrefix(chromosome);
        }
        else if (this == V37)
        {
            return RefGenomeFunctions.stripChrPrefix(chromosome);
        }
        else
        {
            LOGGER.warn("Unrecognized ref genome version for making chromosome ref genome specific: {}", this);
            return chromosome;
        }
    }

    @NotNull
    public String addVersionToFilePath(@NotNull String filePath)
    {
        String modifiedFilePath = filePath;
        if (filePath.endsWith(GZIP_EXTENSION))
        {
            modifiedFilePath = filePath.substring(0, filePath.indexOf(GZIP_EXTENSION));
        }

        if (!modifiedFilePath.contains("."))
        {
            throw new IllegalStateException("Cannot include ref genome version in file path that has no proper extension: " + filePath);
        }

        int extensionStart = modifiedFilePath.lastIndexOf(".");
        String versionedFilePath =
                modifiedFilePath.substring(0, extensionStart) + "." + this.mIdentifier + modifiedFilePath.substring(extensionStart);

        if (filePath.endsWith(GZIP_EXTENSION))
        {
            versionedFilePath = versionedFilePath + GZIP_EXTENSION;
        }

        return versionedFilePath;
    }
}
