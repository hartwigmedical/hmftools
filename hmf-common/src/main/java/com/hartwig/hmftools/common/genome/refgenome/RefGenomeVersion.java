package com.hartwig.hmftools.common.genome.refgenome;

import org.jetbrains.annotations.NotNull;

public enum RefGenomeVersion {
    RG_37("37"),
    RG_38("38"),
    RG_19("37"); // included to distinguish from GRCh37 since has has the 'chr' prefix

    // config option
    public static final String REF_GENOME_VERSION = "ref_genome_version";

    public static final String CHR_PREFIX = "chr";

    private static final String GZIP_EXTENSION = ".gz";

    public static boolean is37(final RefGenomeVersion version) {
        return version == RG_37 || version == RG_19;
    }

    @NotNull
    public static RefGenomeVersion from(final String version) {
        if (version.equals(RG_37.toString()) || version.equals("37") || version.equals("HG37")) {
            return RG_37;
        }

        if (version.equals(RG_38.toString()) || version.equals("38") || version.equals("HG38")) {
            return RG_38;
        }

        if (version.equals(RG_19.toString()) || version.equals("19") || version.equals("HG19")) {
            return RG_19;
        }

        return RG_37;
    }

    @NotNull
    public static String refGenomeChromosome(@NotNull final String chromosome, @NotNull RefGenomeVersion version) {
        if ((version == RG_38 || version == RG_19) && !chromosome.contains(CHR_PREFIX)) {
            return CHR_PREFIX + chromosome;
        } else if (version == RG_37) {
            return stripChromosome(chromosome);
        } else {
            return chromosome;
        }
    }

    @NotNull
    public static String stripChromosome(@NotNull final String chromosome) {
        if (chromosome.startsWith(CHR_PREFIX)) {
            return chromosome.substring(CHR_PREFIX.length());
        }

        return chromosome;
    }

    @NotNull
    private final String identifier;

    RefGenomeVersion(@NotNull final String identifier) {
        this.identifier = identifier;
    }

    @NotNull
    public String addVersionToFilePath(@NotNull String filePath) {
        String modifiedFilePath = filePath;
        if (filePath.endsWith(GZIP_EXTENSION)) {
            modifiedFilePath = filePath.substring(0, filePath.indexOf(GZIP_EXTENSION));
        }

        if (!modifiedFilePath.contains(".")) {
            throw new IllegalStateException("Cannot include ref genome version in file path that has no proper extension: " + filePath);
        }

        int extensionStart = modifiedFilePath.lastIndexOf(".");
        String versionedFilePath =
                modifiedFilePath.substring(0, extensionStart) + "." + this.identifier + modifiedFilePath.substring(extensionStart);

        if (filePath.endsWith(GZIP_EXTENSION)) {
            versionedFilePath = versionedFilePath + GZIP_EXTENSION;
        }

        return versionedFilePath;
    }
}
