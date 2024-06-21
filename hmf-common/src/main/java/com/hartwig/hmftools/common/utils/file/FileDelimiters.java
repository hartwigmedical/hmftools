package com.hartwig.hmftools.common.utils.file;

import java.nio.file.Files;
import java.nio.file.Paths;

public final class FileDelimiters
{
    public static final String TSV_DELIM = "\t";
    public static final String CSV_DELIM = ",";
    public static final String ITEM_DELIM = ";";

    public static final String CSV_EXTENSION = ".csv";
    public static final String TSV_EXTENSION = ".tsv";
    public static final String TSV_ZIP_EXTENSION = ".tsv.gz";
    public static final String ZIP_EXTENSION = ".gz";
    public static final String VCF_ZIP_EXTENSION = ".vcf.gz";
    public static final String BAM_EXTENSION = ".bam";
    public static final String BAM_INDEX_EXTENSION = ".bai";

    public static String inferFileDelimiter(final String filename)
    {
        return filename.endsWith(CSV_EXTENSION) ? CSV_DELIM : TSV_DELIM;
    }

    public static String inferHeaderDelimiter(final String header)
    {
        return header.contains(CSV_DELIM) ? CSV_DELIM : TSV_DELIM;
    }

    public static String switchFilenameExtension(final String filename, final String requiredExtension)
    {
        if(filename.endsWith(requiredExtension))
            return filename;

        String otherExtension = requiredExtension.equals(CSV_EXTENSION) ? TSV_EXTENSION : CSV_EXTENSION;

        if(filename.endsWith(otherExtension))
            return filename.substring(0, filename.length() - otherExtension.length()) + requiredExtension;
        else
            return filename;
    }

    public static String checkFileExtensionRename(final String filename)
    {
        if(Files.exists(Paths.get(filename)))
            return filename;

        String switchedFilename = filename.endsWith(CSV_EXTENSION) ?
                switchFilenameExtension(filename, TSV_EXTENSION) : switchFilenameExtension(filename, CSV_EXTENSION);

        if(Files.exists(Paths.get(switchedFilename)))
            return switchedFilename;

        return filename;
    }
}
