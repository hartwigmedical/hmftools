package com.hartwig.hmftools.common.pipeline;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getValue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

public class PipelineToolDirectoriesFile
{
    public static final String AMBER_DIR = "amberDir";
    public static final String CHORD_DIR = "chordDir";
    public static final String CIDER_DIR = "ciderDir";
    public static final String COBALT_DIR = "cobaltDir";
    public static final String CUPPA_DIR = "cuppaDir";
    public static final String ESVEE_DIR = "esveeDir";
    public static final String ISOFOX_DIR = "isofoxDir";
    public static final String GERMLINE_FLAGSTAT_DIR = "germlineFlagstatDir";
    public static final String GERMLINE_METRICS_DIR = "germlineMetricsDir";
    public static final String LILAC_DIR = "lilacDir";
    public static final String LINX_GERMLINE_DIR = "linxGermlineDir";
    public static final String LINX_SOMATIC_DIR = "linxSomaticDir";
    public static final String ORANGE_DIR = "orangeDir";
    public static final String PAVE_GERMLINE_DIR = "paveGermlineDir";
    public static final String PAVE_SOMATIC_DIR = "paveSomaticDir";
    public static final String PEACH_DIR = "peachDir";
    public static final String PURPLE_DIR = "purpleDir";
    public static final String SAGE_GERMLINE_DIR = "sageGermlineDir";
    public static final String SAGE_SOMATIC_DIR = "sageSomaticDir";
    public static final String SIGS_DIR = "sigsDir";
    public static final String SNP_GENOTYPE_DIR = "snpGenotypeDir";
    public static final String TEAL_DIR = "tealDir";
    public static final String TUMOR_FLAGSTAT_DIR = "tumorFlagstatDir";
    public static final String TUMOR_METRICS_DIR = "tumorMetricsDir";
    public static final String V_CHORD_DIR = "vChordDir";
    public static final String VIRUS_BREAKEND_DIR = "virusBreakendDir";
    public static final String VIRUS_INTERPRETER_DIR = "virusInterpreterDir";

    public static final String DEFAULT_DIR = "";


    public static PipelineToolDirectories read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    private static PipelineToolDirectories fromLines(final List<String> lines)
    {
        return new PipelineToolDirectories(
                getValue(lines, AMBER_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, CHORD_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, CIDER_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, COBALT_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, CUPPA_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, ESVEE_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, GERMLINE_FLAGSTAT_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, GERMLINE_METRICS_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, ISOFOX_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, LILAC_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, LINX_GERMLINE_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, LINX_SOMATIC_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, ORANGE_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, PAVE_GERMLINE_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, PAVE_SOMATIC_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, PEACH_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, PURPLE_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, SAGE_GERMLINE_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, SAGE_SOMATIC_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, SIGS_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, SNP_GENOTYPE_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, TEAL_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, TUMOR_FLAGSTAT_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, TUMOR_METRICS_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, V_CHORD_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, VIRUS_BREAKEND_DIR, DEFAULT_DIR, TSV_DELIM),
                getValue(lines, VIRUS_INTERPRETER_DIR, DEFAULT_DIR, TSV_DELIM)
        );
    }
}
