package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdFile;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jetbrains.annotations.NotNull;

public class ImpactConfig
{
    public final List<String> SampleIds;

    public final String VcfFile;
    public final RefGenomeVersion RefGenVersion;
    public final DatabaseAccess DbAccess;

    public final String OutputDir;
    public final boolean WriteCsv;
    public final boolean WriteTranscriptCsv;

    private static final String SAMPLE = "sample";
    public static final String SAMPLE_ID_FILE = "sample_id_file";
    private static final String VCF_FILE = "vcf_file";
    private static final String WRITE_CSV = "write_csv";
    private static final String WRITE_TRANSCRIPT_CSV = "write_transcript_csv";

    public static final String REF_GENOME = "ref_genome";

    public ImpactConfig(final CommandLine cmd)
    {
        SampleIds = Lists.newArrayList();

        if(cmd.hasOption(SAMPLE))
        {
            SampleIds.add(cmd.getOptionValue(SAMPLE));
            VcfFile = cmd.getOptionValue(VCF_FILE);
            DbAccess = null;
        }
        else
        {
            SampleIds.addAll(loadSampleIdFile(cmd.getOptionValue(SAMPLE_ID_FILE)));
            DbAccess = createDatabaseAccess(cmd);
            VcfFile = null;
        }

        RefGenVersion = cmd.hasOption(REF_GENOME_VERSION) ? RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION)) : V37;


        WriteCsv = cmd.hasOption(WRITE_CSV);
        WriteTranscriptCsv = cmd.hasOption(WRITE_TRANSCRIPT_CSV);

        OutputDir = parseOutputDir(cmd);
    }

    public boolean singleSample() { return SampleIds.size() == 1; }

    public boolean singleSampleValid()
    {
        if(SampleIds.isEmpty())
            return false;

        if(VcfFile == null)
            return false;

        return true;
    }

    public boolean multiSampleValid()
    {
        if(SampleIds.isEmpty())
            return false;

        if(DbAccess == null)
            return false;

        return true;
    }

    @NotNull
    public static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Name of sample");
        options.addOption(SAMPLE_ID_FILE, true, "Sample ID file");
        options.addOption(VCF_FILE, true, "VCF input file");

        options.addOption(REF_GENOME, true, "Path to ref genome fasta file");
        options.addOption(REF_GENOME_VERSION, true, "Ref genome version: V37(default) or V38");
        options.addOption(ENSEMBL_DATA_DIR, true, "Name of sample");
        options.addOption(DRIVER_GENE_PANEL_OPTION, true, DRIVER_GENE_PANEL_OPTION_DESC);

        options.addOption(WRITE_CSV, false, "Write variant impacts to CSV");
        options.addOption(WRITE_TRANSCRIPT_CSV, false, "Write variant impacts per transcript to CSV");

        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}
