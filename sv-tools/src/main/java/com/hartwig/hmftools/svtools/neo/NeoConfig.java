package com.hartwig.hmftools.svtools.neo;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.svtools.common.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.svtools.neo.NeoEpitopeAnnotator.IM_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

public class NeoConfig
{
    public final List<String> SampleIds;

    public final RefGenomeInterface RefGenome;

    public final List<String> RestrictedGeneIds;

    public final int RequiredAminoAcids;
    public final boolean WriteTransData;
    public final String OutputDir;

    public static final String SAMPLE = "sample";

    public static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";
    public static final String SV_FUSION_FILE = "sv_fusion_file";
    public static final String GENE_ID_FILE = "gene_id_file";
    public static final String REF_GENOME = "ref_genome";
    public static final String REQ_AMINO_ACIDS = "req_amino_acids";
    public static final String WRITE_TRANS_DATA = "write_trans_data";

    public static final int DEFAULT_AMINO_ACID_REF_COUNT = 18;

    public NeoConfig(final CommandLine cmd)
    {
        SampleIds = Lists.newArrayList();

        final String sampleIdConfig = cmd.getOptionValue(SAMPLE);

        if(sampleIdConfig.contains(".csv"))
        {
            loadSampleIdsFile(sampleIdConfig);
        }
        else if(sampleIdConfig.contains(";"))
        {
            Arrays.stream(sampleIdConfig.split(";")).forEach(x -> SampleIds.add(x));
        }
        else
        {
            SampleIds.add(sampleIdConfig);
        }

        final String refGenomeFilename = cmd.getOptionValue(REF_GENOME);
        RefGenome = loadRefGenome(refGenomeFilename);

        RestrictedGeneIds = Lists.newArrayList();
        OutputDir = parseOutputDir(cmd);

        RequiredAminoAcids = Integer.parseInt(cmd.getOptionValue(REQ_AMINO_ACIDS, String.valueOf(DEFAULT_AMINO_ACID_REF_COUNT)));
        WriteTransData = cmd.hasOption(WRITE_TRANS_DATA);

        if(cmd.hasOption(GENE_ID_FILE))
        {
            loadGeneIdsFile(cmd.getOptionValue(GENE_ID_FILE), RestrictedGeneIds);
        }
    }

    public NeoConfig(
            final List<String> sampleIds, final RefGenomeInterface refGenome, final List<String> restrictedGeneIds,
            final int requiredAminoAcids)
    {
        SampleIds = sampleIds;
        RefGenome = refGenome;
        RestrictedGeneIds = restrictedGeneIds;
        RequiredAminoAcids = requiredAminoAcids;
        OutputDir = "";
        WriteTransData = false;
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE, true, "Sample - Id(s) separated by ';' or CSV file");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Ensembl data cache directory");
        options.addOption(GENE_ID_FILE, true, "Restrict to specific genes");
        options.addOption(REF_GENOME, true, "Ref genome");
        options.addOption(SV_FUSION_FILE, true, "SV fusion file (single sample or cohort)");
        options.addOption(WRITE_TRANS_DATA, false, "Write transcript data for each neo-epitope");
        options.addOption(REQ_AMINO_ACIDS, true, "Number of amino acids in neo-epitopes (default: 18)");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        addDatabaseCmdLineArgs(options);
    }

    private void loadGeneIdsFile(final String filename, final List<String> geneIdList)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            IM_LOGGER.warn("invalid gene ID file({})", filename);
            return;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            if (fileContents.get(0).contains("GeneId"))
            {
                // check for header row
                fileContents.remove(0);
            }

            geneIdList.addAll(fileContents.stream()
                    .filter(x -> !x.isEmpty())
                    .filter(x -> !x.contains("GeneId"))
                    .filter(x -> !x.startsWith("#"))
                    .map(x -> x.split(",")[0]).collect(Collectors.toList()));
        }
        catch (IOException e)
        {
            IM_LOGGER.warn("failed to load gene ID file({}): {}", filename, e.toString());
        }
    }

    private void loadSampleIdsFile(final String filename)
    {
        if (!Files.exists(Paths.get(filename)))
        {
            IM_LOGGER.warn("invalid sampleId file({})", filename);
            return;
        }

        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            if(fileContents.isEmpty())
                return;

            if (fileContents.get(0).contains("SampleId"))
                fileContents.remove(0);

            SampleIds.addAll(fileContents);
        }
        catch (IOException e)
        {
            IM_LOGGER.warn("failed to load sampleId file({}): {}", filename, e.toString());
        }
    }


}
