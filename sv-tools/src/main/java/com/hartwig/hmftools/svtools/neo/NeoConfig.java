package com.hartwig.hmftools.svtools.neo;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
import static com.hartwig.hmftools.svtools.neo.NeoEpitopeAnnotator.IM_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

public class NeoConfig
{
    public final List<String> SampleIds;

    public final RefGenomeInterface RefGenome;

    public final List<String> RestrictedGeneIds;

    public static final String SAMPLE = "sample";

    public static final String GENE_TRANSCRIPTS_DIR = "gene_transcripts_dir";
    public static final String COHORT_FUSION_FILE = "cohort_fusion_file";
    public static final String GENE_ID_FILE = "gene_id_file";
    public static final String REF_GENOME = "ref_genome";

    public static final int AMINO_ACID_REF_COUNT = 18;

    public NeoConfig(final CommandLine cmd)
    {
        SampleIds = Lists.newArrayList();

        final String refGenomeFilename = cmd.getOptionValue(REF_GENOME);
        RefGenome = loadRefGenome(refGenomeFilename);

        RestrictedGeneIds = Lists.newArrayList();

        if(cmd.hasOption(GENE_ID_FILE))
        {
            loadGeneIdsFile(cmd.getOptionValue(GENE_ID_FILE), RestrictedGeneIds);
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE, true, "Path to the Linx cohort SVs file");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "External LINE data sample counts");
        options.addOption(GENE_ID_FILE, true, "Path to write results");
        options.addOption(REF_GENOME, false, "Log verbose");
        options.addOption(COHORT_FUSION_FILE, false, "Log verbose");
        options.addOption(OUTPUT_DIR, false, "Log verbose");
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


}
