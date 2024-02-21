package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.CommonUtils.logVersion;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.File;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadCanonicalTranscripts
{
    private static final String ENSEMBL_DATA_CACHE_ROOT_DIR = "ensembl_data_root_dir";

    public static void main(@NotNull String[] args) throws ParseException, SQLException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addPath(ENSEMBL_DATA_CACHE_ROOT_DIR, true, "HMF common resources Ensembl data cache root directory");

        configBuilder.checkAndParseCommandLine(args);

        try(DatabaseAccess dbAccess = databaseAccess(configBuilder))
        {
            final String ensemblRootDir = configBuilder.getValue(ENSEMBL_DATA_CACHE_ROOT_DIR);

            LOGGER.info("Persisting transcripts to database");
            loadCanonicalTranscripts(dbAccess, ensemblRootDir, RefGenomeVersion.V37);
            loadCanonicalTranscripts(dbAccess, ensemblRootDir, RefGenomeVersion.V38);

            LOGGER.info("Complete");
        }
        catch(Exception e)
        {
            LOGGER.error("Failed to load canonical transcripts", e);
            System.exit(1);
        }
    }

    private static void loadCanonicalTranscripts(final DatabaseAccess dbAccess, final String ensemblRootDir, final RefGenomeVersion refGenomeVersion)
    {
        String ensemblDir = ensemblRootDir + File.separator + refGenomeVersion.identifier();
        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);
        ensemblDataCache.setRequiredData(true, false, false, true);
        ensemblDataCache.load(false);

        List<GeneData> geneDataList = new ArrayList<>();
        List<TranscriptData> transDataList = new ArrayList<>();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            for(GeneData geneData :  ensemblDataCache.getChrGeneDataMap().get(refGenomeVersion.versionedChromosome(chromosome.toString())))
            {
                TranscriptData transData = ensemblDataCache.getTranscriptData(geneData.GeneId, "");
                geneDataList.add(geneData);
                transDataList.add(transData);
            }
        }

        LOGGER.info("Loading {} canonical transcripts to database for refGenome({})", geneDataList.size(), refGenomeVersion);

        String refGenomeStr = refGenomeVersion.is37() ? "GRCh37" : "GRCh38";
        dbAccess.writeCanonicalTranscripts(refGenomeStr, geneDataList, transDataList);
    }

}
