package com.hartwig.hmftools.amber.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSitesFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class AmberSitesWriter
{
    private static final String INPUT_GERMLINE_HET_FILE = "input_sites_file";
    private static final String SNP_CHECKS_FILE = "snp_check_sites_file";
    private static final String SOURCE_REF_GEN_VERSION = "source_ref_genome_version";

    private final String mInputFile;
    private final String mOutputDir;
    private final String mSnpCheckFile;
    private final GenomeLiftoverCache mGenomeLiftoverCache;
    private final RefGenomeVersion mSourceRefGenVersion;

    public AmberSitesWriter(final ConfigBuilder configBuilder)
    {
        mInputFile = configBuilder.getValue(INPUT_GERMLINE_HET_FILE);
        mOutputDir = parseOutputDir(configBuilder);
        mSnpCheckFile = configBuilder.getValue(SNP_CHECKS_FILE);
        mSourceRefGenVersion = RefGenomeVersion.from(configBuilder.getValue(SOURCE_REF_GEN_VERSION));

        mGenomeLiftoverCache = new GenomeLiftoverCache(true);
    }

    public static String amberSitesFilename(final RefGenomeVersion version)
    {
        return "AmberGermlineSites." + version.identifier() + TSV_EXTENSION + ".gz";
    }

    public void run()
    {
        AMB_LOGGER.info("generating Amber sites");

        try
        {
            ListMultimap<Chromosome,AmberSite> amberSitesMapMM = AmberSitesFile.sites(mInputFile);

            Map<Chromosome,List<AmberSite>> amberSitesMap = Maps.newHashMap();

            for(Chromosome chromosome : amberSitesMapMM.keySet())
            {
                List<AmberSite> amberSites = amberSitesMapMM.get(chromosome);

                if(amberSites == null)
                    continue;

                // clear any previous SnpCheck flags
                amberSites.forEach(x -> x.setSnpCheck(false));
                amberSitesMap.put(chromosome, amberSitesMapMM.get(chromosome));
            }

            ListMultimap<Chromosome,AmberSite> snpCheckSitesMap = AmberSitesFile.loadVcf(mSnpCheckFile);

            AMB_LOGGER.info("loaded {} SnpCheck sites from {}", snpCheckSitesMap.size(), mSnpCheckFile);

            // check that all required SNP check sites are present
            if(!validateSnpCheckSites(amberSitesMap, snpCheckSitesMap))
            {
                System.exit(1);
            }

            AMB_LOGGER.info("validated all SnpCheck sites");

            Map<RefGenomeVersion,BufferedWriter> writers = Maps.newHashMap();

            for(RefGenomeVersion version : RefGenomeVersion.values())
            {
                String amberSitesFile = mOutputDir + amberSitesFilename(version);

                BufferedWriter writer = createBufferedWriter(amberSitesFile);
                writers.put(version, writer);

                writer.write(AmberSitesFile.header());
                writer.newLine();
            }

            for(RefGenomeVersion version : RefGenomeVersion.values())
            {
                BufferedWriter writer = writers.get(version);
                int writeCount = 0;

                for(HumanChromosome chromosome : HumanChromosome.values())
                {
                    List<AmberSite> amberSites = amberSitesMap.get(chromosome);

                    if(amberSites == null)
                        continue;

                    // ensure output is sorted
                    Collections.sort(amberSites);

                    for(AmberSite amberSite : amberSites)
                    {
                        writeVariant(writer, version, amberSite);

                        ++writeCount;

                        if((writeCount % 1000000) == 0)
                        {
                            AMB_LOGGER.debug("version {}: {} variants written, current site({}:{})",
                                    version, writeCount, amberSite.Chromosome, amberSite.Position);
                        }
                    }
                }
            }

            writers.values().forEach(x -> closeBufferedWriter(x));

            AMB_LOGGER.info("Amber site files generated");
        }
        catch(IOException e)
        {
            e.printStackTrace();
            AMB_LOGGER.error("failed to load Amber sites file: {}");
            System.exit(1);
        }

        AMB_LOGGER.info("Amber site generation complete");
    }

    private boolean validateSnpCheckSites(
            final Map<Chromosome,List<AmberSite>> amberSitesMap, final ListMultimap<Chromosome,AmberSite> snpCheckSitesMap)
    {
        for(Chromosome chromosome : snpCheckSitesMap.keySet())
        {
            List<AmberSite> amberSites = amberSitesMap.get(chromosome);

            if(amberSites == null)
            {
                AMB_LOGGER.error("SnpCheck failed: missing Amber sites for chromosome({})", chromosome);
                return false;
            }

            for(AmberSite snpCheckSite : snpCheckSitesMap.get(chromosome))
            {
                AmberSite matchedSite = amberSites.stream().filter(x -> x.matches(snpCheckSite)).findFirst().orElse(null);

                if(matchedSite != null)
                {
                    AMB_LOGGER.debug("matched SnpCheck site({})", snpCheckSite);
                    matchedSite.setSnpCheck(true);
                }
                else
                {
                    AMB_LOGGER.warn("adding missing SnpCheck site({})", snpCheckSite);
                    amberSites.add(snpCheckSite);
                }
            }
        }

        return true;
    }

    private void writeVariant(final BufferedWriter writer, final RefGenomeVersion version, final AmberSite site) throws IOException
    {
        int position = site.Position;

        if(version != mSourceRefGenVersion)
        {
            position = mGenomeLiftoverCache.convertPosition(site.Chromosome, site.Position, version);

            if(position == UNMAPPED_POSITION)
            {
                AMB_LOGGER.warn("unmapped site({}:{} {}>{})", site.Chromosome, site.Position, site.Ref, site.Alt);
                return;
            }
        }

        String destChr = version.versionedChromosome(site.Chromosome);

        writer.write(format("%s\t%d\t%s\t%s\t%s", destChr, position, site.Ref, site.Alt, site.snpCheck()));
        writer.newLine();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addPath(INPUT_GERMLINE_HET_FILE, true, "Input germline locations file");
        configBuilder.addPath(SNP_CHECKS_FILE, true, "Input germline locations file");
        configBuilder.addConfigItem(SOURCE_REF_GEN_VERSION, true, "Ref genome version to convert to V37 or 38)");
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        AmberSitesWriter application = new AmberSitesWriter(configBuilder);
        application.run();
    }
}
