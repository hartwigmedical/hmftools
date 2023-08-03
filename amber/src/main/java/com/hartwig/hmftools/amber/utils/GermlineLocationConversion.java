package com.hartwig.hmftools.amber.utils;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberUtils.logVersion;
import static com.hartwig.hmftools.common.amber.AmberSiteFactory.FLD_ALT;
import static com.hartwig.hmftools.common.amber.AmberSiteFactory.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.amber.AmberSiteFactory.FLD_POSITION;
import static com.hartwig.hmftools.common.amber.AmberSiteFactory.FLD_REF;
import static com.hartwig.hmftools.common.amber.AmberSiteFactory.FLD_SNP_CHECK;
import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSiteFactory;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class GermlineLocationConversion
{
    private static final String INPUT_GERMLINE_HET_FILE = "input_germline_het_file";
    private static final String OUTPUT_GERMLINE_HET_FILE = "output_germline_het_file";
    private static final String SNP_CHECKS_FILE = "snp_checks_file";
    private static final String DESTINATION_REF_GEN_VERSION = "dest_ref_genome_version";

    private final String mInputFile;
    private final String mOutputFile;
    private final String mSnpCheckFile;
    private final GenomeLiftoverCache mGenomeLiftoverCache;
    private final RefGenomeVersion mDestRefGenVersion;

    public GermlineLocationConversion(final ConfigBuilder configBuilder)
    {
        mInputFile = configBuilder.getValue(INPUT_GERMLINE_HET_FILE);
        mOutputFile = configBuilder.getValue(OUTPUT_GERMLINE_HET_FILE);
        mSnpCheckFile = configBuilder.getValue(SNP_CHECKS_FILE);
        mGenomeLiftoverCache = new GenomeLiftoverCache(true);
        mDestRefGenVersion = RefGenomeVersion.from(configBuilder.getValue(DESTINATION_REF_GEN_VERSION));
    }

    public void run()
    {
        try
        {
            ListMultimap<Chromosome,AmberSite> chrSites = AmberSiteFactory.sites(mInputFile);

            Map<String,Set<Integer>> chrSnpChecks = Maps.newHashMap();

            if(mSnpCheckFile != null)
            {
                ListMultimap<Chromosome,AmberSite> existingSnpCheckSites = AmberSiteFactory.sites(mSnpCheckFile);

                for(HumanChromosome chromosome : HumanChromosome.values())
                {
                    String chrStr = mDestRefGenVersion.versionedChromosome(chromosome.toString());
                    Set<Integer> snpCheckPositions = Sets.newHashSet();
                    chrSnpChecks.put(chrStr, snpCheckPositions);

                    if(!existingSnpCheckSites.containsKey(chromosome))
                        continue;

                    existingSnpCheckSites.get(chromosome).stream()
                            .filter(x -> x.snpCheck()).forEach(x -> snpCheckPositions.add(x.position()));
                }

                AMB_LOGGER.info("applying {} existing SNP-check sites from {}",
                        chrSnpChecks.values().stream().mapToInt(x -> x.size()).sum(), mSnpCheckFile);

            }

            BufferedWriter writer = createBufferedWriter(mOutputFile);

            StringJoiner header = new StringJoiner(TSV_DELIM);
            header.add(FLD_CHROMOSOME).add(FLD_POSITION).add(FLD_REF).add(FLD_ALT).add(FLD_SNP_CHECK);
            writer.write(header.toString());
            writer.newLine();

            int unmappedPositions = 0;
            int totalSites = chrSites.size();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                Set<Integer> snpCheckPositions = chrSnpChecks.get(mDestRefGenVersion.versionedChromosome(chromosome.toString()));

                for(AmberSite site : chrSites.get(chromosome))
                {
                    if(!writeVariant(writer, site, snpCheckPositions))
                        ++unmappedPositions;
                }
            }

            writer.close();

            AMB_LOGGER.info("Amber site conversion complete - mapped {}, unmapped {}",
                    totalSites - unmappedPositions, unmappedPositions);
        }
        catch(IOException e)
        {
            e.printStackTrace();
            AMB_LOGGER.error("failed to load Amber sites file: {}");
            System.exit(1);
        }
    }

    private boolean writeVariant(final BufferedWriter writer, final AmberSite site, final Set<Integer> snpCheckPositions) throws IOException
    {
        int convertedPos = mGenomeLiftoverCache.convertPosition(site.chromosome(), site.position(), mDestRefGenVersion);

        if(convertedPos == UNMAPPED_POSITION)
        {
            AMB_LOGGER.trace("unmapped site({}:{} {}>{})", site.chromosome(), site.position(), site.ref(), site.alt());
            return false;
        }

        boolean snpCheckSite = (snpCheckPositions != null && snpCheckPositions.contains(convertedPos)) || site.snpCheck();

        AMB_LOGGER.trace("converted site({}:{} {}>{}) new position({})",
                site.chromosome(), site.position(), site.ref(), site.alt(), convertedPos);

        String destChr = mDestRefGenVersion.versionedChromosome(site.chromosome());

        StringJoiner data = new StringJoiner(TSV_DELIM);
        data.add(destChr);
        data.add(String.valueOf(convertedPos));
        data.add(site.ref());
        data.add(site.alt());
        data.add(String.valueOf(snpCheckSite));

        writer.write(data.toString());
        writer.newLine();
        return true;
    }

    public static void main(@NotNull final String[] args)
    {
        logVersion();

        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addPath(INPUT_GERMLINE_HET_FILE, true, "Input germline locations file");
        configBuilder.addPath(SNP_CHECKS_FILE, false, "Input germline locations file");
        configBuilder.addConfigItem(OUTPUT_GERMLINE_HET_FILE, true, "Output germline locations file");
        configBuilder.addConfigItem(DESTINATION_REF_GEN_VERSION, true, "Ref genome version to convert to V37 or 38)");
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);
        setLogLevel(configBuilder);

        GermlineLocationConversion application = new GermlineLocationConversion(configBuilder);
        application.run();
    }
}
