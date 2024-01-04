package com.hartwig.hmftools.amber.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.utils.Mappability.MAPPABILITY_BED;
import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSitesFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GcCalcs;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineLocationConversion
{
    private static final String INPUT_GERMLINE_HET_FILE = "input_germline_het_file";
    private static final String OUTPUT_GERMLINE_HET_FILE = "output_germline_het_file";
    private static final String SNP_CHECKS_FILE = "snp_checks_file";
    private static final String DESTINATION_REF_GEN_VERSION = "dest_ref_genome_version";

    private static final int GC_SEQUENCE_LENGTH = 120;

    private final String mInputFile;
    private final String mOutputFile;
    private final String mSnpCheckFile;
    private final GenomeLiftoverCache mGenomeLiftoverCache;
    private final RefGenomeVersion mDestRefGenVersion;
    private final Mappability mMappability;
    private final RefGenomeSource mRefGenome;

    public GermlineLocationConversion(final ConfigBuilder configBuilder)
    {
        mInputFile = configBuilder.getValue(INPUT_GERMLINE_HET_FILE);
        mOutputFile = configBuilder.getValue(OUTPUT_GERMLINE_HET_FILE);
        mSnpCheckFile = configBuilder.getValue(SNP_CHECKS_FILE);
        mDestRefGenVersion = RefGenomeVersion.from(configBuilder.getValue(DESTINATION_REF_GEN_VERSION));
        mGenomeLiftoverCache = new GenomeLiftoverCache(true);
        mMappability = new Mappability(configBuilder);

        String refGenomeFile = configBuilder.getValue(REF_GENOME);
        mRefGenome = loadRefGenome(refGenomeFile);
    }

    public void run()
    {
        try
        {
            ListMultimap<Chromosome, AmberReferenceSite> chrSites = loadReferenceAmberSites(mInputFile);

            Map<String,Set<Integer>> chrSnpChecks = Maps.newHashMap();

            if(mSnpCheckFile != null)
            {
                ListMultimap<Chromosome, AmberSite> existingSnpCheckSites = AmberSitesFile.sites(mSnpCheckFile);

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

            writer.write(format("%s\tGnomadFreq\tMappability\tGcRatio", AmberSitesFile.header()));
            writer.newLine();

            int unmappedPositions = 0;
            int totalSites = chrSites.size();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                Set<Integer> snpCheckPositions = chrSnpChecks.get(mDestRefGenVersion.versionedChromosome(chromosome.toString()));

                for(AmberReferenceSite site : chrSites.get(chromosome))
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

    private boolean writeVariant(final BufferedWriter writer, final AmberReferenceSite site, final Set<Integer> snpCheckPositions) throws IOException
    {
        int convertedPos = mGenomeLiftoverCache.convertPosition(site.Chromosome, site.Position, mDestRefGenVersion);

        if(convertedPos == UNMAPPED_POSITION)
        {
            AMB_LOGGER.trace("unmapped site({}:{} {}>{})", site.Chromosome, site.Position, site.Ref, site.Alt);
            return false;
        }

        boolean snpCheckSite = snpCheckPositions != null && snpCheckPositions.contains(convertedPos);
        String destChr = mDestRefGenVersion.versionedChromosome(site.Chromosome);

        double mappability = mMappability.getMappability(destChr, convertedPos);

        double gcRatio = calcsGcContent(destChr, convertedPos, site.Ref, site.Alt);

        AMB_LOGGER.trace("converted site({}:{} {}>{}) new position({})",
                site.Chromosome, site.Position, site.Ref, site.Alt, convertedPos);

        StringJoiner data = new StringJoiner(TSV_DELIM);
        data.add(destChr);
        data.add(String.valueOf(convertedPos));
        data.add(site.Ref);
        data.add(site.Alt);
        data.add(String.valueOf(snpCheckSite));
        data.add(String.valueOf(site.GnomadFrequency));
        data.add(String.valueOf(mappability));
        data.add(format("%.3f", gcRatio));

        writer.write(data.toString());
        writer.newLine();
        return true;
    }

    private double calcsGcContent(final String chromosome, final int position, final String ref, final String alt)
    {
        int altLength = alt.length();
        int refLength = ref.length();
        int startLength = GC_SEQUENCE_LENGTH / 2 - altLength / 2;
        int startPos = position - startLength;

        String basesStart = mRefGenome.getBaseString(chromosome, startPos, position - 1);
        int endBaseLength = GC_SEQUENCE_LENGTH - basesStart.length() - altLength;

        int postPosition = position + refLength;
        String basesEnd = mRefGenome.getBaseString(chromosome, postPosition, postPosition + endBaseLength - 1);

        String sequence = basesStart + alt + basesEnd;
        return GcCalcs.calcGcPercent(sequence);
    }

    private class AmberReferenceSite
    {
        public final String Chromosome;
        public final int Position;
        public final String Ref;
        public final String Alt;
        public final double GnomadFrequency;

        public AmberReferenceSite(final String chromosome, final int position, final String ref, final String alt, final double gnomadFrequency)
        {
            Chromosome = chromosome;
            Position = position;
            Ref = ref;
            Alt = alt;
            GnomadFrequency = gnomadFrequency;
        }
    }

    private static final String GNOMAD_AF = "AF";

    private ListMultimap<Chromosome, AmberReferenceSite> loadReferenceAmberSites(final String vcfFile) throws IOException
    {
        final ListMultimap<Chromosome, AmberReferenceSite> result = ArrayListMultimap.create();

        VcfFileReader reader = new VcfFileReader(vcfFile);

        if(!reader.fileValid())
            throw new IOException("invalid Amber sites file");

        for(VariantContext variant : reader.iterator())
        {
            if(variant.isFiltered())
                continue;

            if(!HumanChromosome.contains(variant.getContig()))
                continue;

            HumanChromosome chromosome = HumanChromosome.fromString(variant.getContig());

            result.put(
                    chromosome,
                    new AmberReferenceSite(variant.getContig(), variant.getStart(),
                            variant.getReference().getBaseString(), variant.getAlternateAllele(0).getBaseString(),
                            variant.getAttributeAsDouble(GNOMAD_AF, 0)));
        }

        AMB_LOGGER.info("loaded {} Amber germline sites", result.size());
        return result;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addPath(INPUT_GERMLINE_HET_FILE, true, "Input germline locations file");
        configBuilder.addPath(SNP_CHECKS_FILE, false, "Input germline locations file");
        configBuilder.addConfigItem(OUTPUT_GERMLINE_HET_FILE, true, "Output germline locations file");
        configBuilder.addConfigItem(DESTINATION_REF_GEN_VERSION, true, "Ref genome version to convert to V37 or 38)");
        configBuilder.addPath(MAPPABILITY_BED, false, "Mappability file");
        configBuilder.addPath(REF_GENOME, false, REF_GENOME_CFG_DESC);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GermlineLocationConversion application = new GermlineLocationConversion(configBuilder);
        application.run();
    }
}
