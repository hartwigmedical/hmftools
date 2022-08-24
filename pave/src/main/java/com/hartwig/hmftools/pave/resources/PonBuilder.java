package com.hartwig.hmftools.pave.resources;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.variant.VariantVcfTags.TIER;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.variant.VariantTier.HOTSPOT;
import static com.hartwig.hmftools.pave.PaveConfig.PON_FILE;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PonAnnotation.PON_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.pave.PonAnnotation;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class PonBuilder
{
    private final List<String> mSampleIds;
    private final String mVcfPath;
    private final String mOutputDir;

    private final int mQualCutoff;
    private final int mMinSamples;
    private final RefGenomeVersion mRefGenomeVersion;
    private final Map<String,List<VariantData>> mChrVariantsMap;

    private final PonAnnotation mExistingPon;

    private int mLastIndex;
    private String mLastChromosome;

    private static final String SAMPLE_ID_FILE = "sample_id_file";
    private static final String VCF_PATH = "vcf_path";
    private static final String QUAL_CUTOFF = "qual_cutoff";
    private static final String MIN_SAMPLES = "min_samples";
    private static final String MANUAL_ENTRIES = "manual_entries";

    public PonBuilder(final CommandLine cmd)
    {
        mSampleIds = loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE));
        mVcfPath = cmd.getOptionValue(VCF_PATH);
        mOutputDir = parseOutputDir(cmd);

        mQualCutoff = Integer.parseInt(cmd.getOptionValue(QUAL_CUTOFF));
        mMinSamples = Integer.parseInt(cmd.getOptionValue(MIN_SAMPLES));

        mRefGenomeVersion = RefGenomeVersion.from(cmd.getOptionValue(REF_GENOME_VERSION));

        mExistingPon = new PonAnnotation(cmd.getOptionValue(PON_FILE), true);

        mChrVariantsMap = Maps.newHashMap();

        mLastChromosome = "";
        mLastIndex = 0;

        if(cmd.hasOption(MANUAL_ENTRIES))
        {
            String[] entries = cmd.getOptionValue(MANUAL_ENTRIES).split(";", -1);

            for(String entry : entries)
            {
                String[] varValues = entry.split(":", 5);
                addVariant(varValues[0], Integer.parseInt(varValues[1]), varValues[2], varValues[3], Integer.parseInt(varValues[4]));
            }
        }
    }

    public void run()
    {
        if(mSampleIds.isEmpty() || mVcfPath == null)
        {
            PV_LOGGER.error("missing sample IDs file or VCF path in config");
            System.exit(1);
        }

        PV_LOGGER.info("generating PON from {} samples, minSamples({}) qualCutoff({})",
                mSampleIds.size(), mMinSamples, mQualCutoff);

        for(String sampleId : mSampleIds)
        {
            loadSampleVcf(sampleId);
        }

        writePon();
    }

    private void loadSampleVcf(final String sampleId)
    {
        String vcfFilename = mVcfPath.replaceAll("\\*", sampleId);

        AbstractFeatureReader<VariantContext, LineIterator> vcfReader = AbstractFeatureReader.getFeatureReader(vcfFilename, new VCFCodec(), false);

        try
        {
            int varCount = 0;

            int lastPosition = 0;
            String lastChromosome = "";
            String lastRef = "";
            String lastAlt = "";

            for(VariantContext variantContext : vcfReader.iterator())
            {
                double qual = variantContext.getPhredScaledQual();

                if(qual < mQualCutoff)
                    continue;

                String tier = variantContext.getAttributeAsString(TIER, "");
                if(tier.equals(HOTSPOT.toString()))
                    continue;

                int position = variantContext.getStart();
                String chromosome = variantContext.getContig();
                String ref = variantContext.getReference().getBaseString();
                String alt = variantContext.getAlternateAlleles().get(0).toString();

                // ignore duplicates (eg with different read-contexts)
                if(chromosome.equals(lastChromosome) && lastPosition == position && lastRef.equals(ref) && lastAlt.equals(alt))
                    continue;

                addVariant(chromosome, position, ref, alt, 1);
                ++varCount;

                lastChromosome = chromosome;
                lastPosition = position;
                lastRef = ref;
                lastAlt = alt;

                if((varCount % 100000) == 0)
                {
                    PV_LOGGER.debug("sample({}) loaded {} variants", sampleId, varCount);
                }
            }

            PV_LOGGER.info("sample({}) read {} variants from vcf({})", sampleId, varCount, vcfFilename);
        }
        catch(IOException e)
        {
            PV_LOGGER.error("error reading vcf file: {}", e.toString());
        }
    }

    private void writePon()
    {
        try
        {
            String fileName = format("%s/somatic_pon_%d_samples.%s.tsv", mOutputDir, mSampleIds.size(), mRefGenomeVersion.identifier());

            PV_LOGGER.info("writing PON file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName);

            StringJoiner sj = new StringJoiner(PON_DELIM);
            sj.add("Chromosome");
            sj.add("Position");
            sj.add("Ref");
            sj.add("Alt");
            sj.add("SampleCount");
            writer.write(sj.toString());
            writer.newLine();

            for(HumanChromosome chromosome : HumanChromosome.values())
            {
                String chrStr = mRefGenomeVersion.versionedChromosome(chromosome.toString());

                List<VariantData> variants = mChrVariantsMap.get(chrStr);

                if(variants.isEmpty())
                    continue;

                for(VariantData variant : variants)
                {
                    if(variant.SampleCount < mMinSamples)
                        continue;

                    if(mExistingPon.hasEntry(variant.Chromosome, variant.Position, variant.Ref, variant.Alt))
                        continue;

                    sj = new StringJoiner(PON_DELIM);
                    sj.add(variant.Chromosome);
                    sj.add(String.valueOf(variant.Position));
                    sj.add(variant.Ref);
                    sj.add(variant.Alt);
                    sj.add(String.valueOf(variant.SampleCount));
                    writer.write(sj.toString());
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to initialise output file: {}", e.toString());
            return;
        }
    }

    private void addVariant(final String chromosome, final int position, final String ref, final String alt, int sampleCount)
    {
        List<VariantData> variants = mChrVariantsMap.get(chromosome);

        if(variants == null)
        {
            variants = Lists.newArrayList();
            mChrVariantsMap.put(chromosome, variants);
        }

        if(!mLastChromosome.equals(chromosome))
        {
            mLastChromosome = chromosome;
            mLastIndex = 0;
        }

        // start from the last inserted index since each VCF is ordered
        int index = mLastIndex;
        while(index < variants.size())
        {
            VariantData variant = variants.get(index);

            if(position > variant.Position)
            {
                ++index;
                continue;
            }

            if(position < variant.Position)
                break;

            if(variant.Ref.equals(ref) && variant.Alt.equals(alt))
            {
                ++variant.SampleCount;
                return;
            }

            ++index;
        }

        VariantData newVariant = new VariantData(chromosome, position, ref, alt);
        newVariant.SampleCount = sampleCount;
        variants.add(index, newVariant);
        mLastIndex = index;
    }

    private class VariantData
    {
        public final String Chromosome;
        public final int Position;
        public final String Ref;
        public final String Alt;
        public int SampleCount;

        public VariantData(final String chromosome, final int position, final String ref, final String alt)
        {
            Chromosome = chromosome;
            Position = position;
            Ref = ref;
            Alt = alt;
            SampleCount = 0;
        }

        public String toString() { return format("var(%s:%d %s>%s) samples(%d)", Chromosome, Position, Ref, Alt, SampleCount); }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(SAMPLE_ID_FILE, true, "Sample IDs file");
        options.addOption(VCF_PATH, true, "VCF path for samples");
        options.addOption(MIN_SAMPLES, true, "Min samples for variant to be included in PON");
        options.addOption(QUAL_CUTOFF, true, "Qual cut-off for variant inclusion");
        options.addOption(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        options.addOption(MANUAL_ENTRIES, true, "Manual PON entries in form Chr:Pos:Ref:Alt separated by ';'");
        options.addOption(PON_FILE, true, "PON entries");

        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PonBuilder ponBuilder = new PonBuilder(cmd);
        ponBuilder.run();

        PV_LOGGER.info("Pave PON building from VCFs complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
