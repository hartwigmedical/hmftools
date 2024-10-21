package com.hartwig.hmftools.pave.resources;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MAP_QUAL_FACTOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.TIER;
import static com.hartwig.hmftools.common.variant.VariantTier.HOTSPOT;
import static com.hartwig.hmftools.pave.PaveConfig.PON_FILE;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.APP_NAME;


import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.pave.annotation.PonAnnotation;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class PonBuilder
{
    private final List<String> mSampleIds;
    private final String mVcfPath;
    private final String mOutputDir;

    private final int mQualCutoff;
    private final int mMqfCutoff;
    private final int mMinSamples;
    private final RefGenomeVersion mRefGenomeVersion;
    private final Map<String,List<VariantData>> mChrVariantsMap;

    private final PonAnnotation mExistingPon;

    private int mLastIndex;
    private String mLastChromosome;

    private static final String VCF_PATH = "vcf_path";
    private static final String QUAL_CUTOFF = "qual_cutoff";
    private static final String MQF_CUTOFF = "mqf_cutoff";
    private static final String MIN_SAMPLES = "min_samples";
    private static final String MANUAL_ENTRIES = "manual_entries";

    private static final int DEFAULT_MIN_SAMPLES = 3;
    private static final int DEFAULT_MIN_MAP_QUAL = -10;
    private static final int DEFAULT_MIN_QUAL = 40;

    public PonBuilder(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder);
        mVcfPath = configBuilder.getValue(VCF_PATH);
        mOutputDir = parseOutputDir(configBuilder);

        mQualCutoff = configBuilder.getInteger(QUAL_CUTOFF);
        mMqfCutoff = configBuilder.getInteger(MQF_CUTOFF);
        mMinSamples = configBuilder.getInteger(MIN_SAMPLES);

        mRefGenomeVersion = RefGenomeVersion.from(configBuilder);

        mExistingPon = new PonAnnotation(configBuilder.getValue(PON_FILE), false);

        mChrVariantsMap = Maps.newHashMap();

        mLastChromosome = "";
        mLastIndex = 0;

        if(configBuilder.hasValue(MANUAL_ENTRIES))
        {
            String[] entries = configBuilder.getValue(MANUAL_ENTRIES).split(ITEM_DELIM, -1);

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

        VcfFileReader vcfReader = new VcfFileReader(vcfFilename);

        int varCount = 0;

        int lastPosition = 0;
        String lastChromosome = "";
        String lastRef = "";
        String lastAlt = "";

            for(VariantContext variantContext : vcfReader.iterator())
            {
                double qual = variantContext.getPhredScaledQual();
                double mqf = variantContext.getAttributeAsDouble(MAP_QUAL_FACTOR, 0);

                if(qual < mQualCutoff || mqf < mMqfCutoff)
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

    private void writePon()
    {
        try
        {
            String fileName = format("%s/somatic_pon_%d_samples.%s.tsv", mOutputDir, mSampleIds.size(), mRefGenomeVersion.identifier());

            PV_LOGGER.info("writing PON file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
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

                    sj = new StringJoiner(TSV_DELIM);
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

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        addSampleIdFile(configBuilder, true);
        configBuilder.addConfigItem(VCF_PATH, true, "VCF path for samples");
        configBuilder.addInteger(MIN_SAMPLES, "Min samples for variant to be included in PON", DEFAULT_MIN_SAMPLES);
        configBuilder.addInteger(QUAL_CUTOFF, "Qual cut-off for variant inclusion", DEFAULT_MIN_QUAL);
        configBuilder.addInteger(MQF_CUTOFF, "MQF cut-off for variant inclusion", DEFAULT_MIN_MAP_QUAL);
        configBuilder.addConfigItem(REF_GENOME_VERSION, true, REF_GENOME_VERSION_CFG_DESC);
        configBuilder.addConfigItem(MANUAL_ENTRIES, false, "Manual PON entries in form Chr:Pos:Ref:Alt separated by ';'");
        configBuilder.addPath(PON_FILE, false, "PON entries");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PonBuilder ponBuilder = new PonBuilder(configBuilder);
        ponBuilder.run();

        PV_LOGGER.info("Pave PON building from VCFs complete");
    }
}
