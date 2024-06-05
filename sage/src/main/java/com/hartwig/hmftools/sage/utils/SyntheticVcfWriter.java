package com.hartwig.hmftools.sage.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.sage.ReferenceData.loadRefGenome;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.vcf.VariantContextFactory.NO_CALL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_BASE_QUAL;
import static com.hartwig.hmftools.sage.vcf.VcfTags.AVG_MAP_QUALITY;
import static com.hartwig.hmftools.sage.vcf.VcfTags.FRAG_STRAND_BIAS;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_IMPROPER_PAIR;
import static com.hartwig.hmftools.sage.vcf.VcfTags.READ_CONTEXT_JITTER;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileReaderUtils;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.vcf.CandidateSerialisation;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.common.VariantTier;
import com.hartwig.hmftools.sage.vcf.VariantVCF;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;

public class SyntheticVcfWriter
{
    private final String mSampleId;
    private final String mOutputVcfFile;
    private final String mInputFilename;

    private final IndexedFastaSequenceFile mRefGenome;
    private final List<VariantData> mVariants;
    private final VariantContextWriter mWriter;

    private static final String OUTPUT_VCF_FILE = "output_vcf_file";
    private static final String INPUT_FILE = "input_file";

    public SyntheticVcfWriter(final ConfigBuilder configBuilder)
    {
        mOutputVcfFile = configBuilder.getValue(OUTPUT_VCF_FILE);
        mInputFilename = configBuilder.getValue(INPUT_FILE);
        mSampleId = configBuilder.getValue(SAMPLE);

        mRefGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));
        mVariants = Lists.newArrayList();

        mWriter = initialiseVcf();
    }

    public void run()
    {
        loadVariants();

        if(mVariants.isEmpty())
            return;

        SG_LOGGER.info("build ref context for {} variants", mVariants.size());

        for(VariantData variant : mVariants)
        {
            buildAndWriteVariant(variant);
        }

        mWriter.close();

        SG_LOGGER.info("compilation complete");
    }

    private class VariantData implements Comparable<VariantData>
    {
        public final String Chromosome;
        public final int Position;
        public final String Ref;
        public final String Alt;

        public VariantData(final String chromosome, final int position, final String ref, final String alt)
        {
            Chromosome = chromosome;
            Position = position;
            Ref = ref;
            Alt = alt;
        }

        @Override
        public int compareTo(@NotNull final VariantData other)
        {
            if(Chromosome.equals(other.Chromosome))
            {
                if(Position != other.Position)
                    return Position < other.Position ? -1 : 1;
                else
                    return 0;
            }

            int rank1 = HumanChromosome.chromosomeRank(Chromosome);
            int rank2 = HumanChromosome.chromosomeRank(other.Chromosome);

            if(rank1 == rank2)
                return 0;

            return rank1 < rank2 ? -1 : 1;
        }

        public String toString()
        {
            return format("%s:%d %s>%s", Chromosome, Position, Ref, Alt);
        }
    }

    private void loadVariants()
    {
        if(mInputFilename.endsWith(TSV_EXTENSION))
        {
            try
            {
                List<String> lines = Files.readAllLines(Paths.get(mInputFilename));

                String header = lines.get(0);
                lines.remove(0);

                Map<String,Integer> fieldIndexMap = FileReaderUtils.createFieldsIndexMap(header, TSV_DELIM);

                int chrIndex = fieldIndexMap.get(FLD_CHROMOSOME);
                int posIndex = fieldIndexMap.get(FLD_POSITION);
                int refIndex = fieldIndexMap.get(FLD_REF);
                int altIndex = fieldIndexMap.get(FLD_ALT);

                for(String line : lines)
                {
                    String[] values = line.split(TSV_DELIM, -1);

                    mVariants.add(new VariantData(values[chrIndex], Integer.parseInt(values[posIndex]), values[refIndex], values[altIndex]));
                }

            }
            catch(Exception e)
            {
                SG_LOGGER.error("failed to read variant output file: {}", e.toString());
                System.exit(1);
            }
        }
        else
        {
            VcfFileReader vcfFileReader = new VcfFileReader(mInputFilename);

            if(!vcfFileReader.fileValid())
            {
                System.exit(1);
                return;
            }

            for(VariantContext variantContext : vcfFileReader.iterator())
            {
                mVariants.add(new VariantData(
                        variantContext.getContig(), variantContext.getStart(),
                        variantContext.getReference().getBaseString(), variantContext.getAlternateAlleles().get(0).toString()));
            }
        }

        // ensure variants are sorted
        Collections.sort(mVariants);
    }

    private VariantContextWriter initialiseVcf()
    {
        SAMSequenceDictionary sequenceDictionary = mRefGenome.getSequenceDictionary();

        VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(mOutputVcfFile)
                .setReferenceDictionary(sequenceDictionary)
                .build();

        VersionInfo version = fromAppName(APP_NAME);
        VCFHeader header = VariantVCF.createHeader(version.version(), List.of(mSampleId), false);

        // CLEAN-UP: something needs to set the tri-nuc, MH and Repeat headers and values if necessary
        // mRefContextEnrichment.appendHeader(header);

        final SAMSequenceDictionary condensedDictionary = new SAMSequenceDictionary();
        for(SAMSequenceRecord sequence : sequenceDictionary.getSequences())
        {
            if(HumanChromosome.contains(sequence.getContig()))
            {
                condensedDictionary.addSequence(sequence);
            }
        }

        header.setSequenceDictionary(condensedDictionary);
        writer.writeHeader(header);

        return writer;
    }

    private VariantReadContext buildReadContext(final VariantData variant)
    {
        VariantReadContextBuilder readContextBuilder = new VariantReadContextBuilder(DEFAULT_FLANK_LENGTH);

        int readFlankLength = 25;
        int refLocationStart = variant.Position - readFlankLength;
        int refLocationEnd = variant.Position + variant.Ref.length() - 1 + readFlankLength;
        String refBasesBefore = mRefGenome.getSubsequenceAt(variant.Chromosome, refLocationStart, variant.Position - 1).getBaseString();

        String refBasesAfter = mRefGenome.getSubsequenceAt(variant.Chromosome, variant.Position, refLocationEnd).getBaseString();

        ChrBaseRegion refRegion = new ChrBaseRegion(variant.Chromosome, refLocationStart, refLocationEnd);
        RefSequence refSequence = new RefSequence(refRegion, mRefGenome);

        // form synthetic read bases
        String readBases = refBasesBefore + variant.Alt + refBasesAfter.substring(variant.Ref.length());

        int readIndex = readFlankLength;

        VariantReadContext readContext = null;

        /* CLEAN=UP
        if(variant.isDelete())
            readContext =  readContextBuilder.createIndelContext(variant.Ref, variant.Position, readIndex, readBases.getBytes(), refSequence);
        else if(variant.isInsert())
            readContext = readContextBuilder.createInsertContext(variant.Alt, variant.Position, readIndex, readBases.getBytes(), refSequence);
        else
            readContext = readContextBuilder.createMNVContext(variant.Position, readIndex, variant.Ref.length(), readBases.getBytes(), refSequence);

        if(readContext.hasIncompleteCore() || readContext.hasIncompleteFlanks())
        {
            SG_LOGGER.error("var({}) incomplete core: core({}) indexed({})",
                    variant, readContext.toString(), readContext.indexedBases().toString());
            return null;
        }
        */

        return readContext;
    }

    private void buildAndWriteVariant(final VariantData variant)
    {
        Genotype genotype = new GenotypeBuilder(mSampleId)
                .DP(0)
                .AD(new int[] { 0, 0 })
                .attribute(READ_CONTEXT_QUALITY, 0)
                .attribute(READ_CONTEXT_COUNT, 0)
                .attribute(READ_CONTEXT_IMPROPER_PAIR, 0)
                .attribute(READ_CONTEXT_JITTER, 0)
                .attribute(FRAG_STRAND_BIAS, 0)
                .attribute(AVG_BASE_QUAL, new int[] { 0, 0 })
                .attribute(AVG_MAP_QUALITY, new int[] { 0, 0 })
                .attribute(VCFConstants.ALLELE_FREQUENCY_KEY, 0)
                .alleles(NO_CALL)
                .make();

        // SimpleVariant simpledVariant = new SimpleVariant(variant.Chromosome, variant.Position, variant.Ref, variant.Alt);

        VariantReadContext readContext = buildReadContext(variant);

        if(readContext == null)
            return;

        Candidate candidate = new Candidate(VariantTier.PANEL, readContext, 0, 0);

        VariantContextBuilder builder = CandidateSerialisation.toContext(candidate)
                .log10PError(0)
                .genotypes(genotype)
                .filters(PASS);

        VariantContext variantContext = builder.make();

        // CLEAN-UP: set those 3 fields or ignore if only for passing SNVs??

        // mRefContextEnrichment.processVariant(variantContext);
        mWriter.add(variantContext);

        SG_LOGGER.debug("variant({}) readContext({} - {} - {}) index({})",
                variant, readContext.leftFlankStr(), readContext.coreStr(), readContext.rightFlankStr(), readContext.toString());
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        configBuilder.addConfigItem(OUTPUT_VCF_FILE, true, "Output filename");
        configBuilder.addPath(INPUT_FILE, true, "Input VCF or TSV");
        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        configBuilder.addPath(REF_GENOME, true, REF_GENOME_CFG_DESC);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        SyntheticVcfWriter syntheticVcfWriter = new SyntheticVcfWriter(configBuilder);
        syntheticVcfWriter.run();
    }
}
