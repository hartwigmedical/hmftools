package com.hartwig.hmftools.geneutils.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache.UNMAPPED_POSITION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.formPairedAltString;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.fromRefAlt;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.sv.VariantAltInsertCoords;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfLiftover
{
    private final String mInputVcf;
    private final String mOutputVcf;
    private final boolean mConvertSVs;
    private final boolean mWriteLiftoverTag;
    private final RefGenomeVersion mSourceVersion;
    private final RefGenomeVersion mDestVersion;

    private final GenomeLiftoverCache mLiftoverCache;

    private List<SAMSequenceRecord> mNewSeqDictionary;
    private final Map<String,List<VariantContext>> mChrVariantMap;

    private static final String INPUT_VCF = "input_vcf";
    private static final String OUTPUT_VCF = "output_vcf";
    private static final String CONVERT_SVS = "convert_svs";
    private static final String SOURCE_REF_GENOME_VERSION = "source_ref_genome_version";
    private static final String WRITE_LIFTOVER_TAG = "write_liftover_tag";

    private static final String LIFTOVER_TAG = "LIFTOVER";
    private static final String LIFTOVER_TAG_DESC = "Genomic coordination liftover information";

    public VcfLiftover(final ConfigBuilder configBuilder)
    {
        mInputVcf = configBuilder.getValue(INPUT_VCF);
        mOutputVcf = configBuilder.getValue(OUTPUT_VCF);
        mConvertSVs = configBuilder.hasFlag(CONVERT_SVS);
        mWriteLiftoverTag = configBuilder.hasFlag(WRITE_LIFTOVER_TAG);

        mSourceVersion = RefGenomeVersion.from(configBuilder.getValue(SOURCE_REF_GENOME_VERSION));
        mDestVersion = mSourceVersion.is37() ? V38 : V37;

        mLiftoverCache = new GenomeLiftoverCache(true);

        mChrVariantMap = Maps.newHashMap();
        mNewSeqDictionary = Lists.newArrayList();
    }

    public void run()
    {
        GU_LOGGER.info("lifting over VCF({}) to ref genome version({})", mInputVcf, mDestVersion);

        long startTimeMs = System.currentTimeMillis();

        try
        {
            VcfFileReader vcfReader = new VcfFileReader(mInputVcf);
            VCFHeader vcfHeader = vcfReader.vcfHeader();

            VCFHeader newHeader = new VCFHeader(vcfHeader.getMetaDataInInputOrder(), vcfHeader.getGenotypeSamples());

            newHeader.addMetaDataLine(new VCFHeaderLine("liftover", format("%s_to_%s", mSourceVersion.toString(), mDestVersion.toString())));

            if(mWriteLiftoverTag)
            {
                newHeader.addMetaDataLine(new VCFInfoHeaderLine(LIFTOVER_TAG, 1, VCFHeaderLineType.String, LIFTOVER_TAG_DESC));
            }

            RefGenomeCoordinates newCoordinates = mDestVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

            // build new dictionary
            for(SAMSequenceRecord sequenceRecord : vcfHeader.getSequenceDictionary().getSequences())
            {
                String chromosome = sequenceRecord.getSequenceName();
                String newChromosome = mDestVersion.versionedChromosome(chromosome);

                int newLength = newCoordinates.length(newChromosome);

                if(newLength == 0)
                {
                    newLength = sequenceRecord.getSequenceLength();
                }

                mNewSeqDictionary.add(new SAMSequenceRecord(newChromosome, newLength));
            }

            VariantContextWriter writer = new VariantContextWriterBuilder()
                    .setReferenceDictionary(new SAMSequenceDictionary(mNewSeqDictionary))
                    .setOutputFile(mOutputVcf)
                    .modifyOption(Options.INDEX_ON_THE_FLY, true)
                    .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                    .build();

            writer.writeHeader(newHeader);

            int liftoverCount = 0;
            int failedLiftoverCount = 0;

            String currentChromosome = "";
            List<VariantContext> variants = null;

            for(VariantContext variantContext : vcfReader.iterator())
            {
                VariantContext newContext = convertVariantContext(variantContext);

                if(newContext == null)
                {
                    ++failedLiftoverCount;
                    continue;
                }

                ++liftoverCount;

                if(!currentChromosome.equals(newContext.getContig()))
                {
                    currentChromosome = newContext.getContig();
                    variants = mChrVariantMap.get(currentChromosome);

                    if(variants == null)
                    {
                        variants = Lists.newArrayList();
                        mChrVariantMap.put(currentChromosome, variants);
                    }
                }

                variants.add(newContext);
            }

            // write variants in order based on new contigs and positions
            for(SAMSequenceRecord seqRecord : mNewSeqDictionary)
            {
                List<VariantContext> newVariants = mChrVariantMap.get(seqRecord.getSequenceName());

                if(newVariants == null)
                    continue;

                Collections.sort(newVariants, Comparator.comparingInt(x -> x.getStart()));

                for(VariantContext variantContext : newVariants)
                {
                    writer.add(variantContext);
                }
            }

            mChrVariantMap.clear();

            writer.close();

            if(failedLiftoverCount > 0)
            {
                GU_LOGGER.warn("lifted over {} variants, failed {} entries", liftoverCount, failedLiftoverCount);
            }
            else
            {
                GU_LOGGER.info("lifted over {} variants", liftoverCount);
            }
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to convert VCF: {}", e.toString());
            System.exit(1);
        }

        GU_LOGGER.info("VCF liftover complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private VariantContext convertVariantContext(final VariantContext variantContext)
    {
        String chromosome = variantContext.getContig();
        int origPosition = variantContext.getStart();

        int newPosition = mLiftoverCache.convertPosition(chromosome, origPosition, mDestVersion);

        if(newPosition == UNMAPPED_POSITION)
        {
            GU_LOGGER.debug("skipped writing unmapped location({}:{})", chromosome, origPosition);
            return null;
        }

        List<String> newAlleles = variantContext.getAlleles().stream().map(x -> x.getDisplayString()).collect(Collectors.toList());
        VariantContextBuilder builder = new VariantContextBuilder(variantContext);

        builder.start(newPosition);

        String newChromosome = mDestVersion.versionedChromosome(chromosome);
        builder.chr(newChromosome);

        if(mWriteLiftoverTag)
        {
            builder.attribute(LIFTOVER_TAG, format("OrigPos=%d", origPosition));
        }

        if(mConvertSVs)
        {
            // convert the mate coordinates as well
            boolean isSgl = StructuralVariantFactory.isSingleBreakend(variantContext);

            if(!isSgl)
            {
                String ref = variantContext.getAlleles().get(0).getDisplayString();
                VariantAltInsertCoords altInsertCoords = fromRefAlt(variantContext.getAlleles().get(1).getDisplayString(), ref);

                String newMateChromosome = mDestVersion.versionedChromosome(altInsertCoords.OtherChromsome);
                int newMatePosition = mLiftoverCache.convertPosition(altInsertCoords.OtherChromsome, altInsertCoords.OtherPosition, mDestVersion);

                if(newMatePosition == UNMAPPED_POSITION)
                {
                    GU_LOGGER.debug("skipped writing unmapped mate SV location({}:{})",
                            altInsertCoords.OtherChromsome, altInsertCoords.OtherPosition);
                    return null;
                }

                String newMateStr = formPairedAltString(
                        altInsertCoords.Alt, altInsertCoords.InsertSequence, newMateChromosome, newMatePosition,
                        altInsertCoords.Orient, altInsertCoords.OtherOrient);

                newAlleles.set(1, newMateStr);
            }
        }

        builder.alleles(newAlleles);

        return builder.make();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(INPUT_VCF, true, "Input VCF to be lifted over");
        configBuilder.addConfigItem(OUTPUT_VCF, true, "Output VCF from lift-over");

        configBuilder.addConfigItem(SOURCE_REF_GENOME_VERSION, true, "Source ref genome version");
        configBuilder.addFlag(CONVERT_SVS, "Sort output by chromosome and position");
        configBuilder.addFlag(WRITE_LIFTOVER_TAG, "Write VCF tag with liftover details");

        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        VcfLiftover vcfLiftover = new VcfLiftover(configBuilder);
        vcfLiftover.run();
    }
}
