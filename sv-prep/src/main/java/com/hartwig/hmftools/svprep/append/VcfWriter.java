package com.hartwig.hmftools.svprep.append;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.READ_CONTEXT_QUALITY;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;

public class VcfWriter
{
    private final AppendConfig mConfig;

    public VcfWriter(final AppendConfig config)
    {
        mConfig = config;
    }

    public void writeBreakends(final Map<String,List<BreakendData>> chrBreakendMap)
    {
        VariantContextWriter writer = initialiseVcf();

        if(writer == null)
            return;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            List<BreakendData> breakends = chrBreakendMap.get(chrStr);

            if(breakends == null)
                continue;

            for(BreakendData breakendData : breakends)
            {
                writeBreakend(writer, breakendData);
            }
        }

        writer.close();
    }

    private VariantContextWriter initialiseVcf()
    {
        IndexedFastaSequenceFile refGenome = null;

        try
        {
            refGenome = new IndexedFastaSequenceFile(new File(mConfig.RefGenomeFile));
        }
        catch(Exception e)
        {
            SV_LOGGER.error("failed to open ref genome: {}", e.toString());
            System.exit(1);
        }

        VcfFileReader vcfFileReader = new VcfFileReader(mConfig.InputVcf);

        final SAMSequenceDictionary sequenceDictionary = refGenome.getSequenceDictionary();

         VariantContextWriter writer = new VariantContextWriterBuilder().setOutputFile(mConfig.OutputVcf)
                .modifyOption(Options.INDEX_ON_THE_FLY, true)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .setReferenceDictionary(sequenceDictionary)
                .build();

         /*
        final SAMSequenceDictionary condensedDictionary = new SAMSequenceDictionary();
        for(SAMSequenceRecord sequence : sequenceDictionary.getSequences())
        {
            if(HumanChromosome.contains(sequence.getContig()) || MitochondrialChromosome.contains(sequence.getContig()))
            {
                condensedDictionary.addSequence(sequence);
            }
        }
         */

        VCFHeader newHeader = new VCFHeader(vcfFileReader.vcfHeader());
        newHeader.getGenotypeSamples().add(mConfig.SampleId);

        // header.setSequenceDictionary(condensedDictionary);
        writer.writeHeader(newHeader);

        return writer;
    }

    private static final List<Allele> NO_CALL = Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL);

    private void writeBreakend(final VariantContextWriter writer, final BreakendData breakendData)
    {
        final VariantContextBuilder builder = new VariantContextBuilder(breakendData.variant());
        final List<Genotype> genotypes = Lists.newArrayList(breakendData.variant().getGenotypes());

        GenotypeBuilder gBuilder = new GenotypeBuilder(mConfig.SampleId);

        int depth = breakendData.depth();
        int junctionSupport = breakendData.totalSupport();
        int refSupport = max(depth - junctionSupport, 0);

        gBuilder.DP(depth)
                .AD(new int[] { refSupport, junctionSupport })
                .alleles(NO_CALL);

        genotypes.add(gBuilder.make());

        VariantContext newVariant = builder.genotypes(genotypes).make();

        writer.add(newVariant);
    }
}
