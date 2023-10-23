package com.hartwig.hmftools.svprep.append;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REF_READ_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.SGL_FRAGMENT_COUNT;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.SV_FRAGMENT_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS_DESCRIPTION;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.append.AppendConstants.JUNCTION_FRAGMENTS;
import static com.hartwig.hmftools.svprep.append.AppendConstants.JUNCTION_FRAGMENTS_DESCRIPTION;
import static com.hartwig.hmftools.svprep.append.AppendConstants.JUNCTION_FRAGMENT_TYPE_COUNT;

import java.io.File;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.samtools.UmiReadType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.svprep.reads.ReadType;

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
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

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

        VCFHeader newHeader = new VCFHeader(vcfFileReader.vcfHeader());
        newHeader.getGenotypeSamples().add(mConfig.SampleId);

        newHeader.addMetaDataLine(new VCFFormatHeaderLine(
                UMI_TYPE_COUNTS, UMI_TYPE_COUNT, VCFHeaderLineType.Integer, UMI_TYPE_COUNTS_DESCRIPTION));

        newHeader.addMetaDataLine(new VCFFormatHeaderLine(
                JUNCTION_FRAGMENTS, JUNCTION_FRAGMENT_TYPE_COUNT, VCFHeaderLineType.Integer, JUNCTION_FRAGMENTS_DESCRIPTION));

        writer.writeHeader(newHeader);

        return writer;
    }

    private static final List<Allele> NO_ALLELES = Lists.newArrayList(Allele.NO_CALL);

    private void writeBreakend(final VariantContextWriter writer, final BreakendData breakendData)
    {
        final VariantContextBuilder builder = new VariantContextBuilder(breakendData.variant());
        final List<Genotype> genotypes = Lists.newArrayList(breakendData.variant().getGenotypes());

        GenotypeBuilder gBuilder = new GenotypeBuilder(mConfig.SampleId);

        int depth = breakendData.depth();
        int junctionSupport = breakendData.totalSupport();
        int refSupport = max(depth - junctionSupport, 0);

        // Gridss does not write AD & DP, so instead for now write ref support into REF and junction support into VF and BVF
        // and capture UMI counts if available
        Map<String, Object> attributes = Maps.newHashMap();
        attributes.put(REF_READ_COVERAGE, refSupport);

        if(breakendData.IsSingle)
        {
            attributes.put(SGL_FRAGMENT_COUNT, junctionSupport);
            attributes.put(SV_FRAGMENT_COUNT, 0);
        }
        else
        {
            attributes.put(SGL_FRAGMENT_COUNT, 0);
            attributes.put(SV_FRAGMENT_COUNT, junctionSupport);
        }

        int[] umiTypeCounts = new int[UMI_TYPE_COUNT];

        for(int i = 0; i < breakendData.umiTypeCounts().length; ++i)
        {
            umiTypeCounts[i + 3] = breakendData.umiTypeCounts()[i];
        }

        attributes.put(UMI_TYPE_COUNTS, umiTypeCounts);

        int[] juncFragCounts = new int[JUNCTION_FRAGMENT_TYPE_COUNT];
        juncFragCounts[0] = breakendData.readTypeSupport()[ReadType.JUNCTION.ordinal()];
        juncFragCounts[1] = breakendData.readTypeSupport()[ReadType.EXACT_SUPPORT.ordinal()];
        juncFragCounts[2] = breakendData.readTypeSupport()[ReadType.SUPPORT.ordinal()];
        attributes.put(JUNCTION_FRAGMENTS, juncFragCounts);

        gBuilder.attributes(attributes);
        gBuilder.alleles(NO_ALLELES);

        /*
        gBuilder.DP(depth)
                .AD(new int[] { refSupport, junctionSupport })
                .alleles(NO_CALL);
        */

        genotypes.add(gBuilder.make());

        VariantContext newVariant = builder.genotypes(genotypes).make();

        writer.add(newVariant);
    }
}
